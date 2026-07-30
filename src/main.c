/* ResonarSat command-line driver. */

#include "resonarsat/focus.h"
#include "resonarsat/geom.h"
#include "resonarsat/ccd.h"
#include "resonarsat/microm.h"
#include "resonarsat/raster.h"
#include "resonarsat/simulate.h"
#include "resonarsat/geocode.h"
#include "resonarsat/readers.h"
#include "resonarsat/subaperture.h"
#include "resonarsat/tomo.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Print the top-level usage summary. */
static void rs_usage(void)
{
    printf(
"resonarsat -- SAR Doppler tomography and micro-motion analysis\n"
"\n"
"  feasibility   compute the observable vibration band and what it costs\n"
"  info          print a product's geometry and timing\n"
"  focus         form an image from phase history (backprojection)\n"
"  mmotion       track sub-looks and extract vibration spectra\n"
"  tomo          focus vibration observations into a depth profile\n"
"  sweep         vary the assumed constants and see where depths move\n"
"\n"
"Run 'resonarsat <command>' with no arguments for command-specific help.\n");
}

/* Warn when a grid cell is coarser than the resolution the collect supports.
 *
 * Backprojecting onto an under-sampled grid aliases: each scatterer acquires
 * ghost peaks that look like additional targets. Silence here would mean
 * shipping artefacts that resemble structure, which is the one failure mode
 * this project cannot afford, so the warning is unconditional and names the
 * cell size that would be safe. */
static void rs_warn_sampling(const rs_cphd_t *cphd, size_t pulse_count, double cell)
{
    const double res = rs_focus_azimuth_resolution(cphd, pulse_count);
    if (res > 0.0 && cell > res) {
        fprintf(stderr,
            "warning: grid cell %.3f m is coarser than the %.3f m azimuth\n"
            "         resolution this collect supports. The image will alias:\n"
            "         each scatterer gains ghost peaks that resemble real\n"
            "         targets. Use --cell %.3f or finer, or treat the result\n"
            "         as deliberately multi-looked.\n",
            cell, res, res);
    }
}

/* Look up a --name value pair in an argument vector, returning the value or
 * NULL. Keeps the subcommand parsers short; there are few enough options that a
 * real parser would be more machinery than the problem needs. */
static const char *rs_opt(int argc, char **argv, const char *name)
{
    for (int i = 0; i < argc - 1; i++) {
        if (strcmp(argv[i], name) == 0) return argv[i + 1];
    }
    return NULL;
}

/* True when a valueless flag such as --no-detrend is present.
 *
 * Separate from rs_opt() because that one scans to argc-1, looking for a value
 * after the name; a flag written last would be missed. */
static int rs_opt_flag(int argc, char **argv, const char *name)
{
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], name) == 0) return 1;
    }
    return 0;
}

/* As rs_opt(), but parsing the value as a double and falling back to a default
 * when the option is absent. */
static double rs_opt_double(int argc, char **argv, const char *name, double fallback)
{
    const char *v = rs_opt(argc, argv, name);
    return v ? atof(v) : fallback;
}

/* Resolve --no-optimize and say on stderr exactly what it did and did not do.
 *
 * One flag rather than three because the audit is only meaningful if the whole
 * chain is in the same mode, and a caller who has to remember three names will
 * one day set two of them. The library keeps the settings separate per stage --
 * rs_focus_opts_t, rs_subap_params_t, rs_microm_params_t -- so nothing here
 * depends on global mutable state and each stage still documents its own
 * behaviour; this function is the single place that ties them together.
 *
 * The notice is printed, and printed in this much detail, because the flag's name
 * invites a claim it does not support. "Unoptimized baseline" suggests the
 * numbers that come out are more trustworthy, and for the backprojection half
 * that is not merely unsupported but false: the output is bitwise identical, so a
 * run made with the flag is the same measurement, not a better one. Only the
 * correlation peak search can move a number. Saying so at the point of use costs
 * four lines and keeps a comparison from being read as a confirmation. */
static int rs_opt_no_optimize(int argc, char **argv)
{
    if (!rs_opt_flag(argc, argv, "--no-optimize")) return 0;

    fprintf(stderr,
        "--no-optimize: unoptimised reference mode.\n"
        "  correlation peak search  WHOLE zero-padded surface, global maximum,\n"
        "                           instead of one pixel about the integer peak.\n"
        "                           This is the only change that can move a number.\n"
        "  backprojection           single-threaded, ascending cell order. The\n"
        "                           samples are BITWISE IDENTICAL to a threaded run:\n"
        "                           each cell accumulates privately over pulses in\n"
        "                           chronological order and no accumulator is shared,\n"
        "                           so there is no threading drift here to remove.\n"
        "                           This half is a reproducibility check, not a fix.\n"
        "  tracking loop            serial, so window visitation order is fixed.\n"
        "                           Also cannot change a result; windows are\n"
        "                           independent.\n"
        "  Measured cost: the exhaustive search is 1.7-3.2x the optimised one per\n"
        "  call (2.5x at these defaults), and running serially costs about 4x more\n"
        "  on eight cores. Single digits overall, not orders of magnitude -- this\n"
        "  is affordable on a full-scale scene.\n"
        "  A result from this mode is a second measurement by a slower route, to be\n"
        "  compared with the first -- not a more trustworthy one. Neither passes a\n"
        "  null test alone.\n");
    return 1;
}

/* Report the vibration band an acquisition can observe, and what each choice
 * costs in azimuth resolution.
 *
 * This is the cheapest checkpoint in the whole project and belongs before any
 * large download: it answers whether the structure of interest is observable at
 * all, from geometry alone, in about a millisecond. */
static int rs_cmd_feasibility(int argc, char **argv)
{
    if (argc < 2) {
        printf("usage: resonarsat feasibility --dwell S --range M --velocity MS "
               "--wavelength M [--n N] [--tomo-velocity MS --tomo-frequency HZ]\n");
        return 1;
    }

    const double dwell  = rs_opt_double(argc, argv, "--dwell", 10.0);
    const double range  = rs_opt_double(argc, argv, "--range", 500000.0);
    const double vp     = rs_opt_double(argc, argv, "--velocity", 7500.0);
    const double lambda = rs_opt_double(argc, argv, "--wavelength", 0.031);

    if (dwell <= 0.0 || range <= 0.0 || vp <= 0.0 || lambda <= 0.0) {
        fprintf(stderr, "feasibility: dwell, range, velocity and wavelength "
                        "must all be positive\n");
        return 1;
    }

    printf("Acquisition\n");
    printf("  dwell              %.3f s\n", dwell);
    printf("  slant range        %.1f km\n", range / 1000.0);
    printf("  platform speed     %.1f m/s\n", vp);
    printf("  wavelength         %.4f m\n", lambda);
    printf("  full aperture      %.1f km\n", rs_orbital_aperture(vp, dwell) / 1000.0);
    printf("  full-aperture res  %.2f m\n",
           rs_azimuth_resolution(lambda, range, vp, dwell));

    printf("\nVibration band against sub-look resolution\n");
    printf("  %8s %10s %10s %12s\n", "N_looks", "dt (s)", "f_max (Hz)", "d_az (m)");
    const int counts[] = { 8, 16, 32, 64, 128, 300 };
    for (size_t i = 0; i < sizeof counts / sizeof counts[0]; i++) {
        const int n = counts[i];
        const double dt = dwell / (double)n;
        printf("  %8d %10.4f %10.2f %12.2f\n",
               n, dt, rs_vibration_fmax(n, dwell),
               rs_azimuth_resolution(lambda, range, vp, dt));
    }
    printf("\n  Resolution degrades as the square of the frequency reach.\n"
           "  Independently validated results for this technique sit at 1-4 Hz.\n");

    /* The depth chain is only meaningful once the operator supplies the two
     * assumptions, so it is printed only when both are given. */
    const double tv = rs_opt_double(argc, argv, "--tomo-velocity", 0.0);
    const double tf = rs_opt_double(argc, argv, "--tomo-frequency", 0.0);
    if (tv > 0.0 && tf > 0.0) {
        const double lam_ac = rs_acoustic_wavelength(tv, tf, RS_WAVELEN_PAPER);
        const double A = rs_orbital_aperture(vp, dwell);
        printf("\nDepth chain (ASSUMED constants, not measured)\n");
        printf("  assumed wave speed  %.1f m/s\n", tv);
        printf("  assumed frequency   %.1f Hz\n", tf);
        printf("  acoustic wavelength %.4f m   (paper convention, v/2f)\n", lam_ac);
        printf("  ... patent          %.4f m   (v/f -- differs by 2x)\n",
               rs_acoustic_wavelength(tv, tf, RS_WAVELEN_PATENT));
        printf("  depth resolution    %.3f m\n", rs_tomo_resolution(lam_ac, range, A));
        printf("\n  These depths are scaled by the two assumptions above. Nothing in\n"
               "  the pipeline measures either of them.\n");
    }

    return 0;
}

/* Write one horizontal slice of a tomogram as an image.
 *
 * This is the figure format the published work presents: both display axes are
 * image coordinates -- azimuth and range -- and the depth is not an axis at all
 * but a label attached to whichever slice you chose to look at. Understanding
 * that is most of understanding what a tomogram of this kind shows.
 *
 * It follows that every shape in such a picture was inherited from the surface
 * image geometry. Nothing about an outline was measured from below; the depth
 * stage only decides which sheet of the stack is displayed, and the metres
 * printed on that sheet come from the two assumed constants.
 *
 * Which makes this a demonstration as well as an export. Render the same
 * normalised depth under two different assumed wavelengths and the images come
 * out identical to the byte, while the captions differ by the ratio of the
 * assumptions. That is the sweep's slope-of-one result expressed as a picture,
 * and it can be checked with cmp(1) rather than a regression.
 *
 * 'depth_m' selects the nearest available slice; the depth actually used is
 * written back through 'used' and its index through 'index' when non-NULL. */
static resonarsat_status_t rs_tomo_write_slice(const rs_tomo_t *t, double depth_m,
                                               const char *path, double *used,
                                               size_t *index,
                                              rs_palette_t palette)
{
    if (!t || !t->profile || !t->depth || !path || t->n_depth == 0) return RS_ERR_ARG;

    size_t k = 0;
    double best = fabs(t->depth[0] - depth_m);
    for (size_t i = 1; i < t->n_depth; i++) {
        const double d = fabs(t->depth[i] - depth_m);
        if (d < best) { best = d; k = i; }
    }

    double *plane = malloc(t->n_win * sizeof *plane);
    if (!plane) return RS_ERR_ALLOC;
    for (size_t w = 0; w < t->n_win; w++) plane[w] = t->profile[w * t->n_depth + k];

    const resonarsat_status_t st =
        rs_raster_write_map(plane, t->n_win_az, t->n_win_rg, path, 0.0, 0.0, palette);
    free(plane);

    if (used)  *used = t->depth[k];
    if (index) *index = k;
    return st;
}

/* Write a vertical section of a tomogram along a line between two windows.
 *
 * This is the other figure the published work presents, and the one its claims
 * actually rest on: an operator draws a line across the SAR image, and the
 * section beneath that line is rendered with depth down the vertical axis and
 * distance along the line across the horizontal.
 *
 * Two properties of that presentation are worth being able to reproduce, because
 * they are easier to demonstrate than to argue.
 *
 * The horizontal axis is a count of samples along the line, not a distance, so a
 * section can be stretched or squeezed horizontally at will. Shape and aspect
 * ratio are therefore not readable from such a figure: whether a feature looks
 * like a narrow vertical shaft or a broad blob depends on a display choice, and
 * a claim about a shaft is a claim about shape. 'n_samp' is exposed for exactly
 * this reason -- render the same section at two widths and the geometry of every
 * feature changes while the data does not.
 *
 * The vertical axis is oversampled whenever the depth grid is finer than the
 * resolution the geometry supports. A target that is not resolved in depth
 * smears along that axis into a vertical streak, which is what an unresolved
 * point looks like and also what a shaft looks like. Since rs_tomo_params_check()
 * refuses grids finer than the resolution, a section rendered through this
 * pipeline cannot manufacture that appearance by interpolation -- which makes it
 * a fair control for a figure that can.
 *
 * The line runs from (az0,rg0) to (az1,rg1) in window indices. Samples are taken
 * at nearest-neighbour window positions; interpolating between windows would
 * invent smoothness the measurement does not have. */
static resonarsat_status_t rs_tomo_write_section(const rs_tomo_t *t,
                                                 double az0, double rg0,
                                                 double az1, double rg1,
                                                 size_t n_samp, const char *path,
                                                rs_palette_t palette)
{
    if (!t || !t->profile || !path || t->n_depth == 0 || n_samp < 2) return RS_ERR_ARG;
    if (t->n_win_az == 0 || t->n_win_rg == 0) return RS_ERR_ARG;

    /* Rows are depth, columns are position along the line, so the image comes
     * out oriented the way the published sections are. */
    double *img = malloc(t->n_depth * n_samp * sizeof *img);
    if (!img) return RS_ERR_ALLOC;

    for (size_t i = 0; i < n_samp; i++) {
        const double f = (double)i / (double)(n_samp - 1);
        long a = lround(az0 + f * (az1 - az0));
        long r = lround(rg0 + f * (rg1 - rg0));
        if (a < 0) a = 0;
        if (r < 0) r = 0;
        if ((size_t)a >= t->n_win_az) a = (long)t->n_win_az - 1;
        if ((size_t)r >= t->n_win_rg) r = (long)t->n_win_rg - 1;

        const size_t w = (size_t)a * t->n_win_rg + (size_t)r;
        for (size_t k = 0; k < t->n_depth; k++) {
            img[k * n_samp + i] = t->profile[w * t->n_depth + k];
        }
    }

    const resonarsat_status_t st =
        rs_raster_write_map(img, t->n_depth, n_samp, path, 0.0, 0.0, palette);
    free(img);
    return st;
}

/* Load phase history from either the project's interchange format or a real
 * CPHD 1.x product, chosen by inspecting the file rather than by a flag.
 *
 * The two are told apart by their first bytes: a CPHD product begins with the
 * ASCII text "CPHD/", the interchange format with a binary magic number. Every
 * command that takes phase history goes through here, so a real collect works
 * anywhere a simulated one does and the commands themselves stay unaware of the
 * distinction.
 *
 * Real products are read with a range window, because they are not sized for
 * memory: the reference collect is 3.8 GB of signal over a 7 km swath, of which
 * the sub-aperture stage uses a few hundred metres. See rs_cphd_read_opts_t.
 *
 * 'pulse_stride' keeps every nth pulse. It was for a time accepted on the
 * command line, documented in the out-of-memory message below, and then
 * discarded here in favour of a hard-coded 0 -- so runs that passed it read the
 * whole collect anyway, at full azimuth resolution, and aliased when that
 * resolution was rendered onto a coarse grid. Striding is a real cost and is
 * reported rather than applied quietly: it lowers the effective PRF, and with it
 * the vibration frequency the sub-aperture stage can reach without aliasing. */
static resonarsat_status_t rs_load_cphd(rs_cphd_t *c, const char *path, size_t rbin_window,
                                        size_t pulse_stride, size_t max_pulses)
{
    char probe[5] = { 0 };
    FILE *f = fopen(path, "rb");
    if (!f) {
        rs_set_error("cannot open %s", path);
        return RS_ERR_IO;
    }
    const size_t n = fread(probe, 1, sizeof probe, f);
    fclose(f);

    if (n == sizeof probe && memcmp(probe, "CPHD/", 5) == 0) {
        const rs_cphd_read_opts_t o = { .rbin_window = rbin_window,
                                        .pulse_stride = pulse_stride,
                                        .max_pulses = max_pulses };
        const resonarsat_status_t st = rs_read_cphd(path, &o, c);
        if (st == RS_OK && pulse_stride > 1) {
            fprintf(stderr,
                "warning: --pulse-stride %zu keeps every %zuth pulse, lowering the\n"
                "         effective PRF to %.2f Hz and the observable vibration\n"
                "         band with it. Acceptable for a registration image;\n"
                "         not for a measurement run.\n",
                pulse_stride, pulse_stride, c->prf);

            /* And it bounds how far from the motion-compensation reference the
             * grid may extend, which is the part that bites silently.
             *
             * Azimuth sampling at PRF is unambiguous only over
             * lambda*R*PRF/(2V) of ground. Beyond that the spectrum folds and
             * the image does not degrade gracefully -- it fills with aliased
             * energy that looks exactly like speckle, so a grid placed past the
             * limit comes back plausible and empty. A whole-scene registration
             * sweep is precisely where this happens, because the natural
             * instinct on a slow backprojection is to raise the stride and
             * widen the grid at the same time, and the two limits move in
             * opposite directions. The figure is printed here so the caller can
             * compare it against the grid they are about to ask for. */
            const double R = (c->r_near > 0.0) ? c->r_near : 0.0;
            /* Mean speed straight off the recorded track, rather than an
             * assumed orbital value: the whole point of the figure is that it
             * describes THIS collect. */
            double V = 0.0;
            if (c->n_pulse > 1 && c->pos && c->t) {
                const size_t last = c->n_pulse - 1;
                const double dt = c->t[last] - c->t[0];
                const double dx = c->pos[3 * last + 0] - c->pos[0];
                const double dy = c->pos[3 * last + 1] - c->pos[1];
                const double dz = c->pos[3 * last + 2] - c->pos[2];
                if (dt > 0.0) V = sqrt(dx * dx + dy * dy + dz * dz) / dt;
            }
            if (c->lambda > 0.0 && R > 0.0 && V > 0.0 && c->prf > 0.0) {
                const double extent = c->lambda * R * c->prf / (2.0 * V);
                fprintf(stderr,
                    "         At this PRF the unambiguous azimuth extent is\n"
                    "         %.0f m (+/-%.0f m from the scene reference). A grid\n"
                    "         wider than that aliases into something that still\n"
                    "         looks like speckle.\n", extent, 0.5 * extent);
            }
        }
        return st;
    }
    if (pulse_stride > 1) {
        fprintf(stderr, "warning: --pulse-stride applies to CPHD products only; "
                        "ignored for this input\n");
    }
    return rs_cphd_read(c, path);
}

/* Parse --reference and --b-shift, which together select how master and slave
 * sub-bands are paired.
 *
 * These are read before the stack is built, because 'pair' decides whether the
 * slave bands are synthesised at all. Returns the reference mode; on an unknown
 * name it reports and returns -1.
 *
 * Shared by mmotion and tomo so the two cannot drift apart on a choice that
 * determines what the displacement series is. */
static int rs_parse_reference(int argc, char **argv, rs_subap_params_t *sp)
{
    sp->b_shift_hz = rs_opt_double(argc, argv, "--b-shift", 0.0);

    int ref = RS_MICROM_REF_FIRST;
    const char *name = rs_opt(argc, argv, "--reference");
    if (name) {
        if      (strcmp(name, "first") == 0)    ref = RS_MICROM_REF_FIRST;
        else if (strcmp(name, "adjacent") == 0) ref = RS_MICROM_REF_ADJACENT;
        else if (strcmp(name, "pair") == 0)     ref = RS_MICROM_REF_PAIR;
        else {
            fprintf(stderr, "unknown --reference '%s'; expected first, "
                            "adjacent or pair\n", name);
            return -1;
        }
    }
    sp->pair = (ref == RS_MICROM_REF_PAIR);

    if (sp->b_shift_hz > 0.0 && !sp->pair) {
        fprintf(stderr,
            "warning: --b-shift %g Hz sets the master-slave separation, which\n"
            "         only has an effect under --reference pair. Ignored here.\n",
            sp->b_shift_hz);
    }
    return ref;
}

/* Fill the tomographic geometry from the collect that was actually loaded.
 *
 * rs_tomo_params_default() carries the simulator's geometry -- 500 km slant,
 * 75 km aperture, 35 degrees -- and before this existed those defaults survived
 * into every run, including runs on real products whose geometry is nothing
 * like them. On simulated input the defaults happened to match the simulator, so
 * nothing disagreed and the substitution stayed invisible. The consequence on a
 * real collect is that the depth resolution, the unambiguous extent and every
 * absolute depth are computed for a scene that was not observed.
 *
 * All three quantities are already implied by the phase history. The slant range
 * is the reference range at the middle of the dwell. The aperture is the
 * straight-line distance the phase centre travelled between the first and last
 * pulse. The incidence angle follows from the platform position expressed in the
 * scene frame, whose z axis is the reference surface normal: the angle between
 * the line of sight and that normal.
 *
 * Explicit --range, --aperture and --incidence still win, so a caller can force
 * a geometry deliberately, but they are no longer required to avoid a wrong one
 * silently. */
static void rs_tomo_geometry_from_cphd(rs_tomo_params_t *tp, const rs_cphd_t *c)
{
    if (c && c->n_pulse > 1 && c->pos && c->r_ref) {
        const size_t mid = c->n_pulse / 2;
        if (c->r_ref[mid] > 0.0) tp->slant_range = c->r_ref[mid];

        const double *a = c->pos, *b = c->pos + 3 * (c->n_pulse - 1);
        const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
        const double arc = sqrt(dx * dx + dy * dy + dz * dz);
        if (arc > 0.0) tp->aperture = arc;

        const double *m = c->pos + 3 * mid;
        const double rm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
        if (rm > 0.0) {
            double ct = m[2] / rm;
            if (ct < -1.0) ct = -1.0;
            if (ct > 1.0) ct = 1.0;
            tp->incidence = acos(fabs(ct));
        }
    }

}

/* Apply explicit geometry overrides from the command line, after
 * rs_tomo_geometry_from_cphd() has taken what the collect implies. */
static void rs_tomo_geometry_overrides(rs_tomo_params_t *tp, int argc, char **argv)
{
    tp->slant_range = rs_opt_double(argc, argv, "--range", tp->slant_range);
    tp->aperture    = rs_opt_double(argc, argv, "--aperture", tp->aperture);
    const double inc_deg = rs_opt_double(argc, argv, "--incidence", -1.0);
    if (inc_deg > 0.0) tp->incidence = inc_deg * M_PI / 180.0;
}

/* Print a focused product's geometry and timing, and flag anything implausible.
 *
 * Azimuth sampling rate and transmit PRF are printed as separate labelled
 * fields, deliberately: conflating them is the specific error the data model
 * exists to prevent, and showing both makes a substitution visible. */
static int rs_cmd_info(int argc, char **argv)
{
    if (argc < 2) {
        printf("usage: resonarsat info --uavsar SLC --ann ANN\n"
               "       resonarsat info --cphd FILE\n"
               "       resonarsat info --sicd FILE.nitf\n");
        return 1;
    }

    const char *sicd_path = rs_opt(argc, argv, "--sicd");
    if (sicd_path) {
        rs_slc_t img;
        const resonarsat_status_t sst = rs_read_sicd(sicd_path, &img);
        if (sst != RS_OK) { rs_report_error("info", sst); return 1; }
        printf("source              %s (focused image)\n", img.source);
        printf("dimensions          %zu azimuth x %zu range\n", img.n_az, img.n_rg);
        printf("carrier             %.4f GHz\n", img.fc / 1e9);
        printf("wavelength          %.4f m\n", img.lambda);
        printf("azimuth spacing     %.4f m\n", img.az_spacing_m);
        printf("range spacing       %.4f m\n", img.rg_spacing_m);
        printf("slant range         %.1f km\n", img.r0 / 1000.0);
        printf("incidence           %.2f deg\n", img.incidence * 180.0 / M_PI);
        printf("platform speed      %.1f m/s\n", img.v_platform);
        printf("collect duration    %.3f s\n", img.t_dwell);
        const resonarsat_status_t vst = rs_slc_validate(&img);
        if (vst != RS_OK) {
            fprintf(stderr, "\nwarning: %s\n", rs_last_error());
        }
        if (img.plane.valid) {
            double la, lo, he;
            if (rs_geo_plane_llh(&img.plane, 0.0, 0.0, &la, &lo, &he) == RS_OK) {
                printf("\nscene corners (WGS 84)%s:\n",
                       img.plane.is_slant ? ", by range-Doppler" : "");
                const struct { const char *name; double a, r; } pts[] = {
                    { "first az, first rg", 0.0, 0.0 },
                    { "last az,  first rg", (double)(img.n_az - 1), 0.0 },
                    { "first az, last rg ", 0.0, (double)(img.n_rg - 1) },
                    { "last az,  last rg ", (double)(img.n_az - 1), (double)(img.n_rg - 1) },
                };
                for (size_t k = 0; k < 4; k++) {
                    double ecf[3];
                    if (rs_geo_plane_point(&img.plane,
                                           pts[k].a * img.az_spacing_m,
                                           pts[k].r * img.rg_spacing_m, ecf) != RS_OK) {
                        continue;
                    }

                    /* A slant-plane grid carries the right range and Doppler but
                     * the wrong height; range-Doppler puts it on the ground. */
                    if (img.plane.is_slant) {
                        double g[3];
                        if (rs_geo_slant_to_ground(img.plane.sensor, img.plane.sensor_vel,
                                                   ecf, img.plane.ref_hae, g) == RS_OK) {
                            memcpy(ecf, g, sizeof ecf);
                        }
                    }
                    if (rs_geo_ecf_to_llh(ecf, &la, &lo, &he) == RS_OK) {
                        printf("  %s  %10.6f, %11.6f   %7.1f m\n",
                               pts[k].name, la, lo, he);
                    }
                }
                double cecf[3];
                if (rs_geo_plane_point(&img.plane,
                                       0.5 * (double)img.n_az * img.az_spacing_m,
                                       0.5 * (double)img.n_rg * img.rg_spacing_m,
                                       cecf) == RS_OK) {
                    if (img.plane.is_slant) {
                        double g[3];
                        if (rs_geo_slant_to_ground(img.plane.sensor, img.plane.sensor_vel,
                                                   cecf, img.plane.ref_hae, g) == RS_OK) {
                            memcpy(cecf, g, sizeof cecf);
                        }
                    }
                    if (rs_geo_ecf_to_llh(cecf, &la, &lo, &he) == RS_OK) {
                        printf("  centre              %10.6f, %11.6f   %7.1f m\n", la, lo, he);
                    }
                }
            }
        }

        printf("\nThis product is already focused. Sub-apertures come from spectral\n"
               "splitting (rs_subaperture_split), which is the route the published\n"
               "method describes: transform, band-pass, inverse transform, track.\n");
        rs_slc_free(&img);
        return 0;
    }

    const char *cphd_path = rs_opt(argc, argv, "--cphd");
    if (cphd_path) {
        rs_cphd_t c;
        const size_t rbin_window = (size_t)rs_opt_double(argc, argv, "--rbins", 0.0);
        resonarsat_status_t st = rs_load_cphd(&c, cphd_path, rbin_window,
                       (size_t)rs_opt_double(argc, argv, "--pulse-stride", 0.0),
                       (size_t)rs_opt_double(argc, argv, "--max-pulses", 0.0));
        if (st != RS_OK) { rs_report_error("info", st); return 1; }

        printf("source              phase history (unfocused)\n");
        printf("pulses              %zu\n", c.n_pulse);
        printf("range bins          %zu\n", c.n_rbin);
        printf("carrier             %.4f GHz\n", c.fc / 1e9);
        printf("wavelength          %.4f m\n", c.lambda);
        printf("transmit PRF        %.2f Hz    (NOT the vibration sampling rate)\n", c.prf);
        printf("dwell               %.3f s\n", c.t[c.n_pulse - 1] - c.t[0]);
        printf("near range          %.1f m\n", c.r_near);
        printf("range bin spacing   %.4f m\n", c.dr);
        printf("\nThis product carries no image. Run 'resonarsat focus' first;\n"
               "quicklook and the processing stages need a focused product.\n");
        rs_cphd_free(&c);
        return 0;
    }

    const char *slc = rs_opt(argc, argv, "--uavsar");
    const char *ann = rs_opt(argc, argv, "--ann");
    if (!slc || !ann) {
        fprintf(stderr, "info: --uavsar requires --ann\n");
        return 1;
    }

    rs_slc_t img;
    resonarsat_status_t st = rs_read_uavsar(slc, ann, &img);
    if (st != RS_OK) { rs_report_error("info", st); return 1; }

    printf("source              %s\n", img.source);
    printf("dimensions          %zu azimuth x %zu range\n", img.n_az, img.n_rg);
    printf("carrier             %.4f GHz\n", img.fc / 1e9);
    printf("wavelength          %.4f m\n", img.lambda);
    printf("azimuth line rate   %.3f Hz    (1/azimuth_time_interval)\n", img.fs_az);
    printf("transmit PRF        %.3f Hz    (diagnostic only)\n", img.pulse_prf);
    printf("azimuth spacing     %.4f m\n", img.az_spacing_m);
    printf("range spacing       %.4f m\n", img.rg_spacing_m);
    printf("platform speed      %.2f m/s\n", img.v_platform);
    printf("dwell               %.3f s\n", img.t_dwell);

    double fdc = 0.0, bw = 0.0;
    if (rs_estimate_doppler_centroid(&img, &fdc) == RS_OK)
        printf("Doppler centroid    %.3f Hz  (measured from data)\n", fdc);
    if (rs_estimate_doppler_bandwidth(&img, 0.5, &bw) == RS_OK)
        printf("Doppler bandwidth   %.3f Hz  (-3 dB, measured)\n", bw);

    rs_slc_free(&img);
    return 0;
}

/* Build a sub-look stack by the route the caller asked for.
 *
 * WHY THIS EXISTS. There are two ways to divide an aperture, they are not
 * interchangeable, and for most of this project's life only one of them was
 * reachable.
 *
 *   "pulse"    Focus n_looks shifted windows of the phase history. The aperture
 *              is divided in TIME, so every timing quantity comes off the
 *              recorded pulse times and none is inferred. Exact by construction.
 *
 *   "uniform"  Focus the full aperture once, then split the focused image's
 *              Doppler spectrum into equal sub-bands.
 *
 *   "paper"    Focus the full aperture once, then apply the decomposition of
 *              Biondi & Malanga (2022) section 3.1: hold out a fraction B_DL of
 *              the Doppler band and step the remainder across the spectrum in
 *              rigid shifts. The paper credits the held-out band with providing
 *              "a sufficient sensitivity to estimate target motions", so it is
 *              not a detail -- it is the stage the method's motion sensitivity
 *              is attributed to.
 *
 * The default remains "pulse". It is the better measurement and it is what the
 * end-to-end test exercises. But defaulting to it while describing this project
 * as a reimplementation of the paper overstated the correspondence: every real
 * result here was produced with a uniform filter bank the paper does not use.
 *
 * A CAVEAT THAT APPLIES TO BOTH SPECTRAL ROUTES. They need the focused image's
 * azimuth sampling rate to map bands to times, and for a BACKPROJECTED image
 * that quantity is nominal: rs_focus_backproject() derives it from the grid as
 * t_dwell/(n_x-1), because a focused image's azimuth axis is spatial rather
 * than temporal. Relative band positions are therefore trustworthy and the
 * absolute time step between sub-looks is not, which makes the recovered
 * vibration FREQUENCY axis uncertain by whatever that mapping is wrong by.
 * Coherence and detection strength do not depend on it; reported frequencies
 * do. The paper starts from a vendor-focused SLC carrying real azimuth timing
 * and does not face this. */
static resonarsat_status_t rs_build_subaps(const rs_cphd_t *c, const rs_grid_t *grid,
                                           rs_subap_params_t *sp, const char *route,
                                           rs_subap_stack_t *stack)
{
    if (!route || strcmp(route, "pulse") == 0) {
        sp->mode = RS_SUBAP_UNIFORM;
        return rs_subaperture_from_cphd(c, grid, sp, stack);
    }

    if (strcmp(route, "paper") == 0) {
        sp->mode = RS_SUBAP_PAPER;
    } else if (strcmp(route, "uniform") == 0) {
        sp->mode = RS_SUBAP_UNIFORM;
    } else {
        rs_set_error("subap route '%s' is not pulse, uniform or paper", route);
        return RS_ERR_ARG;
    }

    /* The spectral routes need a focused image first. Splitting one requires at
     * least twice as many azimuth lines as looks, which is checked there and
     * reported rather than silently reducing the look count. */
    rs_slc_t full;
    resonarsat_status_t st = rs_slc_alloc(&full, grid->n_x, grid->n_y);
    if (st != RS_OK) return st;

    fprintf(stderr, "subap route '%s': focusing the full aperture first ...\n", route);
    {
        const rs_focus_opts_t fopts = { .single_thread = sp->single_thread,
                                        .range_taps = sp->range_taps };
        st = rs_focus_backproject_opts(c, grid, 0, c->n_pulse, &fopts, &full);
    }
    if (st != RS_OK) { rs_slc_free(&full); return st; }

    st = rs_subaperture_split(&full, sp, stack);
    rs_slc_free(&full);
    return st;
}

/* Parse "--offset X,Y" into a grid origin, in metres from the scene centre.
 *
 * WHY THIS EXISTS. The processing grid is smaller than the collect: a 2048-cell
 * grid at 1 m covers 2.0 km of a 4 km scene, and hard-coding the origin at zero
 * silently restricts every command to the middle quarter. That is invisible in
 * the output -- the image looks complete, because it IS a complete image of
 * somewhere -- and it is the failure that makes a target near the scene edge
 * come back empty while appearing to have been processed.
 *
 * Structures worth pointing this software at are rarely at the scene centre,
 * because the scene was framed by whoever tasked the collect and not by whoever
 * is analysing it. So the origin is an input.
 *
 * 'x' runs along the platform track and 'y' across it, matching rs_grid_t; the
 * pair is therefore in the scene's own planar frame rather than in north/east,
 * and a caller converting from geodetic coordinates has to go through the image
 * plane. A NULL or malformed specification leaves the origin at the scene
 * centre, which is the historical behaviour. */
static void rs_parse_offset(const char *spec, double origin[3])
{
    origin[0] = origin[1] = origin[2] = 0.0;
    if (!spec || !*spec) return;

    char buf[128];
    snprintf(buf, sizeof buf, "%s", spec);
    char *comma = strchr(buf, ',');
    if (!comma) {
        fprintf(stderr, "warning: --offset needs X,Y in metres; ignoring '%s'\n", spec);
        return;
    }
    *comma = '\0';
    origin[0] = atof(buf);
    origin[1] = atof(comma + 1);
}

/* Resolve "--at LAT,LON" against the product's own image plane, overriding
 * --offset when present.
 *
 * WHY THE PLAIN OFFSET IS NOT ENOUGH. The frame --offset speaks is the file's:
 * CPHD's uIAX/uIAY, SICD's row and column vectors. Those axes point where the
 * collector chose, so converting a place into them by assuming north/east, or
 * azimuth/ground-range, produces a number that is wrong in a way no output
 * reveals -- the image is a complete, well-focused image of the wrong ground.
 * That happened on this project's own Giza runs, whose manifests recorded a
 * hand-derived Khufu offset roughly 900 m from the pyramid; the geodetic values
 * in data/README.md were right and were overridden as though measured.
 *
 * Giving the coordinates instead makes the product do the conversion, which it
 * can do exactly. Returns non-zero on a malformed specification or a product
 * that carries no plane, so the caller can refuse rather than silently fall
 * back to a grid centred somewhere else. */
static int rs_resolve_at(const char *spec, const rs_geo_plane_t *plane,
                         double origin[3])
{
    if (!spec || !*spec) return 0;

    char buf[128];
    snprintf(buf, sizeof buf, "%s", spec);
    char *comma = strchr(buf, ',');
    if (!comma) {
        fprintf(stderr, "--at needs LAT,LON in degrees; got '%s'\n", spec);
        return 1;
    }
    *comma = '\0';
    const double lat = atof(buf), lon = atof(comma + 1);

    if (!plane || !plane->valid) {
        fprintf(stderr, "--at needs the product's image plane, which this input "
                        "does not carry; use --offset X,Y instead\n");
        return 1;
    }

    double x = 0.0, y = 0.0;
    const resonarsat_status_t st =
        rs_geo_plane_offset(plane, lat, lon, plane->ref_hae, &x, &y);
    if (st != RS_OK) {
        rs_report_error("--at", st);
        return 1;
    }
    origin[0] = x;
    origin[1] = y;
    printf("--at %.6f,%.6f resolves to --offset %.0f,%.0f in this product's "
           "image plane\n", lat, lon, x, y);
    return 0;
}

/* Form an image from phase history. */
static int rs_cmd_focus(int argc, char **argv)
{
    const char *in = rs_opt(argc, argv, "--cphd");
    const char *out = rs_opt(argc, argv, "--out");
    if (!in || !out) {
        printf("usage: resonarsat focus --cphd FILE --out FILE.png\n"
               "                       [--size N] [--cell M] [--offset X,Y | --at LAT,LON]\n"
               "                       [--pulse-start N --pulse-count N] [--raw FILE]\n"
               "                       [--rbins N] [--max-pulses N] [--pulse-stride N]\n"
               "                       [--range-taps N] [--no-optimize]\n"
               "\n"
               "--rbins reads a window of range bins and --max-pulses caps how many\n"
               "pulses are READ, not merely used. A large collect will not fit in\n"
               "memory otherwise: --pulse-count limits focusing but the whole\n"
               "collect is still loaded first.\n"
               "\n"
               "--pulse-stride keeps every nth pulse. It does NOT coarsen azimuth\n"
               "resolution: the aperture still spans the same dwell, so resolution\n"
               "is unchanged and only the sampling within it thins, folding azimuth\n"
               "ambiguities into the image. It lowers the effective PRF and so the\n"
               "observable vibration band, which is why no measurement run may use\n"
               "it.\n"
               "\n"
               "To coarsen resolution instead -- the honest way to match a large\n"
               "grid cell, since focusing at full resolution onto cells far larger\n"
               "than it aliases into ghost peaks that resemble targets -- shorten\n"
               "the aperture with --max-pulses. The warning this command prints\n"
               "names the cell the collect actually supports.\n"
               "\n"
               "The grid is centred on the scene reference point unless --offset\n"
               "moves it, in metres on the product's own image plane. A grid\n"
               "smaller than the collect otherwise silently covers only the\n"
               "middle of it. Output format follows the extension: .png or .pgm.\n"
               "\n"
               "--at LAT,LON places the grid on a known place instead, and is the\n"
               "safer of the two. The axes --offset speaks are the FILE's (CPHD's\n"
               "uIAX/uIAY), which point where the collector chose -- not north and\n"
               "east, and not azimuth and ground range. Converting coordinates into\n"
               "them by hand yields a number that is wrong in a way the output\n"
               "cannot reveal, because a misplaced grid still produces a complete,\n"
               "well-focused image OF THE WRONG GROUND. --at makes the product do\n"
               "the conversion, and prints the offset it resolved to.\n");
        return 1;
    }

    rs_cphd_t c;
    const size_t rbin_window = (size_t)rs_opt_double(argc, argv, "--rbins", 0.0);
    resonarsat_status_t st = rs_load_cphd(&c, in, rbin_window,
                       (size_t)rs_opt_double(argc, argv, "--pulse-stride", 0.0),
                       (size_t)rs_opt_double(argc, argv, "--max-pulses", 0.0));
    if (st != RS_OK) { rs_report_error("focus", st); return 1; }

    const size_t size = (size_t)rs_opt_double(argc, argv, "--size", 256);
    const double cell = rs_opt_double(argc, argv, "--cell", 1.0);

    rs_grid_t grid = { .n_x = size, .n_y = size,
                       .dx = cell, .dy = cell, .height = 0.0 };
    rs_parse_offset(rs_opt(argc, argv, "--offset"), grid.origin);
    if (rs_resolve_at(rs_opt(argc, argv, "--at"), &c.plane, grid.origin)) {
        rs_cphd_free(&c); return 1;
    }

    const size_t p_start = (size_t)rs_opt_double(argc, argv, "--pulse-start", 0);
    const size_t p_count = (size_t)rs_opt_double(argc, argv, "--pulse-count",
                                                 (double)c.n_pulse);

    rs_slc_t img;
    if ((st = rs_slc_alloc(&img, grid.n_x, grid.n_y)) != RS_OK) {
        rs_report_error("focus", st); rs_cphd_free(&c); return 1;
    }

    rs_warn_sampling(&c, p_count, cell);
    printf("backprojecting %zu pulses onto %zux%zu cells of %.2f m ...\n",
           p_count, grid.n_y, grid.n_x, cell);

    const rs_focus_opts_t fopts = {
        .single_thread = rs_opt_no_optimize(argc, argv),
        .range_taps = (int)rs_opt_double(argc, argv, "--range-taps", 0.0)
    };
    if (fopts.range_taps >= 4) {
        printf("range interpolation: %d-tap windowed sinc "
               "(default is 2-tap linear)\n", fopts.range_taps);
    }
    st = rs_focus_backproject_opts(&c, &grid, p_start, p_count, &fopts, &img);
    if (st != RS_OK) {
        rs_report_error("focus", st);
        rs_slc_free(&img); rs_cphd_free(&c);
        return 1;
    }

    /* The quicklook is a 40 dB log stretch clipped at the 99th percentile,
     * which is right for looking at a scene and wrong for measuring one: bright
     * scatterers saturate, so a point-spread function measured off it comes out
     * broader than it is, by an unknown factor. --raw writes the amplitudes
     * themselves so resolution can be measured rather than estimated. */
    const double dyn = rs_opt_double(argc, argv, "--dyn-range", 40.0);
    st = rs_raster_write_quicklook(&img, out, dyn);

    const char *raw = rs_opt(argc, argv, "--raw");
    if (raw) {
        double *amp = malloc(img.n_az * img.n_rg * sizeof *amp);
        if (amp) {
            for (size_t i = 0; i < img.n_az * img.n_rg; i++) {
                amp[i] = (double)cabsf(img.data[i]);
            }
            /* The tag rides in the axis description because that is the only
             * free-text field the .hdr sidecar has. A raw cube is the input to
             * every measurement made off this image, so which arithmetic produced
             * it has to travel with the pixels -- a bare .f32 carries nothing. */
            if (rs_raster_write_cube(amp, img.n_az, img.n_rg, 1, raw,
                                     fopts.single_thread
                                       ? "azimuth, range, amplitude [UNOPTIMIZED]"
                                       : "azimuth, range, amplitude") == RS_OK) {
                printf("wrote %s (%zu x %zu float32 amplitudes) and %s.hdr\n",
                       raw, img.n_az, img.n_rg, raw);
            }
            free(amp);
        }
    }
    if (st != RS_OK) rs_report_error("focus", st);
    else printf("wrote %s\n", out);

    rs_slc_free(&img);
    rs_cphd_free(&c);
    return st == RS_OK ? 0 : 1;
}

/* Warn that a shuffled null floor does not bound a phase measurement.
 *
 * "Only the ordering is destroyed" holds for an observable each look carries on
 * its own. Phase is unwrapped ACROSS looks, so a permutation sets
 * non-consecutive looks side by side -- exactly where the series steps furthest
 * -- and inflates the per-step noise the test exists to hold constant. Measured
 * on the Giza collect at 128 looks and 0.99 overlap: median largest step 0.052
 * rad in order against 1.878 rad shuffled, a factor of 36. A drifting series
 * then beats its own shuffles by construction, and one did, at p = 0.03, while
 * eight motionless simulations reproduced its frequency at 99 percent of its
 * prominence.
 *
 * Printed rather than refused, because the floor is still worth seeing and a
 * caller may want it for comparison. What must not happen is a reader taking it
 * for a bound it is not, which is why this sits next to the number rather than
 * only in the header. See rs_microm_estimator_t and
 * runs/giza/2026-07-30-uniform-phase-khufu/. */
static void rs_warn_shuffle_null_on_phase(rs_microm_estimator_t est)
{
    if (est != RS_MICROM_EST_PHASE) return;
    printf("  WARNING: this is a PHASE measurement, and a shuffle does not bound\n"
           "           one. Unwrapping runs ACROSS looks, so reordering inflates\n"
           "           the per-step noise the test is meant to hold fixed -- 36x\n"
           "           on the Giza collect -- and a drifting series beats its own\n"
           "           shuffles whatever it contains. This floor is not a bound.\n"
           "           Use --null-static, which a motionless scene carries the\n"
           "           same overlap, unwrap and detrend through and so cannot\n"
           "           walk over.\n");
}

/* Permute the time order of a sub-aperture stack in place, for a null test.
 *
 * On synthetic data a null test is easy: generate the scene again with nothing
 * moving and see what the pipeline claims anyway. Real data offers no such
 * control -- the motion, whatever it is, cannot be switched off. Shuffling the
 * looks supplies the missing control from the other direction. Scene content,
 * brightness, coherence, speckle, geometry and the number of looks are all
 * exactly preserved; the only thing destroyed is the order in which the looks
 * were taken, and that ordering is the entire basis of a vibration measurement.
 *
 * So a peak that survives the shuffle was never temporal. It is the spectrum of
 * a fixed spatial pattern being read out in some order, and the prominence it
 * scores is the floor that a genuine detection has to clear. This is the same
 * question bug 15 was about, asked where no second scene can be generated.
 *
 * 'seed' selects the permutation so a run can be repeated. The centre times are
 * deliberately left alone: the sampling grid must stay uniform, or the spectrum
 * would change for a second reason and the test would no longer be clean. */
static void rs_shuffle_looks(rs_subap_stack_t *stack, unsigned seed)
{
    if (!stack || stack->n_looks < 2) return;
    unsigned st = seed ? seed : 1u;
    for (size_t i = stack->n_looks - 1; i > 0; i--) {
        /* xorshift keeps the permutation reproducible across platforms, which
         * rand() would not. */
        st ^= st << 13; st ^= st >> 17; st ^= st << 5;
        const size_t j = (size_t)(st % (unsigned)(i + 1));
        const rs_slc_t tmp = stack->look[i];
        stack->look[i] = stack->look[j];
        stack->look[j] = tmp;
    }
}

/* Measure the null floor by repeated shuffling, and report the margin.
 *
 * NOT VALID FOR RS_MICROM_EST_PHASE. This holds everything but the time order
 * constant only for an observable that reads each look independently. A phase
 * series is unwrapped ACROSS looks, so reordering -- which sets non-consecutive
 * looks side by side -- inflates the very per-step noise the test is supposed to
 * preserve: 0.052 rad in order against 1.878 rad shuffled on the Giza collect at
 * 128 looks and 0.99 overlap, a factor of 36. A drifting phase series then beats
 * its own shuffles by construction, and one did, at p = 0.03, while a motionless
 * simulation reproduced its frequency at 99 percent of its prominence. Use
 * rs_null_static() for phase. See rs_microm_estimator_t and
 * runs/giza/2026-07-30-uniform-phase-khufu/.
 *
 * Loads nothing and focuses nothing: the caller's stack is reused, so a trial
 * costs a track and a spectrum rather than another read and range compression
 * of a multi-gigabyte collect. That matters because a floor is only as good as
 * the number of samples behind it, and a floor built from a handful of trials
 * misleads. Thirteen trials of one Melbourne configuration put every null below
 * the real measurement; forty-nine put one above it, and the margin over the
 * worst null went from 1.08 to 0.96. The tail is the whole question, so trials
 * have to be cheap enough to run in bulk.
 *
 * Each trial shuffles again from wherever the previous one left the stack. Any
 * permutation is as good as any other for this purpose, so there is no need to
 * restore the original order between trials.
 *
 * Returns RS_OK and writes the mean, standard deviation and maximum prominence
 * over 'trials' shuffles, plus how many reached 'real'. */
static resonarsat_status_t rs_null_floor(rs_subap_stack_t *stack,
                                         const rs_microm_params_t *mp,
                                         size_t trials, double real, double f_min,
                                         double *mean, double *sd, double *max_out,
                                         size_t *n_ge)
{
    if (!stack || !mp || trials == 0) return RS_ERR_ARG;

    double sum = 0.0, sum2 = 0.0, hi = 0.0;
    size_t n = 0, ge = 0;

    for (size_t i = 0; i < trials; i++) {
        rs_shuffle_looks(stack, (unsigned)(i * 2654435761u + 1u));

        rs_microm_t m;
        if (rs_microm_track(stack, mp, &m) != RS_OK) continue;

        rs_spectrum_t sp;
        if (rs_spectrum_compute_band(&m, RS_SPEC_VELOCITY, f_min, &sp) == RS_OK) {
            size_t best = 0;
            double prom = 0.0;
            const resonarsat_status_t bst = rs_spectrum_best_window(&sp, &best, &prom, NULL);
            /* A trial in which nothing cleared the quantisation floor is a
             * trial that produced NO detection, and that is a sample of the
             * null distribution rather than a missing one. Recording it as zero
             * keeps it in the count; dropping it would retain only the trials
             * that happened to produce something and inflate the null. */
            if (bst == RS_OK || bst == RS_ERR_RANGE) {
                if (bst != RS_OK) prom = 0.0;
                sum += prom;
                sum2 += prom * prom;
                if (prom > hi) hi = prom;
                if (prom >= real) ge++;
                n++;
            }
            rs_spectrum_free(&sp);
        }
        rs_microm_free(&m);
    }

    if (n == 0) return RS_ERR_SINGULAR;
    *mean = sum / (double)n;
    const double var = sum2 / (double)n - (*mean) * (*mean);
    *sd = var > 0.0 ? sqrt(var) : 0.0;
    *max_out = hi;
    *n_ge = ge;
    return RS_OK;
}

/* Measure the null floor from a STATIC SIMULATED SCENE, not from a shuffle.
 *
 * The shuffled floor above answers "is there temporal structure". That is the
 * wrong question for an overlapping decomposition, which manufactures temporal
 * structure regardless of the ground: adjacent sub-looks sharing most of their
 * bandwidth share most of their speckle, so their shifts are correlated before
 * anything moves, the series is smooth by construction, and shuffling destroys
 * exactly that smoothness. The unshuffled series then wins for a reason that is
 * not motion.
 *
 * This floor answers "is there more than a motionless world would give through
 * this same processing". Each trial synthesises phase history for static
 * scatterers over the reference collect's own geometry and runs the identical
 * chain -- same sub-aperture route, same tracker, same estimator -- so overlap,
 * tracker bias and estimator behaviour are all inherited. Only the motion is
 * missing.
 *
 * It costs far more than a shuffle: every trial refocuses. That is the price of
 * a floor that a heavily overlapped decomposition cannot walk over.
 *
 * Returns RS_OK and writes the mean, standard deviation and maximum prominence
 * over 'trials' realisations, plus how many reached 'real'. */
static resonarsat_status_t rs_null_static(const rs_cphd_t *ref, const rs_grid_t *grid,
                                          rs_subap_params_t *sp, const char *route,
                                          const rs_microm_params_t *mp,
                                          size_t trials, double real, double f_min,
                                          size_t sim_rbin, double extent_m,
                                          double *mean, double *sd, double *max_out,
                                          size_t *n_ge)
{
    if (!ref || !grid || !sp || !mp || trials == 0) return RS_ERR_ARG;

    double sum = 0.0, sum2 = 0.0, hi = 0.0;
    size_t n = 0, ge = 0;

    for (size_t i = 0; i < trials; i++) {
        rs_cphd_t sim;
        const double centre[2] = { grid->origin[0], grid->origin[1] };
        if (rs_simulate_static_like(ref, (unsigned)(i + 1), 0, centre, extent_m,
                                    sim_rbin, &sim) != RS_OK) {
            continue;
        }

        rs_subap_stack_t st;
        if (rs_build_subaps(&sim, grid, sp, route, &st) == RS_OK) {
            rs_microm_t m;
            if (rs_microm_track(&st, mp, &m) == RS_OK) {
                rs_spectrum_t spx;
                if (rs_spectrum_compute_band(&m, RS_SPEC_VELOCITY, f_min, &spx) == RS_OK) {
                    size_t best = 0;
                    double prom = 0.0;
                    const resonarsat_status_t bst =
                        rs_spectrum_best_window(&spx, &best, &prom, NULL);
                    /* As above: nothing resolved is a null sample of zero, not
                     * a trial to discard. */
                    if (bst == RS_OK || bst == RS_ERR_RANGE) {
                        if (bst != RS_OK) prom = 0.0;
                        sum += prom; sum2 += prom * prom;
                        if (prom > hi) hi = prom;
                        if (prom >= real) ge++;
                        n++;
                        if (bst == RS_OK) {
                            fprintf(stderr, "  static trial %zu/%zu: prominence %.1f "
                                            "at %.3f Hz\n",
                                    i + 1, trials, prom, spx.dominant_freq[best]);
                        } else {
                            fprintf(stderr, "  static trial %zu/%zu: nothing resolved "
                                            "above the quantisation floor\n",
                                    i + 1, trials);
                        }
                    }
                    rs_spectrum_free(&spx);
                }
                rs_microm_free(&m);
            }
            rs_subap_stack_free(&st);
        }
        rs_cphd_free(&sim);
    }

    if (n == 0) return RS_ERR_SINGULAR;
    *mean = sum / (double)n;
    const double var = sum2 / (double)n - (*mean) * (*mean);
    *sd = var > 0.0 ? sqrt(var) : 0.0;
    *max_out = hi;
    *n_ge = ge;
    return RS_OK;
}

/* Run the full micro-motion chain on a phase-history file.
 *
 * Focus, decompose into sub-looks, track, and estimate spectra. Prints the
 * observable band and the resolution it cost on the same line, since reporting
 * either alone is how the trade gets misread. */
static int rs_cmd_mmotion(int argc, char **argv)
{
    const char *in = rs_opt(argc, argv, "--cphd");
    if (!in) {
        printf("usage: resonarsat mmotion --cphd FILE [--n N] [--overlap F]\n"
               "                          [--offset X,Y | --at LAT,LON]\n"
               "                          [--subap pulse|uniform|paper]\n"
               "                          [--no-detrend] [--null-static N]\n"
               "                          [--estimator correlation|phase|splitband]\n"
               "                          [--shuffle-looks SEED] [--null-trials N]\n"
               "                          [--fmin HZ]\n"
               "                          [--reference first|adjacent|pair]\n"
               "                          [--b-shift HZ] [--shifts FILE.csv]\n"
               "                          [--upsample N]\n"
               "                          [--size N] [--cell M] [--win N]\n"
               "                          [--coherence F] [--out PREFIX]\n"
               "                          [--ccd-out PREFIX] [--ccd-win N]\n"
               "                          [--ccd-loading F]\n"
               "                          [--no-optimize]\n"
               "\n"
               "--ccd-out runs the scale-invariant change-detection LOCATOR over\n"
               "the same sub-aperture stack and writes a map. It answers 'where in\n"
               "this scene is something moving', not 'at what frequency' -- a\n"
               "different question from everything else this command does, and one\n"
               "that does not need the tracker to succeed. The statistic's\n"
               "no-change value is 1.0; a bright but STATIONARY target scores near\n"
               "1.0 too, which is the point of it. It has no detection threshold:\n"
               "compare a map against one from a motionless scene before reading\n"
               "structure into it. --ccd-win sets the sliding window (default 5)\n"
               "and --ccd-loading the noise floor as a fraction of mean scene\n"
               "power (default 1e-3; zero shows the unregularised behaviour).\n"
               "\n"
               "--no-optimize is an audit baseline, not a better measurement. It\n"
               "searches the WHOLE upsampled correlation surface for the peak instead\n"
               "of the neighbourhood of the integer peak, and runs serially. Only the\n"
               "first of those can change a number; backprojection is bitwise identical\n"
               "either way (see rs_focus_opts_t). Measured cost is 2.5x per correlator\n"
               "call at these defaults plus about 4x for losing the threads -- single\n"
               "digits, so it is affordable on a full-scale scene.\n"
               "\n"
               "--coherence masks windows whose sub-looks do not correlate (default\n"
               "0.4). Isolated point targets on an empty scene score below that even\n"
               "when tracking perfectly; pass 0 to inspect an unmasked result.\n"
               "\n"
               "--reference selects which images each correlation is taken between.\n"
               "'first' (default) compares every look to look 0. 'pair' compares each\n"
               "look's slave to its own master, the two held --b-shift apart and swept\n"
               "together, which is what WO2024008365A1 and the Giza paper describe.\n"
               "'adjacent' accumulates consecutive differences.\n"
               "\n"
               "NEITHER 'pair' NOR 'adjacent' RECOVERS A FREQUENCY on the synthetic\n"
               "single-target fixture: both return the lowest spectral bin whatever is\n"
               "injected. 'pair' is exposed because it is what the sources describe and\n"
               "it should be testable, not because it works. Do not read a measurement\n"
               "out of either. See rs_microm_ref_t.\n"
               "\n"
               "--b-shift sets the master-slave separation in Hz and needs --subap\n"
               "paper and --reference pair. Zero derives it from the sweep step. The\n"
               "gap is a time lag, so a larger --b-shift observes a LOWER mechanical\n"
               "frequency. 'pair' measures a difference across that lag, which\n"
               "high-passes the series: amplitudes are attenuated by |2 sin(pi f dt)|\n"
               "and are not comparable to 'first'.\n"
               "\n"
               "AT THE DEFAULT --b-shift THE PAIR IS DEGENERATE: the separation equals\n"
               "the sweep step, so each slave IS its neighbour's master, sample for\n"
               "sample, and 'pair' reduces to differencing consecutive looks. The band\n"
               "layout holds one step of headroom, so any --b-shift that makes the\n"
               "slave a distinct band must be SMALLER than the step -- a shorter lag,\n"
               "and a weaker difference. Larger values are refused.\n"
               "\n"
               "--shifts writes the raw per-look shift series to CSV before detrending\n"
               "or any spectral estimation. Use it with sim_cphd --clutter to tell a\n"
               "real low-frequency motion from correlator bias: the spectrum cannot\n"
               "separate them, the series can.\n"
               "\n"
               "--upsample N locates the correlation peak to 1/N of a pixel, which is\n"
               "also the QUANTISATION of the reported series. A motion whose excursion\n"
               "is under one step comes back as a two-level series whose energy sits at\n"
               "low frequency regardless of what drove it. Check the peak-to-peak in\n"
               "--shifts against 1/N before believing any lowest-bin result.\n");
        return 1;
    }

    const char *prefix = rs_opt(argc, argv, "--out");

    rs_cphd_t c;
    const size_t rbin_window = (size_t)rs_opt_double(argc, argv, "--rbins", 0.0);
    resonarsat_status_t st = rs_load_cphd(&c, in, rbin_window,
                       (size_t)rs_opt_double(argc, argv, "--pulse-stride", 0.0),
                       (size_t)rs_opt_double(argc, argv, "--max-pulses", 0.0));
    if (st != RS_OK) { rs_report_error("mmotion", st); return 1; }

    const size_t size = (size_t)rs_opt_double(argc, argv, "--size", 128);
    const double cell = rs_opt_double(argc, argv, "--cell", 1.0);
    rs_grid_t grid = { .n_x = size, .n_y = size,
                       .dx = cell, .dy = cell, .height = 0.0 };
    rs_parse_offset(rs_opt(argc, argv, "--offset"), grid.origin);
    if (rs_resolve_at(rs_opt(argc, argv, "--at"), &c.plane, grid.origin)) {
        rs_cphd_free(&c); return 1;
    }

    rs_warn_sampling(&c, c.n_pulse, cell);

    /* Resolved before the sub-apertures are built, because that is the stage it
     * first has to reach. */
    const int no_optimize = rs_opt_no_optimize(argc, argv);

    rs_subap_params_t sp;
    rs_subap_params_default(&sp);
    sp.n_looks = (size_t)rs_opt_double(argc, argv, "--n", 16);
    sp.overlap = rs_opt_double(argc, argv, "--overlap", 0.40);
    sp.single_thread = no_optimize;
    sp.range_taps = (int)rs_opt_double(argc, argv, "--range-taps", 0.0);
    if (sp.range_taps >= 4) {
        printf("range interpolation: %d-tap windowed sinc\n", sp.range_taps);
    }
    const int ref_mode = rs_parse_reference(argc, argv, &sp);
    if (ref_mode < 0) { rs_cphd_free(&c); return 1; }
    const char *subap_route = rs_opt(argc, argv, "--subap");
    if (!subap_route && ref_mode == RS_MICROM_REF_PAIR) subap_route = "paper";

    /* Pair mode is defined by two Doppler filters and therefore selects the
     * paper route when the caller did not choose a route explicitly. */
    rs_subap_stack_t stack;
    if ((st = rs_build_subaps(&c, &grid, &sp, subap_route, &stack)) != RS_OK) {
        rs_report_error("mmotion", st); rs_cphd_free(&c); return 1;
    }

    rs_microm_params_t mp;
    rs_microm_params_default(&mp);
    mp.reference = (rs_microm_ref_t)ref_mode;
    mp.no_optimize = no_optimize;
    {
        /* --estimator selects WHAT is measured, not merely how well. Phase and
         * correlation live in different regimes; see rs_microm_estimator_t.
         *
         * Read before the shuffle below rather than after it, because whether a
         * shuffled floor bounds anything depends on the answer, and the caveat
         * belongs beside the notice rather than several screens later. */
        const char *est = rs_opt(argc, argv, "--estimator");
        if (est) {
            if (strcmp(est, "phase") == 0)            mp.estimator = RS_MICROM_EST_PHASE;
            else if (strcmp(est, "correlation") == 0) mp.estimator = RS_MICROM_EST_CORRELATION;
            else if (strcmp(est, "splitband") == 0)   mp.estimator = RS_MICROM_EST_SPLITBAND;
            else fprintf(stderr, "warning: unknown --estimator '%s'; "
                                 "using correlation\n", est);
        }
    }

    const unsigned shuffle = (unsigned)rs_opt_double(argc, argv, "--shuffle-looks", 0.0);
    if (shuffle) {
        rs_shuffle_looks(&stack, shuffle);
        printf("NULL TEST: sub-look time order shuffled with seed %u. Scene, coherence\n"
               "  and geometry are unchanged; only the ordering is destroyed. Whatever\n"
               "  prominence appears below is the floor, not a detection.\n", shuffle);
        rs_warn_shuffle_null_on_phase(mp.estimator);
    }

    printf("sub-apertures: %zu looks, dt %.4f s\n", stack.n_looks, stack.dt);
    printf("  observable band  f_max %.2f Hz   AT sub-look resolution %.2f m\n",
           stack.f_max, stack.az_resolution);

    /* The change-detection locator, run here rather than beside the other
     * outputs at the end of this function.
     *
     * Placement is deliberate. It needs only the sub-aperture stack, and every
     * Giza run so far has exited at the spectrum stage with RS_ERR_RANGE --
     * honestly, because no window resolved motion -- which is upstream of the
     * --out block. A locator written down there would never run on precisely
     * the scenes it exists to say something about. See rs_ccd_t. */
    {
        const char *ccd_out = rs_opt(argc, argv, "--ccd-out");
        if (ccd_out) {
            rs_ccd_params_t cp;
            rs_ccd_params_default(&cp);
            const double cw = rs_opt_double(argc, argv, "--ccd-win", 0.0);
            if (cw > 0.0) cp.win = (size_t)cw;
            const double cl = rs_opt_double(argc, argv, "--ccd-loading", -1.0);
            if (cl >= 0.0) cp.loading = cl;

            rs_ccd_t ccd;
            if (rs_ccd_locate(&stack, &cp, &ccd) != RS_OK) {
                rs_report_error("mmotion", RS_ERR_ARG);
            } else {
                /* Reported rather than only written, because the map's absolute
                 * level is the whole of its meaning: 1.0 is the no-change value,
                 * and a map whose median sits there has found nothing however
                 * structured its picture looks. */
                double lo = 0.0, hi = 0.0, sum = 0.0;
                size_t n = 0;
                for (size_t p = 0; p < ccd.n_row * ccd.n_col; p++) {
                    const double v = ccd.map[p];
                    /* The border is not computed and holds zero; seeding the
                     * extremes from it reported a minimum of 0.000 for every
                     * scene, which reads as a pixel where the statistic
                     * collapsed rather than one that was never evaluated. */
                    if (v <= 0.0) continue;
                    if (n == 0 || v < lo) lo = v;
                    if (n == 0 || v > hi) hi = v;
                    sum += v; n++;
                }
                printf("CCD locator: %zux%zu window over %zu sub-aperture triples\n",
                       cp.win, cp.win, ccd.n_triples);
                printf("  statistic  min %.3f  mean %.3f  max %.3f   "
                       "(1.000 is the no-change value)\n",
                       lo, n ? sum / (double)n : 0.0, hi);
                printf("  A MAP IS NOT EVIDENCE WITHOUT A FLOOR. The source method\n"
                       "  implements no detection threshold; compare against a\n"
                       "  motionless scene through this same chain before reading\n"
                       "  structure into it.\n");

                char path[512];
                snprintf(path, sizeof path, "%s_ccd.png", ccd_out);
                rs_raster_write_map(ccd.map, ccd.n_row, ccd.n_col,
                                    path, 0.0, 0.0, RS_PALETTE_VIRIDIS);
                snprintf(path, sizeof path, "%s_ccd.f32", ccd_out);
                rs_raster_write_cube(ccd.map, 1, ccd.n_row, ccd.n_col, path,
                                     "axes row (grid x), col (grid y); "
                                     "scale-invariant CCD statistic, 1 = no change");
                printf("wrote %s_ccd.png and %s_ccd.f32\n", ccd_out, ccd_out);
                rs_ccd_free(&ccd);
            }
        }
    }

    mp.win_az = mp.win_rg = (size_t)rs_opt_double(argc, argv, "--win", 32);
    mp.stride_az = mp.stride_rg = mp.win_az / 2;
    mp.coherence_min = rs_opt_double(argc, argv, "--coherence", 0.4);

    /* Sub-pixel refinement factor, exposed because it is a FLOOR on what can be
     * measured and not merely a precision setting.
     *
     * The correlation peak is located to 1/upsample of a pixel, so the shift
     * series is quantised at that step. Any motion whose apparent excursion is
     * smaller than one step comes back as a series that is constant, or that
     * flips between two adjacent levels -- and a two-level series has its
     * energy at low frequency whatever drove the flips. That is
     * indistinguishable, in the spectrum, from a drift.
     *
     * It binds hardest on --reference pair. That observable is a difference
     * across one fixed lag, so its excursion is smaller than the full
     * displacement by roughly 2*sin(pi*f*dt), and at the default upsampling a
     * pair series can sit entirely inside one quantisation step while the same
     * scene measured against a fixed reference does not. Raising this is the
     * first thing to try before concluding the pair carries no signal.
     *
     * The cost is quadratic in the refinement window, so it is left at the
     * published default rather than raised globally. */
    {
        const double up = rs_opt_double(argc, argv, "--upsample", 0.0);
        if (up > 0.0) {
            mp.upsample_az = (size_t)up;
            mp.upsample_rg = (size_t)up;
            printf("sub-pixel refinement: 1/%zu px "
                   "(default 1/%d azimuth, 1/%d range)\n",
                   mp.upsample_az, 10, 20);
        }
    }

    rs_microm_t m;
    if ((st = rs_microm_track(&stack, &mp, &m)) != RS_OK) {
        rs_report_error("mmotion", st);
        rs_subap_stack_free(&stack); rs_cphd_free(&c);
        return 1;
    }
    size_t n_pass = 0;
    for (size_t w = 0; w < m.n_win; w++) if (m.quality[w] >= mp.coherence_min) n_pass++;
    printf("tracked %zu windows (%zu x %zu); %zu pass the %.2f coherence mask\n",
           m.n_win, m.n_win_az, m.n_win_rg, n_pass, mp.coherence_min);
    if (n_pass == 0) {
        printf("  no window is coherent enough to carry a measurement. Any\n"
               "  frequency reported below is tracking noise, not a mode.\n");
    }

    /* Raw per-look shifts, dumped BEFORE any spectral estimation.
     *
     * WHY THIS IS SEPARATE FROM EVERY OTHER OUTPUT. Everything else this command
     * prints has been through rs_spectrum_compute(): detrended, windowed,
     * periodogrammed, and reduced to one dominant frequency per window. That
     * chain is where a weak observable becomes indistinguishable from a biased
     * one -- a slowly varying correlator bias and a genuine low-frequency motion
     * both end up as energy in the lowest bins, and no amount of staring at the
     * spectrum separates them.
     *
     * The shift series itself does separate them. A real sinusoid at f is
     * visible in the raw series as a sinusoid; correlator bias is monotone or
     * smoothly curved across the sweep, because it follows the sub-look's own
     * point response rather than the scene. Writing the series out is what makes
     * that a matter of looking rather than of argument.
     *
     * This is the measurement RS_MICROM_REF_PAIR needs. Its samples are first
     * differences across a fixed lag, so they are small by construction, and
     * whether they are dominated by bias is precisely the open question in
     * rs_microm_ref_t. Pair it with sim_cphd --clutter, which removes the
     * empty-window confound that makes the bias largest. */
    {
        const char *shift_out = rs_opt(argc, argv, "--shifts");
        if (shift_out) {
            FILE *sf = fopen(shift_out, "w");
            if (!sf) {
                rs_report_error("mmotion", RS_ERR_IO);
            } else {
                fprintf(sf, "# raw tracked shifts, before detrend or spectrum\n");
                fprintf(sf, "# reference=%s b_shift_hz=%.12g pair_lag_s=%.12g "
                            "dt_s=%.12g looks=%zu\n",
                        (mp.reference == RS_MICROM_REF_PAIR) ? "pair" :
                        (mp.reference == RS_MICROM_REF_ADJACENT) ? "adjacent" : "first",
                        stack.b_shift_hz, stack.pair_lag_s, stack.dt, stack.n_looks);
                /* disp_los and phase are written too, because without them this
                 * file cannot describe a phase run at all: that estimator
                 * leaves disp_az and disp_rg identically zero, so the three
                 * original columns showed nothing but the differenced velocity
                 * and a diagnostic dump could not see the primary observable. */
                fprintf(sf, "window,look,centre_time_s,disp_az_px,disp_rg_px,"
                            "vel_los_ms,disp_los_m,phase_rad,quality\n");
                for (size_t w = 0; w < m.n_win; w++) {
                    for (size_t k = 0; k < m.n_looks; k++) {
                        const size_t idx = w * m.n_looks + k;
                        fprintf(sf, "%zu,%zu,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.6f\n",
                                w, k,
                                (k < stack.n_looks) ? stack.centre_time[k] : 0.0,
                                m.disp_az[idx], m.disp_rg[idx], m.vel_los[idx],
                                m.disp_los[idx], m.phase[idx],
                                m.quality[w]);
                    }
                }
                fclose(sf);
                printf("wrote %s (%zu windows x %zu looks, pre-spectrum)\n",
                       shift_out, m.n_win, m.n_looks);
            }
        }
    }

    const double f_min = rs_opt_double(argc, argv, "--fmin", 0.0);
    if (f_min > 0.0) {
        printf("band floor: ignoring bins below %.3f Hz when picking the peak\n", f_min);
    }

    /* WHICH OBSERVABLE THE SPECTRUM IS TAKEN OF DEPENDS ON THE ESTIMATOR, and
     * getting it wrong costs a whole spectrum's worth of meaning.
     *
     * Correlation measures WHERE a patch sits, and a radially moving target is
     * displaced in azimuth by dx = R*v/V -- so the tracked shift is
     * proportional to VELOCITY, and 'vel_los' is the natural observable.
     *
     * Phase measures displacement directly, d = -psi*lambda/(4*pi). Feeding the
     * spectrum 'vel_los' there differentiates it first, and a derivative
     * multiplies each Fourier component by its own frequency -- which turns
     * flat noise into blue noise and puts the peak at the top of the band for
     * no physical reason. Measured on the Giza control: the same series gives a
     * median dominant frequency of 20.7 Hz through velocity, with 47 of 49
     * windows above 12 Hz, and 0.53 Hz through displacement, with none above
     * 12 Hz. The 20.7 Hz was an artefact of differencing. */
    const rs_spectrum_source_t src = (mp.estimator == RS_MICROM_EST_PHASE)
                                   ? RS_SPEC_DISPLACEMENT : RS_SPEC_VELOCITY;
    if (mp.estimator == RS_MICROM_EST_PHASE) {
        printf("spectrum taken of line-of-sight DISPLACEMENT, which is what the "
               "phase estimator measures directly\n");
    }

    rs_spectrum_t spec;
    if ((st = rs_spectrum_compute_opts(&m, src, f_min,
                                 rs_opt_flag(argc, argv, "--no-detrend")
                                     ? RS_DETREND_NONE : RS_DETREND_LINEAR,
                                 &spec)) != RS_OK) {
        rs_report_error("mmotion", st);
        rs_microm_free(&m); rs_subap_stack_free(&stack);
        rs_cphd_free(&c);
        return 1;
    }

    printf("spectra: %zu bins, %.4f Hz resolution\n", spec.n_freq, spec.df);

    /* Report the window with the most prominent spectral peak.
     *
     * Not the largest displacement excursion: that selects the NOISIEST window,
     * because noise excursions exceed real ones, and it is why this command
     * previously returned the same answer whatever motion was injected. Not the
     * highest tracking coherence either: that selects static ground. Prominence
     * asks whether a window's motion concentrates at one frequency, which is the
     * question being posed. */
    size_t best = 0;
    double prom = 0.0;
    size_t n_cand = 0;
    {
        /* The real status, not a substituted one. This reported RS_ERR_ARG
         * whatever went wrong, which hid the message the failing call had
         * already written -- and the interesting case now is RS_ERR_RANGE,
         * meaning no window resolved motion above the tracker's own
         * resolution. That is a result about the scene, so it is worth saying
         * plainly rather than as a generic argument error. */
        const resonarsat_status_t bst =
            rs_spectrum_best_window(&spec, &best, &prom, &n_cand);
        if (bst != RS_OK) {
            rs_report_error("mmotion", bst);
            if (bst == RS_ERR_RANGE) {
                fprintf(stderr,
                    "  No frequency is reported because none is supported by the\n"
                    "  data, which is a different statement from finding none.\n"
                    "  Raise --upsample to resolve a smaller excursion, or check\n"
                    "  with --shifts whether the series moves at all.\n");
            }
            rs_spectrum_free(&spec); rs_microm_free(&m);
            rs_subap_stack_free(&stack); rs_cphd_free(&c);
            return 1;
        }
    }

    double pp = 0.0;
    {
        double lo = m.vel_los[best * m.n_looks], hi = lo;
        for (size_t k = 1; k < m.n_looks; k++) {
            const double v = m.vel_los[best * m.n_looks + k];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        pp = hi - lo;
    }

    printf("strongest peak in window %zu: %.3f Hz, prominence %.1f, "
           "quality %.3f, peak-to-peak velocity %.1f mm/s\n",
           best, spec.dominant_freq[best], prom, spec.quality[best], pp * 1e3);

    /* How many windows were eligible, printed beside the winner rather than
     * left implicit.
     *
     * The floor is three sigma for ONE window and is applied to every window
     * independently, so the chance that something crosses it on quantisation
     * noise alone rises with the grid. On the Giza collect the identical chain
     * refused honestly at 225 windows and reported "0.183 Hz, prominence 29.9"
     * at 961, off two crossings -- the observable was the same in both, only
     * the number of tries changed. A reader given the count sees that; a reader
     * given the prominence alone does not. See rs_spectrum_best_window(). */
    /* Name the gates that actually ran. The count is taken after BOTH the
     * relative coherence gate and the quantisation floor, and the floor is
     * skipped entirely when 'quant_px' is zero -- which is the phase
     * estimator's value, since it has no correlation surface. Reporting "126 of
     * 225 cleared the quantisation floor" on a phase run would name the one
     * gate that did not run: those 99 exclusions were the coherence gate. */
    printf("  %zu of %zu windows were eligible for selection (coherence gate%s)\n",
           n_cand, spec.n_win,
           (spec.quant_px > 0.0) ? " and quantisation floor"
                                 : "; the floor does not apply to this estimator");
    if (spec.quant_px > 0.0 && n_cand > 0 && n_cand < 4) {
        printf("  WARNING: windows overlap at the tracking stride, so a target\n"
               "           large enough to resolve falls in a 2x2 block of them\n"
               "           at minimum. Fewer than four qualifying windows cannot\n"
               "           describe a spatially resolved mode, and is the\n"
               "           signature of a threshold crossed by chance across\n"
               "           many tries rather than of motion.\n");
    }

    if (prom < 3.0) {
        printf("  WARNING: prominence %.1f is close to the flat-spectrum value of\n"
               "           1.0. No window holds motion concentrated at a single\n"
               "           frequency; the figure above is not a detection.\n", prom);
    }

    /* The observation ratio can only be checked once a frequency exists, since
     * it is the measurement that supplies the period. */
    const double eta = rs_spectrum_observation_ratio(stack.t_sap,
                                                     spec.dominant_freq[best]);
    const double resp = rs_spectrum_subaperture_response(stack.t_sap,
                                                         spec.dominant_freq[best]);
    printf("  sub-aperture response %.4f (%.1f dB) at an observation ratio of %.2f\n",
           resp, 20.0 * log10(resp > 1e-12 ? resp : 1e-12), eta);

    /* Amplitude attenuation, not invalidity. An earlier version called any
     * ratio above 0.5 unreliable, which would reject two accelerometer-
     * validated measurements in the literature at ratios of 18 and 36. */
    if (resp < 0.1) {
        printf("     Amplitudes here are attenuated by the sub-aperture window\n"
               "     and are underestimates; frequencies are not affected.\n"
               "     Note the response above assumes displacement averaging,\n"
               "     which is what this tracker measures. Phase-based micro-\n"
               "     Doppler is not bound by it: the literature recovers 36 Hz\n"
               "     at an observation ratio of exactly 18, where this model\n"
               "     predicts zero.\n");
    }

    /* The comb of blind frequencies, which a bare ratio hides entirely. */
    {
        const double near = fabs(eta - floor(eta + 0.5));
        if (eta > 0.5 && near < 0.15) {
            printf("     WARNING: that frequency sits within %.2f of an integer\n"
                   "     observation ratio, where the sub-aperture averages a whole\n"
                   "     number of cycles and the response is exactly zero. The\n"
                   "     blind frequencies here are multiples of %.3f Hz. A peak\n"
                   "     beside a null is more likely to be an artefact of the\n"
                   "     window than a mode of the scene -- for THIS observable.\n",
                   near, 1.0 / stack.t_sap);
        }
    }

    if (spec.dominant_freq[best] > 0.9 * stack.f_max) {
        printf("  WARNING: that is close to the %.2f Hz Nyquist limit; the true\n"
               "           frequency may be higher and aliased. Raise --n.\n",
               stack.f_max);
    }
    /* The static-scene floor, which an overlapping decomposition cannot walk
     * over the way it can walk over a shuffle. Costs a refocus per trial. */
    const size_t s_trials = (size_t)rs_opt_double(argc, argv, "--null-static", 0.0);
    if (s_trials > 0) {
        double nm = 0.0, nsd = 0.0, nmax = 0.0;
        size_t nge = 0;
        const double extent = 0.5 * (double)grid.n_x * grid.dx;
        printf("\nSTATIC-SCENE NULL FLOOR from %zu simulated motionless collects\n"
               "with this collect's own geometry, through the identical chain:\n",
               s_trials);
        if (rs_null_static(&c, &grid, &sp, rs_opt(argc, argv, "--subap"), &mp,
                           s_trials, prom, f_min, 1024, extent,
                           &nm, &nsd, &nmax, &nge) == RS_OK) {
            printf("  mean %.1f, sd %.1f, worst %.1f\n", nm, nsd, nmax);
            printf("  detection %.1f is %.2fx the mean and %.2fx the worst\n",
                   prom, nm > 0.0 ? prom / nm : 0.0, nmax > 0.0 ? prom / nmax : 0.0);
            printf("  %zu of %zu reached it -- empirical p = %.4f\n",
                   nge, s_trials, (double)(nge + 1) / (double)(s_trials + 1));
            if (nge > 0) {
                printf("  A MOTIONLESS SCENE REACHED THIS MEASUREMENT through the\n"
                       "  same processing. Whatever the peak is, it is not\n"
                       "  evidence of motion.\n");
            } else {
                printf("  No motionless realisation reached it.\n");
            }
        } else {
            rs_report_error("mmotion", RS_ERR_SINGULAR);
        }
    }

    const size_t trials = (size_t)rs_opt_double(argc, argv, "--null-trials", 0.0);
    if (trials > 0) {
        double nm = 0.0, nsd = 0.0, nmax = 0.0;
        size_t nge = 0;
        if (rs_null_floor(&stack, &mp, trials, prom, f_min,
                          &nm, &nsd, &nmax, &nge) == RS_OK) {
            printf("\nNULL FLOOR from %zu shuffles of the sub-look time order:\n", trials);
            printf("  mean %.1f, sd %.1f, worst %.1f\n", nm, nsd, nmax);
            printf("  detection %.1f is %.2fx the mean and %.2fx the worst null\n",
                   prom, nm > 0.0 ? prom / nm : 0.0, nmax > 0.0 ? prom / nmax : 0.0);
            printf("  %zu of %zu nulls reached it -- empirical p = %.4f\n",
                   nge, trials, (double)(nge + 1) / (double)(trials + 1));
            if (nge > 0 || prom < nmax) {
                printf("  A null matched or beat the measurement. At this operating point\n"
                       "  the detection is not distinguishable from shuffled noise.\n");
            }
            /* The warning matters most in the case that looks like success: a
             * phase run clearing every shuffle prints nothing above, and that
             * silence is what a reader would otherwise take for a bound. */
            rs_warn_shuffle_null_on_phase(mp.estimator);
        }
    }

    printf("  (amplitude is QUALITATIVE -- relative amplitudes do not validate\n"
           "   against ground truth even where frequencies do)\n");

    if (prefix) {
        char path[512];
        snprintf(path, sizeof path, "%s_freq.png", prefix);
        rs_raster_write_map(spec.dominant_freq, spec.n_win_az, spec.n_win_rg,
                            path, 0.0, 0.0, RS_PALETTE_VIRIDIS);
        snprintf(path, sizeof path, "%s_quality.png", prefix);
        rs_raster_write_map(spec.quality, spec.n_win_az, spec.n_win_rg,
                            path, 0.0, 1.0, RS_PALETTE_VIRIDIS);
        printf("wrote %s_freq.png and %s_quality.png\n", prefix, prefix);
    }

    rs_spectrum_free(&spec);
    rs_microm_free(&m);
    rs_subap_stack_free(&stack);
    rs_cphd_free(&c);
    return 0;
}

/* Focus vibration observations into a depth profile.
 *
 * Requires --velocity and --frequency with no defaults, and prints the full
 * assumption set alongside the result. Both requirements are deliberate: they
 * are what stops a depth axis being fabricated from library defaults. */
static int rs_cmd_tomo(int argc, char **argv)
{
    const char *in = rs_opt(argc, argv, "--cphd");
    if (!in) {
        printf("usage: resonarsat tomo --cphd FILE --velocity MS [--frequency HZ]\n"
               "                       (--frequency is required for models A, B and C;\n"
               "                        model D measures it from the spectrum)\n"
               "                       [--model A|B|C|D] [--solver dft|lstsq]\n"
               "                       [--convention paper|patent]\n"
               "                       [--depth M] [--cell M] [--offset X,Y | --at LAT,LON]\n"
               "                       [--palette gray|viridis|energy|jet]\n"
               "                       [--subap pulse|uniform|paper] [--no-detrend]\n"
               "                       [--aperture M] [--y shifts|los]\n"
               "                       [--no-window] [--keep-mean] [--regularisation F]\n"
               "                       [--eq22-literal-t VALUE (experimental)]\n"
               "                       [--geocode FILE.csv] [--patent-exact]\n"
               "                       [--null-align N] [--null-align-seed S]\n"
               "                       [--no-optimize]\n"
               "\n"
               "--null-align is THE DEPTH STAGE'S OWN NULL TEST, and the one to\n"
               "run before believing a tomogram. It circularly shifts each\n"
               "window's depth profile by an independent random amount and\n"
               "restacks. Every profile is preserved exactly -- same peaks, same\n"
               "widths, same artefacts of steering and windowing -- so the only\n"
               "thing destroyed is the depth at which each window's profile sits,\n"
               "and therefore whether the windows AGREE. A contrast that survives\n"
               "means the windows independently picked the same depth. One that\n"
               "collapses means the stack was reporting the average shape of\n"
               "per-window artefacts, which every scene has.\n"
               "\n"
               "A shuffle of the sub-look order is the wrong instrument here: it\n"
               "changes each window's profile as well as their agreement, so it\n"
               "cannot say which half carried the result.\n"
               "                       [--out FILE] [--section N --section-out FILE.png]\n"
               "\n"
               "--velocity and --frequency are REQUIRED and have no defaults.\n"
               "They are ASSUMED constants that scale the depth axis; nothing in\n"
               "this pipeline measures either. --frequency is the paper's\n"
               "'investigation frequency', not a measured vibration frequency.\n");
        return 1;
    }

    rs_tomo_params_t tp;
    rs_tomo_params_default(&tp);
    tp.velocity  = rs_opt_double(argc, argv, "--velocity", 0.0);
    tp.frequency = rs_opt_double(argc, argv, "--frequency", 0.0);
    tp.depth_max = rs_opt_double(argc, argv, "--depth", 60.0);
    tp.depth_cell = rs_opt_double(argc, argv, "--cell", 1.0);
    tp.eq22_literal_t = rs_opt_double(argc, argv, "--eq22-literal-t", 0.0);

    const char *model = rs_opt(argc, argv, "--model");
    if (model) {
        if (*model == 'B' || *model == 'b') tp.model = RS_TOMO_MODEL_B;
        else if (*model == 'C' || *model == 'c') tp.model = RS_TOMO_MODEL_C;
        else if (*model == 'D' || *model == 'd') tp.model = RS_TOMO_MODEL_D;
        else tp.model = RS_TOMO_MODEL_A;
    }
    const char *solver = rs_opt(argc, argv, "--solver");
    if (solver && strcmp(solver, "lstsq") == 0) tp.solver = RS_TOMO_SOLVER_LSTSQ;

    const char *conv = rs_opt(argc, argv, "--convention");
    if (conv && strcmp(conv, "patent") == 0) tp.convention = RS_WAVELEN_PATENT;

    /* Reject a missing assumption before reading gigabytes: velocity and
     * frequency have no defaults and their absence is the common mistake.
     *
     * Model D is exempt from the frequency requirement because it measures the
     * frequency rather than assuming it. This check duplicated the one in
     * rs_tomo_params_check() and was missed when that one was corrected, so
     * Model D remained unrunnable while the usage text said otherwise -- which
     * is the hazard of validating the same condition in two places. */
    resonarsat_status_t st;
    if (tp.velocity <= 0.0 ||
        (tp.model != RS_TOMO_MODEL_D && tp.frequency <= 0.0)) {
        rs_report_error("tomo", RS_ERR_ARG);
        return 1;
    }

    rs_cphd_t c;
    const size_t rbin_window = (size_t)rs_opt_double(argc, argv, "--rbins", 0.0);
    if ((st = rs_load_cphd(&c, in, rbin_window,
                       (size_t)rs_opt_double(argc, argv, "--pulse-stride", 0.0),
                       (size_t)rs_opt_double(argc, argv, "--max-pulses", 0.0))) != RS_OK) { rs_report_error("tomo", st); return 1; }

    /* Geometry comes from the collect, not from the simulator's defaults. This
     * has to happen after the load and before the check, because the check
     * tests the depth grid against the resolution this geometry supports. */
    rs_tomo_geometry_from_cphd(&tp, &c);
    /* Kept before the overrides so the report below can name what the collect
     * actually implied. Reading it afterwards, as an earlier version did,
     * printed "was X from the collect" where X had already been replaced by the
     * command line -- an audit trail that agreed with itself no matter what. */
    const double aperture_from_collect = tp.aperture;
    rs_tomo_geometry_overrides(&tp, argc, argv);
    int patent_exact = 0;
    /* --patent-exact: select the unconditioned patent-chain interpretation.
     *
     * The historical flag name does not mean every rendered symbol is literal.
     * Patent Eq. 22 prints an undefined, dimensionally inconsistent 2*pi*t
     * factor and conflicts with Eq. 23 about matrix orientation. Model A uses
     * the conventional k-by-F exp(j*Kz*z) repair; metadata states this.
     *
     * The conditioning steps this project adds all default ON, because each
     * removes a measured artefact -- without mean removal, 89 percent of
     * windows peak at zero depth. But a caller asking for the published method
     * should not have to know four separate flags to get it, and must not get
     * something else by default without being told.
     *
     * THE WAVELENGTH CONVENTION IS PART OF THIS. The patent computes
     * lambda = v/f and the Giza paper computes v/2f, so the two differ by
     * exactly a factor of two in every depth. A switch that turned off the
     * tapering and the regularisation but left the paper's convention in place
     * would produce a depth axis twice the patent's while announcing itself as
     * exact -- which is the precise failure this flag exists to prevent.
     *
     * Applied BEFORE rs_tomo_params_check() below, because the check tests the
     * depth grid against the resolution the geometry supports, and that
     * resolution depends on the wavelength. Validating first would test the
     * grid against one convention and then focus with another.
     *
     * Expect it to be noisier. That is the point: it preserves the source's
     * unconditioned processing choices, not what produces the cleanest picture. */
    if (rs_opt_flag(argc, argv, "--patent-exact")) {
        tp.model = RS_TOMO_MODEL_A;
        tp.window = 0;
        tp.remove_y_mean = 0;
        tp.regularisation = 0.0;
        tp.solver = RS_TOMO_SOLVER_LSTSQ;
        tp.y_source = RS_TOMO_Y_SHIFTS;
        tp.convention = RS_WAVELEN_PATENT;
        tp.patent_exact = 1;
        patent_exact = 1;
        printf("patent-chain interpretation (--patent-exact): raw Y (Eq. 21), "
               "%s steering, no taper, true pseudoinverse "
               "(Eq. 24), no detrend, lambda = v/f, --subap paper, "
               "rectangular subaperture filters, no coherence mask, "
               "master-slave pair swept rigidly (Fig. 0.2, Fig. 0.3)\n",
               tp.eq22_literal_t > 0.0
                 ? "EXPERIMENTAL literal exp(j*2*pi*Kz*t*z)"
                 : "conventional exp(j*Kz*z) (rendered Eq. 22 is inconsistent)");
    }

    if (tp.eq22_literal_t > 0.0) {
        fprintf(stderr,
            "warning: --eq22-literal-t enables the experimental literal Eq. 22 "
            "exponent with t=%g. The patent does not define t in this equation, "
            "and this factor is not supported by conventional TomoSAR literature.\n",
            tp.eq22_literal_t);
    }

    /* Switch off what the source does not specify, for a faithful run. */
    if (rs_opt_flag(argc, argv, "--no-window"))   tp.window = 0;
    if (rs_opt_flag(argc, argv, "--keep-mean"))   tp.remove_y_mean = 0;
    tp.regularisation = rs_opt_double(argc, argv, "--regularisation",
                                      tp.regularisation);

    /* An explicit --convention is still parsed for ordinary runs. In
     * --patent-exact it is validated below so an incompatible override fails
     * visibly rather than silently changing the convention. */
    {
        const char *cv = rs_opt(argc, argv, "--convention");
        if (cv) tp.convention = (strcmp(cv, "patent") == 0)
                              ? RS_WAVELEN_PATENT : RS_WAVELEN_PAPER;
    }
    if (patent_exact) {
        if (model && !(*model == 'A' || *model == 'a')) {
            fprintf(stderr,
                "tomo: --patent-exact implements the patent/Giza Eq. 23-24 "
                "steering-matrix inversion; do not combine it with --model B, "
                "C or D.\n");
            rs_cphd_free(&c);
            return 1;
        }
        if (solver && strcmp(solver, "lstsq") != 0) {
            fprintf(stderr,
                "tomo: --patent-exact requires Eq. 24's pseudoinverse; do not "
                "combine it with --solver dft.\n");
            rs_cphd_free(&c);
            return 1;
        }
        if (tp.regularisation != 0.0) {
            fprintf(stderr,
                "tomo: --patent-exact requires Eq. 24's unregularised "
                "pseudoinverse; do not combine it with nonzero "
                "--regularisation.\n");
            rs_cphd_free(&c);
            return 1;
        }
        if (tp.convention != RS_WAVELEN_PATENT) {
            fprintf(stderr,
                "tomo: --patent-exact requires the patent wavelength convention "
                "lambda = v/f; do not combine it with --convention paper.\n");
            rs_cphd_free(&c);
            return 1;
        }
    }

    if ((st = rs_tomo_params_check(&tp)) != RS_OK) {
        rs_report_error("tomo", st);
        rs_cphd_free(&c);
        return 1;
    }

    const size_t size = (size_t)rs_opt_double(argc, argv, "--size", 128);
    rs_grid_t grid = { .n_x = size, .n_y = size,
                       .dx = rs_opt_double(argc, argv, "--grid-cell", 1.0),
                       .dy = rs_opt_double(argc, argv, "--grid-cell", 1.0),
                       .height = 0.0 };
    rs_parse_offset(rs_opt(argc, argv, "--offset"), grid.origin);
    if (rs_resolve_at(rs_opt(argc, argv, "--at"), &c.plane, grid.origin)) {
        rs_cphd_free(&c); return 1;
    }

    /* Depth products default to the blue-low/red-high energy ramp, which is
     * what a reader expects of a tomogram. --palette jet matches the ramp the
     * published figures use, for putting the two side by side. */
    const rs_palette_t palette =
        rs_palette_from_name(rs_opt(argc, argv, "--palette"), RS_PALETTE_ENERGY);

    /* --aperture overrides the along-orbit extent the steering matrix treats as
     * its tomographic baseline.
     *
     * WHY THIS IS EXPOSED. Model A builds baseline i as aperture*(i/k): a
     * uniform ramp derived from a single scalar. A uniform ramp makes the
     * steering matrix a DFT, and in a DFT the sample positions set the scale of
     * the conjugate axis and nothing else. So if this is genuinely a
     * tomographic baseline, changing it must change what is recovered; if it is
     * only a scale factor, changing it must rescale the depth axis and leave
     * every value untouched. That is a decidable question and this option is
     * what makes it decidable from the command line.
     *
     * Genuine multi-baseline tomography does not have this property, because
     * real perpendicular baselines are irregular and their irregularity is the
     * information. See Model C, which takes measured baselines rather than
     * deriving them from one number. */
    if (rs_opt(argc, argv, "--aperture")) {
        /* Already applied by rs_tomo_geometry_overrides(); this only reports it,
         * against the value captured before that call. */
        printf("aperture overridden: %.1f m (collect implies %.1f m)\n",
               tp.aperture, aperture_from_collect);
    }

    {
        /* Which observable becomes the data vector Y. See rs_tomo_y_source_t. */
        const char *y = rs_opt(argc, argv, "--y");
        if (y) {
            if (strcmp(y, "shifts") == 0)   tp.y_source = RS_TOMO_Y_SHIFTS;
            else if (strcmp(y, "los") == 0) tp.y_source = RS_TOMO_Y_LOS;
            else fprintf(stderr, "warning: unknown --y '%s'; keeping the "
                                 "paper's shifts\n", y);
        }
    }
    if (patent_exact && tp.y_source != RS_TOMO_Y_SHIFTS) {
        fprintf(stderr,
            "tomo: --patent-exact requires Eq. 21's complex coregistrator "
            "shift vector; do not combine it with --y los.\n");
        rs_cphd_free(&c);
        return 1;
    }

    /* Orthogonal to --patent-exact and therefore not rejected alongside it: the
     * flag changes how the correlation peak is searched for, not which model,
     * observable, solver or wavelength convention is used. Every choice
     * --patent-exact pins is untouched, so the sidecar's "exact" certification
     * remains correct and gains an [UNOPTIMIZED] tag beside it. */
    const int no_optimize = rs_opt_no_optimize(argc, argv);

    rs_subap_params_t sp;
    rs_subap_params_default(&sp);
    sp.n_looks = (size_t)rs_opt_double(argc, argv, "--n", 64);
    sp.single_thread = no_optimize;
    if (patent_exact) sp.window = 0;

    int ref_mode = rs_parse_reference(argc, argv, &sp);
    if (ref_mode < 0) { rs_cphd_free(&c); return 1; }
    /* The patent's front end is a master-slave pair swept rigidly, and Eq. 21's
     * Y is the offset between them. --patent-exact therefore selects the pair
     * rather than accepting a common master. */
    if (patent_exact) {
        if (rs_opt(argc, argv, "--reference") && ref_mode != RS_MICROM_REF_PAIR) {
            fprintf(stderr,
                "tomo: --patent-exact tracks the master-slave pair the patent "
                "describes; do not combine it with --reference first or "
                "adjacent.\n");
            rs_cphd_free(&c);
            return 1;
        }
        ref_mode = RS_MICROM_REF_PAIR;
        sp.pair = 1;
        fprintf(stderr,
            "warning: --patent-exact selects the master-slave pair the patent\n"
            "         describes. That pairing does NOT recover a frequency on\n"
            "         this project's synthetic single-target fixture -- it\n"
            "         returns the lowest spectral bin whatever is injected. The\n"
            "         run below is faithful to the source, not validated. See\n"
            "         rs_microm_ref_t.\n");
    }

    /* Pulse windows, as in mmotion and sweep. Splitting a focused image
     * spectrally would additionally require the image to have at least twice as
     * many azimuth lines as looks, which couples the grid size to the look
     * count for no physical reason. */
    /* Blocks 2-6 of the patent's scheme are DFT2, two bandpass filters and two
     * IDFT2 -- the spectral route. Fixing the tomographic algebra while leaving
     * the front-end at the pulse-window default would give an exact Eq. 24
     * applied to sub-looks the patent does not describe, so --patent-exact
     * selects the spectral route and rejects incompatible --subap choices. */
    const char *subap_route = rs_opt(argc, argv, "--subap");
    if (!subap_route && patent_exact) subap_route = "paper";
    if (patent_exact) {
        if (subap_route && strcmp(subap_route, "paper") != 0) {
            fprintf(stderr,
                "tomo: --patent-exact requires --subap paper because the patent "
                "front-end is DFT2, bandpass master/slave filters, then IDFT2.\n");
            rs_cphd_free(&c);
            return 1;
        }
    }

    rs_subap_stack_t stack;
    if ((st = rs_build_subaps(&c, &grid, &sp, subap_route, &stack)) != RS_OK) {
        rs_report_error("tomo", st); rs_cphd_free(&c); return 1;
    }

    rs_microm_params_t mp;
    rs_microm_params_default(&mp);
    mp.reference = (rs_microm_ref_t)ref_mode;
    mp.no_optimize = no_optimize;
    {
        /* --estimator selects WHAT is measured, not merely how well. Phase and
         * correlation live in different regimes; see rs_microm_estimator_t. */
        const char *est = rs_opt(argc, argv, "--estimator");
        if (est) {
            if (strcmp(est, "phase") == 0)            mp.estimator = RS_MICROM_EST_PHASE;
            else if (strcmp(est, "correlation") == 0) mp.estimator = RS_MICROM_EST_CORRELATION;
            else if (strcmp(est, "splitband") == 0)   mp.estimator = RS_MICROM_EST_SPLITBAND;
            else fprintf(stderr, "warning: unknown --estimator '%s'; "
                                 "using correlation\n", est);
        }
    }
    if (patent_exact && mp.estimator != RS_MICROM_EST_CORRELATION) {
        fprintf(stderr,
            "tomo: --patent-exact requires the coregistrator shift observable "
            "from --estimator correlation; phase and splitband measure different "
            "quantities.\n");
        rs_subap_stack_free(&stack);
        rs_cphd_free(&c);
        return 1;
    }

    /* Biondi's Y is defined on the coregistrator's shifts -- Eq. 20's {a,b} are
     * "the instantaneous shifts estimated by the coregistrator" -- so the
     * default Y source reads disp_az and disp_rg. The phase estimator does not
     * produce those; it measures pixel phase and writes disp_los. Combining the
     * two silently yields Y = 0 and an identically zero tomogram, which is the
     * same class of failure as substituting a sub-aperture decomposition
     * without saying so. Refuse it and name the alternative. */
    if (mp.estimator == RS_MICROM_EST_PHASE &&
        tp.y_source == RS_TOMO_Y_SHIFTS) {
        fprintf(stderr,
            "tomo: --estimator phase measures pixel phase, not coregistrator\n"
            "      shifts, so the paper's Y = dx + i*dy would be identically\n"
            "      zero. Biondi Eq. 20 defines Y on the shifts, so this pairing\n"
            "      is not the published method under any reading.\n"
            "\n"
            "      Use --y los to build Y from the phase-derived line-of-sight\n"
            "      displacement instead, and record that it is a departure from\n"
            "      the paper. Or keep the default Y and use --estimator\n"
            "      correlation, which is what the paper specifies.\n");
        rs_cphd_free(&c);
        return 1;
    }

    /* The tracking window sets the tomogram's plan-view resolution: one depth
     * profile per window, not per image pixel. A tomogram is therefore far
     * coarser than the radar image it came from, and making it look otherwise
     * by interpolation would be a misrepresentation. Exposed so the trade
     * against tracking reliability can be made deliberately. */
    mp.win_az = mp.win_rg = (size_t)rs_opt_double(argc, argv, "--win", 32);
    mp.stride_az = mp.stride_rg = (size_t)rs_opt_double(argc, argv, "--stride",
                                     (double)(mp.win_az / 2));

    /* The coherence mask has to be settable here, not left at its default.
     *
     * The default of 0.4 suits a simulated scene of bright point targets. Real
     * distributed clutter tracks well below it -- Melbourne windows score about
     * 0.13 -- so the default masks every window, and a masked window contributes
     * nothing to the inversion. The result was an entirely zero tomogram
     * reported without complaint, which is worse than an error because it looks
     * like a measurement of nothing rather than a configuration mistake. */
    mp.coherence_min = patent_exact ? 0.0 : mp.coherence_min;
    mp.coherence_min = rs_opt_double(argc, argv, "--coherence", mp.coherence_min);
    if (patent_exact && mp.coherence_min != 0.0) {
        fprintf(stderr,
            "tomo: --patent-exact feeds block 8 with raw pixel-tracking "
            "vectors; do not combine it with a nonzero --coherence mask.\n");
        rs_subap_stack_free(&stack);
        rs_cphd_free(&c);
        return 1;
    }

    rs_microm_t m;
    if ((st = rs_microm_track(&stack, &mp, &m)) != RS_OK) {
        rs_report_error("tomo", st);
        rs_subap_stack_free(&stack); rs_cphd_free(&c);
        return 1;
    }

    rs_spectrum_t spec;
    memset(&spec, 0, sizeof spec);
    if (!patent_exact) {
        st = rs_spectrum_compute_opts(&m, RS_SPEC_VELOCITY, 0.0,
                                      rs_opt_flag(argc, argv, "--no-detrend")
                                          ? RS_DETREND_NONE : RS_DETREND_LINEAR,
                                      &spec);
        if (st != RS_OK) {
            rs_report_error("tomo", st);
            rs_microm_free(&m); rs_subap_stack_free(&stack);
            rs_cphd_free(&c);
            return 1;
        }
    }

    /* Model C needs genuine perpendicular baselines. Synthesising them here
     * would defeat the purpose of having an uncontested reference model, so it
     * is refused rather than faked. */
    if (tp.model == RS_TOMO_MODEL_C) {
        fprintf(stderr,
            "tomo: model C requires genuine perpendicular baselines from a\n"
            "multi-pass stack, which this single-collect path cannot supply.\n"
            "Synthesising them would make the uncontested reference model as\n"
            "assumption-laden as the model it is meant to check.\n");
        rs_spectrum_free(&spec); rs_microm_free(&m);
        rs_subap_stack_free(&stack); rs_cphd_free(&c);
        return 1;
    }

    rs_tomo_t tomo;
    tp.subap_window = sp.window;
    tp.coherence_min = mp.coherence_min;
    tp.pair_reference = (mp.reference == RS_MICROM_REF_PAIR);
    tp.no_optimize = no_optimize;
    const char *ref_name = (mp.reference == RS_MICROM_REF_PAIR) ? "pair" :
                           (mp.reference == RS_MICROM_REF_ADJACENT) ? "adjacent" :
                           "first";
    snprintf(tp.provenance, sizeof tp.provenance,
             "subap=%s estimator=%s looks=%zu overlap=%.3f win=%zux%zu "
             "coherence=%.2f offset=%.0f,%.0f cell=%.2f fmin=%.2f "
             "detrend=%s ref=%s b_shift=%.12g pair_lag_s=%.12g",
             subap_route ? subap_route : "pulse",
             mp.estimator == RS_MICROM_EST_PHASE ? "phase" :
             mp.estimator == RS_MICROM_EST_SPLITBAND ? "splitband" : "correlation",
             sp.n_looks, sp.overlap, mp.win_az, mp.win_rg, mp.coherence_min,
             grid.origin[0], grid.origin[1], grid.dx,
             rs_opt_double(argc, argv, "--fmin", 0.0),
             (patent_exact || rs_opt_flag(argc, argv, "--no-detrend"))
                 ? "none" : "linear",
             ref_name, stack.b_shift_hz, stack.pair_lag_s);

    st = rs_tomo_focus(&m, &spec, &tp, NULL, &tomo);
    if (st != RS_OK) {
        rs_report_error("tomo", st);
        rs_spectrum_free(&spec); rs_microm_free(&m);
        rs_subap_stack_free(&stack); rs_cphd_free(&c);
        return 1;
    }

    /* An all-zero profile means every window was discarded, not that the scene
     * has no structure. Saying so is the difference between a configuration
     * mistake and an apparent measurement. */
    double tomo_peak = 0.0;
    for (size_t i = 0; i < tomo.n_win * tomo.n_depth; i++) {
        if (tomo.profile[i] > tomo_peak) tomo_peak = tomo.profile[i];
    }
    if (tomo_peak <= 0.0) {
        fprintf(stderr,
            "warning: every depth cell is zero. No window survived the %.2f\n"
            "         coherence mask, so nothing was inverted. Real distributed\n"
            "         clutter often tracks near 0.1; pass --coherence 0 to see an\n"
            "         unmasked result, and treat what comes back accordingly.\n",
            mp.coherence_min);
    }

    /* The alignment null, the depth stage's own null test.
     *
     * Run here, before any product is written, because a tomogram that cannot
     * clear it should not be exported without the number beside it. See
     * rs_tomo_alignment_null() for why a shuffle of the sub-look order is the
     * wrong instrument at this stage: it changes each window's profile as well
     * as their agreement, so it cannot say which half carried the result. */
    {
        const size_t atrials = (size_t)rs_opt_double(argc, argv, "--null-align", 0.0);
        if (atrials > 0) {
            double real = 0.0, nm = 0.0, nsd = 0.0, nmax = 0.0;
            size_t nge = 0;
            const unsigned aseed =
                (unsigned)rs_opt_double(argc, argv, "--null-align-seed", 1.0);
            if (rs_tomo_alignment_null(&tomo, atrials, aseed,
                                       &real, &nm, &nsd, &nmax, &nge) == RS_OK) {
                printf("\nALIGNMENT NULL from %zu re-alignments of the per-window "
                       "depth profiles:\n", atrials);
                printf("  stacked contrast %.2f;  null mean %.2f, sd %.2f, worst %.2f\n",
                       real, nm, nsd, nmax);
                printf("  %.2fx the mean and %.2fx the worst\n",
                       nm > 0.0 ? real / nm : 0.0, nmax > 0.0 ? real / nmax : 0.0);
                printf("  %zu of %zu reached it -- empirical p = %.4f\n",
                       nge, atrials, (double)(nge + 1) / (double)(atrials + 1));
                if (nge > 0) {
                    printf("  RE-ALIGNING THE SAME PROFILES AT RANDOM DEPTHS REACHED\n"
                           "  THIS CONTRAST. The windows do not agree about a depth,\n"
                           "  so the stack is reporting the average shape of their\n"
                           "  own artefacts. This is not a depth measurement.\n");
                } else {
                    printf("  No re-alignment reached it. The windows agree about a\n"
                           "  depth to a degree random alignment does not reproduce --\n"
                           "  which is necessary for a depth claim and not sufficient\n"
                           "  for one: the depth axis is still set by the assumed\n"
                           "  (v, f) above, which nothing here measures.\n");
                }
            } else {
                rs_report_error("tomo", RS_ERR_ARG);
            }
        }
    }

    /* Per-window diagnostic dump.
     *
     * The reported depth for a scene is a peak taken over many windows, and that
     * hides the question of which windows produced it. A depth drawn mostly from
     * windows with no detectable vibration is not a measurement of anything, and
     * it would explain an unstable depth estimate without implicating either the
     * method or this implementation.
     *
     * Emitting prominence, tracking quality, dominant frequency and each
     * window's own peak depth side by side lets that be checked rather than
     * assumed: if the windows that detected something agree on a depth while the
     * rest scatter, the instability is a signal-to-noise problem. */
    const char *win_out = rs_opt(argc, argv, "--windows");
    if (win_out) {
        FILE *wf = fopen(win_out, "w");
        if (wf) {
            fprintf(wf, "window,quality,prominence,dominant_hz,peak_depth_m,peak_value\n");
            for (size_t w = 0; w < tomo.n_win; w++) {
                const double *prof = tomo.profile + w * tomo.n_depth;
                /* Skip the zero bin, exactly as rs_tomo_peak() does. It carries
                 * the profile's DC term and would otherwise dominate; a
                 * diagnostic that used a different convention from the pipeline
                 * would describe itself rather than the pipeline. */
                size_t best = (tomo.n_depth > 1) ? 1 : 0;
                for (size_t k = 2; k < tomo.n_depth; k++) {
                    if (prof[k] > prof[best]) best = k;
                }
                fprintf(wf, "%zu,%.6f,%.6f,%.6f,%.4f,%.6g\n",
                        w,
                        (w < spec.n_win) ? spec.quality[w] : 0.0,
                        (w < spec.n_win) ? spec.prominence[w] : 0.0,
                        (w < spec.n_win) ? spec.dominant_freq[w] : 0.0,
                        tomo.depth[best], prof[best]);
            }
            fclose(wf);
            printf("\nwrote %s (%zu windows)\n", win_out, tomo.n_win);
        }
    }

    /* Block 11 of the patent's computational scheme: geocode the tomogram into
     * a three-dimensional geographic reference system.
     *
     * Everything upstream works in the scene's planar frame, which is enough to
     * process a product and not enough to say where a feature is. Without this
     * a depth cube is a cube of numbers whose position on the earth is implicit
     * in a grid origin recorded elsewhere; with it, every window carries the
     * latitude and longitude it was measured at.
     *
     * Depth is reported as height above the WGS 84 ellipsoid, downward, from
     * the grid plane's own reference height -- so a feature at 12 m depth under
     * a plane at 73 m sits at 61 m. That is stated rather than assumed, because
     * a depth axis and an elevation axis point in opposite directions and
     * silently mixing them is the kind of error that produces plausible
     * coordinates for the wrong place. */
    {
        const char *geo_out = rs_opt(argc, argv, "--geocode");
        if (geo_out) {
            FILE *gf = fopen(geo_out, "w");
            if (!gf) {
                rs_report_error("tomo", RS_ERR_IO);
            } else if (!c.plane.valid) {
                fprintf(stderr, "tomo: the collect carries no scene plane; "
                                "the tomogram cannot be geocoded\n");
                fclose(gf);
            } else {
                fprintf(gf, "window,lat_deg,lon_deg,plane_hae_m,"
                            "peak_depth_m,peak_hae_m,peak_value\n");
                size_t written = 0;
                for (size_t w = 0; w < tomo.n_win; w++) {
                    const double *prof = tomo.profile + w * tomo.n_depth;
                    /* Skip the zero-depth cell, as rs_tomo_peak() and the
                     * --windows dump both do. Every profile has energy there by
                     * construction and it is not a feature, so including it
                     * geolocates the surface rather than the tomogram.
                     *
                     * The exclusion matters most exactly where fidelity to the
                     * patent is highest. Model A's first steering row has
                     * b_perp = 0 and therefore Kz = 0, so its column is all
                     * ones and a constant offset in Y is precisely the z = 0
                     * steering vector. --patent-exact leaves that offset in --
                     * Eq. 21 says nothing about removing a mean -- so the cell
                     * this now skips is the one carrying the registration bias,
                     * and block 11 was reporting its position instead of the
                     * tomogram's. */
                    size_t best = (tomo.n_depth > 1) ? 1 : 0;
                    for (size_t d = 2; d < tomo.n_depth; d++)
                        if (prof[d] > prof[best]) best = d;
                    if (!(prof[best] > 0.0)) continue;

                    /* Window centre in the scene plane, in metres. */
                    const size_t wa = w / m.n_win_rg, wr = w % m.n_win_rg;
                    const double x = grid.origin[0] +
                        ((double)(wa * m.stride_az + m.win_az / 2) -
                         0.5 * (double)(grid.n_x - 1)) * grid.dx;
                    const double y = grid.origin[1] +
                        ((double)(wr * m.stride_rg + m.win_rg / 2) -
                         0.5 * (double)(grid.n_y - 1)) * grid.dy;

                    double lat = 0.0, lon = 0.0, hae = 0.0;
                    if (rs_geo_plane_llh(&c.plane, x, y, &lat, &lon, &hae) != RS_OK)
                        continue;
                    fprintf(gf, "%zu,%.8f,%.8f,%.3f,%.3f,%.3f,%.6g\n",
                            w, lat, lon, hae,
                            tomo.depth[best], hae - tomo.depth[best], prof[best]);
                    written++;
                }
                fclose(gf);
                printf("\nwrote %s (%zu geocoded windows)\n", geo_out, written);
            }
        }
    }

    rs_tomo_write_metadata(&tomo, stdout);

    const char *sec_out = rs_opt(argc, argv, "--section-out");
    if (sec_out) {
        const char *line = rs_opt(argc, argv, "--section");
        double a0 = 0, r0 = 0, a1 = 0, r1 = 0;
        if (!line || sscanf(line, "%lf,%lf,%lf,%lf", &a0, &r0, &a1, &r1) != 4) {
            /* A diagonal across the whole scene, which is what an operator
             * drawing a line freehand tends to produce. */
            a0 = 0; r0 = 0;
            a1 = (double)tomo.n_win_az - 1;
            r1 = (double)tomo.n_win_rg - 1;
        }
        const size_t n_samp = (size_t)rs_opt_double(argc, argv, "--section-width", 512);
        if (rs_tomo_write_section(&tomo, a0, r0, a1, r1, n_samp, sec_out, palette) == RS_OK) {
            printf("\nwrote %s\n", sec_out);
            printf("  section from window (%.0f,%.0f) to (%.0f,%.0f), %zu samples wide,\n"
                   "  %zu depth cells spanning %.1f m\n",
                   a0, r0, a1, r1, n_samp, tomo.n_depth,
                   tomo.depth[tomo.n_depth - 1]);
            printf("  the horizontal axis counts samples, not metres: rendering the same\n"
                   "  section at another --section-width changes every feature's shape.\n");
            printf("  depth cell %.2f m against a resolution of %.2f m\n",
                   tomo.params.depth_cell, tomo.dT);
        } else {
            fprintf(stderr, "tomo: could not write section to %s\n", sec_out);
        }
    }

    const char *slice_out = rs_opt(argc, argv, "--slice-out");
    if (slice_out) {
        const double want = rs_opt_double(argc, argv, "--slice", 0.0);
        double got = 0.0;
        size_t k = 0;
        if (rs_tomo_write_slice(&tomo, want, slice_out, &got, &k, palette) == RS_OK) {
            printf("\nwrote %s\n", slice_out);
            printf("  slice %zu of %zu, labelled %.4f m\n", k, tomo.n_depth, got);
            printf("  that label is lambda_ac = %.5f m applied to slice %zu; the pixels\n"
                   "  are identical at the same z/lambda_ac under any other assumption.\n",
                   tomo.lambda_ac, k);
        } else {
            fprintf(stderr, "tomo: could not write slice to %s\n", slice_out);
        }
    }

    const char *out = rs_opt(argc, argv, "--out");
    if (out) {
        rs_raster_write_cube(tomo.profile, tomo.n_win_az, tomo.n_win_rg, tomo.n_depth,
                             out, no_optimize
                                    ? "window_az, window_rg, depth [UNOPTIMIZED]"
                                    : "window_az, window_rg, depth");
        char meta[512];
        snprintf(meta, sizeof meta, "%s.meta", out);
        FILE *mf = fopen(meta, "w");
        if (mf) { rs_tomo_write_metadata(&tomo, mf); fclose(mf); }
        printf("\nwrote %s, %s.hdr and %s.meta\n", out, out, out);
    }

    rs_tomo_free(&tomo);
    rs_spectrum_free(&spec);
    rs_microm_free(&m);
    rs_subap_stack_free(&stack);
    rs_cphd_free(&c);
    return 0;
}

/* Parse "--offsets x1,y1:x2,y2:..." into a flat array of metre pairs.
 *
 * Writes at most 'cap' pairs and returns how many were parsed. A NULL or empty
 * specification yields the single offset {0,0}, so the option is additive and
 * callers that do not pass it behave exactly as before. */
static size_t rs_parse_offsets(const char *spec, double *out, size_t cap)
{
    if (!spec || !*spec || cap == 0) {
        if (cap > 0) { out[0] = 0.0; out[1] = 0.0; return 1; }
        return 0;
    }
    char buf[512];
    snprintf(buf, sizeof buf, "%s", spec);

    size_t n = 0;
    for (char *tok = strtok(buf, ":;"); tok && n < cap; tok = strtok(NULL, ":;")) {
        char *comma = strchr(tok, ',');
        if (!comma) continue;
        *comma = '\0';
        out[2 * n] = atof(tok);
        out[2 * n + 1] = atof(comma + 1);
        n++;
    }
    if (n == 0) { out[0] = 0.0; out[1] = 0.0; n = 1; }
    return n;
}

/* Vary the two assumed constants and report where the recovered depth goes.
 *
 * This is the experiment that tells a reader how much of a tomogram's depth
 * axis is measurement and how much is the assumption fed in. It belongs beside
 * the null test in any publication of results from this software, and running
 * it costs seconds. */
/* Run the chain over one collect at several grid offsets, appending sweep rows.
 *
 * The offsets are what make the experiment valid. Merging results by wavelength
 * assumes the realisations differ only in scene content, so that a depth which
 * moves must have moved because the assumed constants moved. Three separate
 * collects violate that: they differ in slant range, incidence, dwell and
 * aperture as well, and geometry variation then arrives disguised as wavelength
 * dependence. Chipping one collect at several positions holds every geometric
 * quantity fixed and varies only what is on the ground, which is the real-data
 * counterpart of translating the scatterers in the synthetic fixture.
 *
 * The collect is read and range compressed once and reused across offsets,
 * because that is the expensive step and it does not depend on where the grid
 * sits.
 *
 * Returns the number of rows appended. An offset that cannot be processed is
 * skipped with a message rather than aborting, since the point of running
 * several is to tolerate one going wrong. */
static size_t rs_sweep_one_scene(const char *path, const rs_tomo_params_t *tp,
                                 size_t n_looks, double overlap, size_t grid,
                                 double coherence, double lo, double hi,
                                 size_t steps, size_t rbin_window, double grid_cell,
                                 const double *offsets, size_t n_off, double f_min,
                                 const char *subap_route,
                                 size_t *per_real, size_t *n_real,
                                 rs_tomo_sweep_row_t *rows, size_t cap)
{
    rs_cphd_t c;
    /* No stride: a sweep varies the assumed constants against a fixed
     * measurement, so the pulses it reads must not vary either. */
    if (rs_load_cphd(&c, path, rbin_window, 0, 0) != RS_OK) {
        rs_report_error("sweep", RS_ERR_IO);
        return 0;
    }

    /* Each collect carries its own geometry; a sweep must not judge them all by
     * the first one's, or by the simulator's defaults. */
    rs_tomo_params_t tps = *tp;
    rs_tomo_geometry_from_cphd(&tps, &c);

    const resonarsat_status_t gst = rs_tomo_params_check(&tps);
    if (gst != RS_OK) {
        rs_report_error("sweep", gst);
        rs_cphd_free(&c);
        return 0;
    }

    size_t n_rows = 0;
    for (size_t k = 0; k < n_off && n_rows < cap; k++) {
        rs_grid_t g = { .origin = { offsets[2 * k], offsets[2 * k + 1], 0.0 },
                        .n_x = grid, .n_y = grid,
                        .dx = grid_cell, .dy = grid_cell, .height = 0.0 };

        rs_subap_params_t sp;
        rs_subap_params_default(&sp);
        sp.n_looks = n_looks;
        sp.overlap = overlap;
        /* Carried on the tomo params rather than as another argument to this
         * already long signature. It is the same flag the caller set from
         * --no-optimize, and the sweep tracks once per offset -- outside the
         * (v, f) loop -- so the exhaustive path costs one pass here, not one per
         * row. */
        sp.single_thread = tp->no_optimize;

        rs_subap_stack_t s;
        /* Through the same route selector the other commands use. A sweep
         * that could only run one decomposition was measuring the assumed
         * constants through a stage the paper does not use. */
        if (rs_build_subaps(&c, &g, &sp, subap_route, &s) != RS_OK) {
            rs_report_error("sweep", RS_ERR_ARG);
            continue;
        }

        rs_microm_params_t mp;
        rs_microm_params_default(&mp);
        mp.win_az = mp.win_rg = 32;
        mp.stride_az = mp.stride_rg = 16;
        mp.coherence_min = coherence;
        mp.no_optimize = tp->no_optimize;

        rs_microm_t m;
        rs_spectrum_t spec;
        size_t got = 0;

        if (rs_microm_track(&s, &mp, &m) == RS_OK) {
            /* The sweep varies the assumed constants and nothing else, so it
             * holds the detrend at the default rather than taking it from the
             * command line: a sweep run with one setting and compared against
             * one run with another would attribute the difference to the
             * constants. */
            if (rs_spectrum_compute_band(&m, RS_SPEC_VELOCITY, f_min, &spec) == RS_OK) {
                /* Sized to the request rather than to a fixed 256, and the
                 * shortfalls reported rather than dropped.
                 *
                 * The previous fixed buffer silently produced NOTHING for an
                 * offset whenever steps*steps exceeded it -- so --steps 17 or
                 * more returned an empty sweep with no error, and the summary
                 * counted it as a realisation that simply found nothing. The
                 * sweep is one of this project's two falsification tests; a
                 * quiet zero there is the worst possible failure mode. */
                const size_t want = steps * steps;
                rs_tomo_sweep_row_t *tmp = malloc(want * sizeof *tmp);
                size_t nt = 0;
                if (!tmp) {
                    fprintf(stderr, "sweep: cannot size %zu rows for offset "
                                    "%+.0f,%+.0f\n",
                            want, offsets[2 * k], offsets[2 * k + 1]);
                } else if (rs_tomo_sweep(&m, &spec, &tps, lo, hi, steps,
                                         tmp, &nt) == RS_OK) {
                    for (size_t i = 0; i < nt && n_rows < cap; i++) {
                        rows[n_rows++] = tmp[i];
                        got++;
                    }
                    if (got < nt) {
                        fprintf(stderr, "sweep: TRUNCATED -- offset %+.0f,%+.0f "
                                        "produced %zu rows but only %zu fitted "
                                        "in the result buffer. The fit below "
                                        "omits the rest.\n",
                                offsets[2 * k], offsets[2 * k + 1], nt, got);
                    }
                }
                free(tmp);
                rs_spectrum_free(&spec);
            }
            rs_microm_free(&m);
        }
        printf("  offset %+.0f,%+.0f m: %zu combinations\n",
               offsets[2 * k], offsets[2 * k + 1], got);
        if (per_real && n_real) per_real[(*n_real)++] = got;
        rs_subap_stack_free(&s);
    }

    rs_cphd_free(&c);
    return n_rows;
}

/* Vary the two assumed constants across several scenes and report where the
 * recovered depth goes.
 *
 * This is the experiment that tells a reader how much of a tomogram's depth axis
 * is measurement and how much is the assumption fed in, and it belongs beside the
 * null test in any publication of results from this software.
 *
 * Several scenes are accepted, comma-separated, and it matters that they are. A
 * single realisation cannot separate a real trend from the depth grid's
 * quantisation: the peak lands in whichever cell it lands in, and one such point
 * looks identical to a measurement. With several scenes the depths are averaged
 * per wavelength, each carries a spread, and the fit is additionally reported
 * over those that reproduce. Given one scene the command says so and withholds
 * the filtered fit rather than presenting a number it cannot support. */
static int rs_cmd_sweep(int argc, char **argv)
{
    const char *in = rs_opt(argc, argv, "--cphd");
    if (!in) {
        printf("usage: resonarsat sweep --cphd FILE[,FILE...] --velocity MS --frequency HZ\n"
               "                        [--n N] [--overlap F] [--depth M] [--cell M]\n"
               "                        [--range LO,HI] [--steps N] [--size N]\n"
               "                        [--grid-cell M] [--offsets X,Y:X,Y:...] [--fmin HZ]\n"
               "                        [--model A|B|C|D] [--subap pulse|uniform|paper]\n"
               "                        [--rbins N] [--no-optimize]\n"
               "\n"
               "Re-runs tomographic focusing over a grid of scale factors applied\n"
               "to --velocity and --frequency, and reports where the strongest\n"
               "feature lands in each case. A slope near 1 means the depth axis is\n"
               "the assumed scaling and carries no independent information.\n"
               "\n"
               "SEVERAL REALISATIONS ARE REQUIRED. Prefer --offsets, which chips\n"
               "one collect at several grid positions: geometry stays identical\n"
               "and only the ground content changes, which is what merging by\n"
               "wavelength assumes. Separate collects also differ in range,\n"
               "incidence, dwell and aperture, so geometry variation arrives\n"
               "disguised as wavelength dependence.\n"
               "\n"
               "PASS SEVERAL SCENES, comma-separated. A single scene cannot answer\n"
               "the question: on a discrete depth grid the peak sometimes lands a\n"
               "cell away from where it lands for another scene, and one such point\n"
               "carries the grid's quantisation rather than the data. With several\n"
               "scenes the depths are averaged per wavelength and the fit is also\n"
               "reported over those that reproduce across them.\n");
        return 1;
    }

    rs_tomo_params_t tp;
    rs_tomo_params_default(&tp);
    tp.velocity   = rs_opt_double(argc, argv, "--velocity", 0.0);
    tp.frequency  = rs_opt_double(argc, argv, "--frequency", 0.0);
    tp.depth_max  = rs_opt_double(argc, argv, "--depth", 10.0);
    tp.depth_cell = rs_opt_double(argc, argv, "--cell", 1.0);
    /* rs_sweep_one_scene() reads this off the params to reach the sub-aperture
     * and tracking stages, since it takes no separate flag argument. */
    tp.no_optimize = rs_opt_no_optimize(argc, argv);

    /* The sweep was hardwired to Model A. Model D maps depth from the measured
     * frequency instead of a baseline, and the two respond to the assumed
     * constants in different ways, so which model is being swept has to be the
     * caller's choice rather than an assumption baked in here. */
    const char *sweep_model = rs_opt(argc, argv, "--model");
    if (sweep_model) {
        if (*sweep_model == 'B' || *sweep_model == 'b') tp.model = RS_TOMO_MODEL_B;
        else if (*sweep_model == 'C' || *sweep_model == 'c') tp.model = RS_TOMO_MODEL_C;
        else if (*sweep_model == 'D' || *sweep_model == 'd') tp.model = RS_TOMO_MODEL_D;
        else tp.model = RS_TOMO_MODEL_A;
    }

    /* Only the assumptions can be validated here. The depth grid is judged
     * against the resolution the geometry supports, and the geometry is not
     * known until a collect has been read -- checking it now would test the
     * defaults and then invert with something else. rs_sweep_one_scene() runs
     * the real check once each collect's geometry is in hand. */
    if (tp.velocity <= 0.0 || tp.frequency <= 0.0) {
        rs_report_error("sweep", RS_ERR_ARG);
        return 1;
    }

    double lo = 0.5, hi = 2.0;
    const char *range = rs_opt(argc, argv, "--range");
    if (range) {
        const char *comma = strchr(range, ',');
        if (comma) { lo = atof(range); hi = atof(comma + 1); }
    }
    const size_t steps = (size_t)rs_opt_double(argc, argv, "--steps", 5);

    /* Split the comma-separated scene list. */
    char scenes[1024];
    snprintf(scenes, sizeof scenes, "%s", in);

    const char *paths[16];
    size_t n_scene = 0;
    for (char *tok = strtok(scenes, ","); tok && n_scene < 16; tok = strtok(NULL, ",")) {
        paths[n_scene++] = tok;
    }

    const size_t size = (size_t)rs_opt_double(argc, argv, "--size", 64);
    const size_t n_looks = (size_t)rs_opt_double(argc, argv, "--n", 192);
    const double overlap = rs_opt_double(argc, argv, "--overlap", 0.4);
    const double coherence = rs_opt_double(argc, argv, "--coherence", 0.0);
    const size_t rbin_window = (size_t)rs_opt_double(argc, argv, "--rbins", 0.0);
    const double grid_cell = rs_opt_double(argc, argv, "--grid-cell", 0.5);

    /* The same band floor the micro-motion measurement needs. A drift across
     * the aperture outranks a genuine peak and inflates the background that
     * prominence is measured against, so a sweep run without it is built on
     * contaminated spectra. See rs_spectrum_compute_band(). */
    const double f_min = rs_opt_double(argc, argv, "--fmin", 0.0);

    /* Grid offsets give several realisations from one collect, holding geometry
     * fixed and varying only scene content. See rs_sweep_one_scene(). */
    double offsets[2 * 16];
    const size_t n_off = rs_parse_offsets(rs_opt(argc, argv, "--offsets"),
                                          offsets, 16);

    /* Sized to what was actually asked for: scenes x offsets x steps^2. The
     * old fixed 16*256 silently truncated whenever the request exceeded it, and
     * the summary that follows would then fit a slope to a subset without
     * saying which subset. */
    const size_t cap = n_scene * n_off * steps * steps;
    rs_tomo_sweep_row_t *rows = malloc((cap ? cap : 1) * sizeof *rows);
    if (!rows) {
        fprintf(stderr, "sweep: cannot size %zu rows (%zu scenes x %zu offsets "
                        "x %zu^2 steps)\n", cap, n_scene, n_off, steps);
        return 1;
    }

    size_t n_rows = 0;
    size_t *per_real = malloc((n_scene * n_off + 1) * sizeof *per_real);
    size_t n_real = 0;
    if (!per_real) { free(rows); return 1; }
    for (size_t i = 0; i < n_scene; i++) {
        const size_t got = rs_sweep_one_scene(paths[i], &tp, n_looks, overlap, size,
                                              coherence, lo, hi, steps, rbin_window,
                                              grid_cell, offsets, n_off, f_min,
                                              rs_opt(argc, argv, "--subap"),
                                              per_real, &n_real,
                                              rows + n_rows, cap - n_rows);
        printf("scene %zu (%s): %zu combinations over %zu offset(s)\n",
               i, paths[i], got, n_off);
        n_rows += got;
    }

    if (n_rows == 0) {
        fprintf(stderr, "sweep: no scene produced a usable result\n");
        free(rows);
        return 1;
    }

    /* Merge across scenes so each wavelength carries a spread. */
    rs_tomo_sweep_row_t *merged = malloc(n_rows * sizeof *merged);
    double *sd = malloc(n_rows * sizeof *sd);
    size_t n_merged = 0;
    if (!merged || !sd || rs_tomo_sweep_merge(rows, n_rows, merged, sd, &n_merged) != RS_OK) {
        rs_report_error("sweep", RS_ERR_ARG);
        free(rows); free(merged); free(sd); free(per_real);
        return 1;
    }

    printf("\n%12s %12s %10s %12s %10s\n",
           "lambda_ac", "mean depth", "sd", "depth/lambda", "");
    for (size_t i = 0; i < n_merged; i++) {
        const int reliable = sd[i] <= 0.10 * merged[i].peak_depth;
        printf("%12.5f %12.3f %10.3f %12.1f %10s\n",
               merged[i].lambda_ac, merged[i].peak_depth, sd[i],
               merged[i].peak_depth / merged[i].lambda_ac,
               reliable ? "" : "unstable");
    }

    double slope = 0.0, corr = 0.0;
    if (rs_tomo_sweep_summary(merged, n_merged, &slope, &corr) == RS_OK) {
        printf("\nlog(peak depth) against log(acoustic wavelength):\n");
        printf("  all %zu wavelengths:  slope %+.3f  correlation %+.3f\n",
               n_merged, slope, corr);
    }

    /* Per-realisation slopes.
     *
     * This is the fair form of the test, and the merged form above is not. The
     * merged fit averages depths across realisations at each wavelength and
     * demands they agree, which is right when the realisations are the same
     * scene translated -- as in the synthetic fixture -- but wrong for real
     * data. Two patches of a city 400 m apart hold different buildings, so
     * there is no reason their depths should match even if the method works
     * perfectly, and a spread between them says nothing about the question.
     *
     * What must agree is the RATE. The claim under test is that depth is
     * proportional to the assumed acoustic wavelength; each patch may sit at
     * its own baseline, but the exponent relating the two should be the same
     * everywhere, because that exponent is the whole question. Fitting inside
     * each realisation and comparing slopes tolerates scenes that genuinely
     * differ, which real ones do. */
    if (n_real > 1) {
        printf("\nper-realisation fits (each patch on its own, so patches may\n"
               "differ in content without penalty -- what must agree is the slope):\n");
        double ssum = 0.0, ssum2 = 0.0;
        size_t nfit = 0, base = 0;
        for (size_t r = 0; r < n_real; r++) {
            const size_t cnt = per_real[r];
            double sl = 0.0, co = 0.0;
            if (cnt >= 2 &&
                rs_tomo_sweep_summary(rows + base, cnt, &sl, &co) == RS_OK) {
                printf("  realisation %zu (%zu rows):  slope %+.3f  correlation %+.3f\n",
                       r, cnt, sl, co);
                ssum += sl;
                ssum2 += sl * sl;
                nfit++;
            } else {
                printf("  realisation %zu (%zu rows):  no fit\n", r, cnt);
            }
            base += cnt;
        }
        if (nfit > 1) {
            const double slope_mean = ssum / (double)nfit;
            const double slope_var = ssum2 / (double)nfit - slope_mean * slope_mean;
            const double slope_sd = slope_var > 0.0 ? sqrt(slope_var) : 0.0;
            printf("  mean slope %+.3f, spread %.3f across %zu realisations\n",
                   slope_mean, slope_sd, nfit);
            if (slope_sd < 0.15 && slope_mean > 0.85 && slope_mean < 1.15) {
                printf("\n  Every patch agrees that depth is proportional to the ASSUMED\n"
                       "  wavelength. The depth axis carries no independent information.\n");
            } else if (slope_sd < 0.15 && fabs(slope_mean) < 0.2) {
                printf("\n  Every patch agrees that depth does NOT follow the assumption.\n"
                       "  Something in the data pins it. Worth pursuing hard.\n");
            } else if (slope_sd < 0.15) {
                /* Agreement on a value that is neither 1 nor 0 is a different
                 * situation from disagreement, and saying "the patches
                 * disagree" when their spread is small would misdescribe the
                 * data. The usual cause is a depth grid too coarse for the
                 * peak to move smoothly, which flattens every fit by the same
                 * amount and so preserves agreement while destroying the
                 * value. */
                printf("\n  The patches agree closely (spread %.3f) on a slope of %+.3f,\n"
                       "  which is neither the 1 that means the depth axis is the\n"
                       "  assumption nor the 0 that means something pins it. Ambiguous.\n"
                       "  Check whether the depth cell is coarse enough to be quantising\n"
                       "  the peak before reading anything into the value.\n",
                       slope_sd, slope_mean);
            } else {
                printf("\n  The patches disagree about the slope (spread %.3f). The depth\n"
                       "  estimate is not stable enough for this test to decide anything,\n"
                       "  and no conclusion should be drawn in either direction.\n", slope_sd);
            }
        }
    }

    /* Model D's depth comes from the measured frequency scaled by the assumed
     * velocity, so lambda_ac is the wrong abscissa for it: the fit above
     * regresses against v/f when only v matters. That does not give noise, it
     * gives a clean wrong answer -- analytically +0.5, which measured as +0.496
     * -- and a reader could easily take it for a real intermediate slope.
     * Regressing against the velocity is the equivalent test for this model. */
    if (tp.model == RS_TOMO_MODEL_D) {
        size_t n = 0;
        double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
        for (size_t i = 0; i < n_rows; i++) {
            if (rows[i].peak_depth <= 0.0 || rows[i].velocity <= 0.0) continue;
            const double x = log(rows[i].velocity), y = log(rows[i].peak_depth);
            sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
            n++;
        }
        if (n > 2) {
            const double dn = (double)n;
            const double cov = sxy / dn - (sx / dn) * (sy / dn);
            const double vx = sxx / dn - (sx / dn) * (sx / dn);
            const double vy = syy / dn - (sy / dn) * (sy / dn);
            if (vx > 0.0 && vy > 0.0) {
                printf("\nMODEL D: depth against the ASSUMED VELOCITY (its own variable):\n");
                printf("  %zu points: slope %+.3f  correlation %+.3f\n",
                       n, cov / vx, cov / sqrt(vx * vy));
                printf("  The fit above, against lambda_ac, is not meaningful for this\n"
                       "  model: depth here does not use the assumed frequency at all.\n");
            }
        }
    }

    size_t n_used = 0;
    if (n_scene * n_off > 1 &&
        rs_tomo_sweep_summary_reliable(merged, sd, n_merged, 0.10,
                                       &slope, &corr, &n_used) == RS_OK) {
        printf("  %zu reproducible:      slope %+.3f  correlation %+.3f\n\n",
               n_used, slope, corr);
        if (corr > 0.95 && slope > 0.85 && slope < 1.15) {
            printf("  Depth is proportional to the ASSUMED acoustic wavelength. The\n"
                   "  depth axis carries no independent depth information: it is the\n"
                   "  operator's (v, f) applied to a fixed spectral feature.\n");
        } else if (fabs(slope) < 0.2 && fabs(corr) > 0.8) {
            printf("  Features hold position while the assumptions vary. Something in\n"
                   "  the data pins them. Worth pursuing.\n");
        } else {
            printf("  Neither a clean slope of 1 nor of 0. Ambiguous, and should be\n"
                   "  reported as such rather than rounded toward a preferred answer.\n");
        }
    } else if (n_scene * n_off == 1) {
        printf("\n  ONE REALISATION ONLY. A single realisation cannot separate a real\n"
               "  trend from the depth grid's quantisation -- pass --offsets, or\n"
               "  several comma-separated scenes, to get a reproducibility-filtered\n"
               "  fit.\n");
    }

    free(merged);
    free(sd);
    free(rows);
    free(per_real);
    return 0;
}

/* Dispatch to a subcommand. */
int main(int argc, char **argv)
{
    if (argc < 2) { rs_usage(); return 1; }

    const char *cmd = argv[1];
    if (strcmp(cmd, "feasibility") == 0) return rs_cmd_feasibility(argc - 1, argv + 1);
    if (strcmp(cmd, "info") == 0)        return rs_cmd_info(argc - 1, argv + 1);
    if (strcmp(cmd, "focus") == 0)       return rs_cmd_focus(argc - 1, argv + 1);
    if (strcmp(cmd, "mmotion") == 0)     return rs_cmd_mmotion(argc - 1, argv + 1);
    if (strcmp(cmd, "tomo") == 0)        return rs_cmd_tomo(argc - 1, argv + 1);
    if (strcmp(cmd, "sweep") == 0)       return rs_cmd_sweep(argc - 1, argv + 1);

    if (strcmp(cmd, "-h") == 0 || strcmp(cmd, "--help") == 0) { rs_usage(); return 0; }

    fprintf(stderr, "unknown command '%s'\n\n", cmd);
    rs_usage();
    return 1;
}
