/* Image formation from phase history by time-domain backprojection.
 * See include/resonarsat/focus.h for the contract and the rationale. */

#include "resonarsat/focus.h"
#include "resonarsat/geom.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define RS_CPHD_MAGIC   0x52534348u  /* "RSCH" */
#define RS_CPHD_VERSION 2u

/* Allocate the three arrays a phase-history container owns. */
resonarsat_status_t rs_cphd_alloc(rs_cphd_t *cphd, size_t n_pulse, size_t n_rbin)
{
    if (!cphd) return RS_ERR_ARG;
    if (n_pulse == 0 || n_rbin == 0) {
        rs_set_error("cphd: dimensions must be positive (got %zu pulses, %zu bins)",
                     n_pulse, n_rbin);
        return RS_ERR_ARG;
    }
    if (n_pulse > SIZE_MAX / n_rbin) {
        rs_set_error("cphd: %zu pulses by %zu bins overflows", n_pulse, n_rbin);
        return RS_ERR_RANGE;
    }

    memset(cphd, 0, sizeof *cphd);
    cphd->signal = calloc(n_pulse * n_rbin, sizeof *cphd->signal);
    cphd->pos    = calloc(3 * n_pulse, sizeof *cphd->pos);
    cphd->t      = calloc(n_pulse, sizeof *cphd->t);
    cphd->r_ref  = calloc(n_pulse, sizeof *cphd->r_ref);
    if (!cphd->signal || !cphd->pos || !cphd->t || !cphd->r_ref) {
        rs_cphd_free(cphd);
        rs_set_error("cphd: cannot allocate %zu pulses of %zu bins", n_pulse, n_rbin);
        return RS_ERR_ALLOC;
    }
    cphd->n_pulse = n_pulse;
    cphd->n_rbin  = n_rbin;
    return RS_OK;
}

/* Release every owned buffer. Safe on a partially constructed container. */
void rs_cphd_free(rs_cphd_t *cphd)
{
    if (!cphd) return;
    free(cphd->signal);
    free(cphd->pos);
    free(cphd->t);
    free(cphd->r_ref);
    cphd->signal = NULL;
    cphd->pos = NULL;
    cphd->t = NULL;
    cphd->r_ref = NULL;
    cphd->n_pulse = cphd->n_rbin = 0;
}

/* Backproject one pulse window onto the ground grid.
 *
 * The inner loop is the pipeline's hot spot. It is written to keep the per-cell
 * work minimal: the grid position is computed once per cell, and per pulse the
 * only operations are a vector difference, a square root, a linear interpolation
 * and a complex multiply-accumulate.
 *
 * Note the sign of the phase term. The range-compressed sample at delay R
 * carries exp(-j*4*pi*R/lambda) from propagation, so the conjugate term is
 * applied here to bring every pulse's contribution into phase at the correct
 * cell. Getting this sign wrong does not produce an obviously broken image --
 * it produces a defocused one, which is far harder to notice, so the simulator
 * test in tests/test_focus.c pins it. */
resonarsat_status_t rs_focus_backproject(const rs_cphd_t *cphd,
                                         const rs_grid_t *grid,
                                         size_t pulse_start,
                                         size_t pulse_count,
                                         rs_slc_t *img)
{
    return rs_focus_backproject_opts(cphd, grid, pulse_start, pulse_count, NULL, img);
}

/* Backproject one pulse window onto the ground grid, with execution options.
 * See rs_focus_backproject() for the arithmetic and rs_focus_opts_t for the
 * options -- which change scheduling only, never the result. */
resonarsat_status_t rs_focus_backproject_opts(const rs_cphd_t *cphd,
                                              const rs_grid_t *grid,
                                              size_t pulse_start,
                                              size_t pulse_count,
                                              const rs_focus_opts_t *opts,
                                              rs_slc_t *img)
{
    if (!cphd || !grid || !img || !cphd->signal || !img->data) return RS_ERR_ARG;

    if (pulse_start >= cphd->n_pulse || pulse_count == 0 ||
        pulse_start + pulse_count > cphd->n_pulse) {
        rs_set_error("focus: pulse window [%zu, %zu) outside available %zu pulses",
                     pulse_start, pulse_start + pulse_count, cphd->n_pulse);
        return RS_ERR_ARG;
    }
    /* The image's azimuth axis is the ALONG-TRACK direction, which is the
     * grid's x. Mapping grid y to image azimuth instead would leave every
     * downstream stage that splits "azimuth" spectra operating on the
     * cross-track axis -- silently, and with plausible-looking output. */
    if (img->n_az != grid->n_x || img->n_rg != grid->n_y) {
        rs_set_error("focus: image %zux%zu (az x rg) does not match grid "
                     "%zu along-track x %zu cross-track",
                     img->n_az, img->n_rg, grid->n_x, grid->n_y);
        return RS_ERR_ARG;
    }
    if (cphd->lambda <= 0.0 || cphd->dr <= 0.0) {
        rs_set_error("focus: phase history lacks wavelength or range spacing");
        return RS_ERR_MISSING_META;
    }

    const double k_phase = 4.0 * M_PI / cphd->lambda;
    const size_t n_x = grid->n_x, n_y = grid->n_y, n_rbin = cphd->n_rbin;
    const double dr = cphd->dr;

    const long total = (long)(n_x * n_y);

    /* One loop, one body, two schedules.
     *
     * The serial mode is expressed as OpenMP's 'if' clause rather than as a
     * second copy of the inner loop. A duplicated body is the standard way to
     * write an "unoptimized reference" and it is the wrong way here: the two
     * copies would drift, and a reference implementation that has diverged from
     * the thing it is checking reports the divergence as a finding. With one body
     * the arithmetic is the same expression tree by construction, so the test
     * that asserts the two modes agree bitwise is testing the scheduler, which is
     * the only thing that differs. */
    const int single_thread = opts && opts->single_thread;

    /* Interpolator width. Anything below 4 is the historical two-tap linear
     * path; see rs_focus_opts_t for why the wider kernel exists. */
    int taps = opts ? opts->range_taps : 0;
    if (taps < 4) taps = 0;
    if (taps > 16) taps = 16;
    const int half_taps = taps / 2;

    /* Kernel coefficients, precomputed on a fine grid of sub-bin positions.
     *
     * Evaluating sin() twice per tap per pulse per cell made the wider kernel
     * three to four times the cost of the linear one, which put the measurement
     * it exists for out of reach on a real collect. The kernel depends only on
     * the fractional bin position, so it is tabulated once here and looked up.
     *
     * RS_TAP_STEPS positions quantise that fraction to 1/2048 of a bin. The
     * error that introduces is bounded by the kernel's slope, of order pi, so
     * about -56 dB on the interpolated sample -- three orders below the 1.4 dB
     * the kernel is correcting. It does NOT touch the phase term, which is
     * computed from the exact range independently of this table. */
    double *ktab = NULL;
    const int RS_TAP_STEPS = 2048;
    if (taps > 0) {
        ktab = malloc((size_t)(RS_TAP_STEPS + 1) * (size_t)taps * sizeof *ktab);
        if (!ktab) return RS_ERR_ALLOC;
        for (int q = 0; q <= RS_TAP_STEPS; q++) {
            const double f = (double)q / (double)RS_TAP_STEPS;
            for (int t = -half_taps + 1; t <= half_taps; t++) {
                const double x = f - (double)t;
                double w;
                if (fabs(x) < 1e-12) {
                    w = 1.0;
                } else {
                    const double sx = sin(M_PI * x) / (M_PI * x);
                    w = sx * 0.5 * (1.0 + cos(M_PI * x / (double)half_taps));
                }
                ktab[(size_t)q * (size_t)taps + (size_t)(t + half_taps - 1)] = w;
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (!single_thread)
#else
    (void)single_thread;
#endif
    for (long cell = 0; cell < total; cell++) {
        /* Row index runs along track (azimuth), column across it (range). */
        const size_t ix = (size_t)cell / n_y;
        const size_t iy = (size_t)cell % n_y;

        const double gx = grid->origin[0] + ((double)ix - 0.5 * (double)(n_x - 1)) * grid->dx;
        const double gy = grid->origin[1] + ((double)iy - 0.5 * (double)(n_y - 1)) * grid->dy;
        const double gz = grid->origin[2] + grid->height;

        double acc_re = 0.0, acc_im = 0.0;

        for (size_t p = pulse_start; p < pulse_start + pulse_count; p++) {
            const double px = cphd->pos[3 * p + 0];
            const double py = cphd->pos[3 * p + 1];
            const double pz = cphd->pos[3 * p + 2];

            const double ddx = px - gx, ddy = py - gy, ddz = pz - gz;
            const double R = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);

            /* Locate the sample in this pulse at two-way delay R and
             * interpolate linearly between neighbouring bins. Ranges are
             * relative to this pulse's reference range, because the receive
             * window follows the scene (see rs_cphd_t.r_ref). Cells whose range
             * falls outside the collected swath contribute nothing. */
            const double fbin = (R - cphd->r_ref[p]) / dr + 0.5 * (double)n_rbin;
            if (fbin < 0.0 || fbin >= (double)(n_rbin - 1)) continue;

            const size_t b = (size_t)fbin;
            const double frac = fbin - (double)b;

            const float complex *row = cphd->signal + p * n_rbin;
            float complex s;
            if (taps == 0 ||
                b < (size_t)half_taps || b + (size_t)half_taps >= n_rbin) {
                /* Linear, and also the fallback within half a kernel of the
                 * swath edge, where the wider kernel would read outside the
                 * collected data. Truncating the kernel there instead would
                 * taper the edge of every image by a different amount depending
                 * on the option, which is worse than a documented seam. */
                s = row[b] * (float)(1.0 - frac) + row[b + 1] * (float)frac;
            } else {
                /* Hann-windowed sinc, from the table above. The window is on the
                 * kernel's own support so the taps fall to zero at the ends
                 * rather than being cut, which keeps the stopband from ringing. */
                const double *w = ktab + (size_t)(int)(frac * (double)RS_TAP_STEPS + 0.5)
                                         * (size_t)taps;
                const float complex *base = row + (size_t)((long)b - half_taps + 1);
                double acc_r = 0.0, acc_i = 0.0;
                for (int t = 0; t < taps; t++) {
                    acc_r += w[t] * (double)crealf(base[t]);
                    acc_i += w[t] * (double)cimagf(base[t]);
                }
                s = (float)acc_r + (float)acc_i * I;
            }

            /* Undo the propagation phase, against whichever reference the data
             * were recorded on (see rs_cphd_t.phase_ref_srp). */
            const double ph = k_phase * (cphd->phase_ref_srp ? (R - cphd->r_ref[p]) : R);
            const double cr = cos(ph), ci = sin(ph);
            const double sr = (double)crealf(s), si = (double)cimagf(s);
            acc_re += sr * cr - si * ci;
            acc_im += sr * ci + si * cr;
        }

        img->data[(size_t)cell] = (float)acc_re + (float)acc_im * I;
    }

    free(ktab);

    /* Fill in the geometry the focused image inherits from the collection.
     *
     * Platform speed, slant range and incidence angle are measured from the
     * recorded pulse positions rather than assumed. Downstream stages need all
     * three to convert a tracked pixel shift into a line-of-sight velocity, and
     * a hardcoded fallback there would silently scale every velocity the
     * pipeline reports. */
    img->fc = cphd->fc;
    img->lambda = cphd->lambda;
    img->az_spacing_m = grid->dx;
    img->rg_spacing_m = grid->dy;
    img->pulse_prf = cphd->prf;

    {
        const size_t p_mid = pulse_start + pulse_count / 2;
        const size_t p_end = pulse_start + pulse_count - 1;

        if (pulse_count >= 2) {
            const double dt_win = cphd->t[p_end] - cphd->t[pulse_start];
            const double dx = cphd->pos[3 * p_end + 0] - cphd->pos[3 * pulse_start + 0];
            const double dy = cphd->pos[3 * p_end + 1] - cphd->pos[3 * pulse_start + 1];
            const double dz = cphd->pos[3 * p_end + 2] - cphd->pos[3 * pulse_start + 2];
            const double flown = sqrt(dx * dx + dy * dy + dz * dz);
            if (dt_win > 0.0) img->v_platform = flown / dt_win;
        }

        img->r0 = (cphd->r_ref && cphd->r_ref[p_mid] > 0.0) ? cphd->r_ref[p_mid]
                                                            : cphd->r_near;

        /* Incidence at the scene reference point: the angle between the
         * line of sight and the local vertical, from the platform height. */
        const double pz = cphd->pos[3 * p_mid + 2] - grid->origin[2];
        if (img->r0 > 0.0 && fabs(pz) <= img->r0) {
            img->incidence = acos(fabs(pz) / img->r0);
        }
    }

    if (cphd->prf > 0.0) {
        /* One "azimuth line" of the focused image corresponds to the dwell
         * divided by the number of grid rows; but for a backprojected image the
         * meaningful timing quantity is the window duration, which is what the
         * sub-aperture stage needs. Record the window's duration and derive a
         * nominal line interval from the grid so downstream validation has
         * something self-consistent to check. */
        const double t_first = cphd->t[pulse_start];
        const double t_last  = cphd->t[pulse_start + pulse_count - 1];
        img->t0 = t_first;
        img->t_dwell = t_last - t_first;
        img->azimuth_time_interval = (n_x > 1 && img->t_dwell > 0.0)
                                   ? img->t_dwell / (double)(n_x - 1)
                                   : 1.0 / cphd->prf;
        img->fs_az = 1.0 / img->azimuth_time_interval;
    }

    snprintf(img->source, sizeof img->source, "backprojected/%s",
             cphd->source[0] ? cphd->source : "cphd");

    return RS_OK;
}

/* Allocate to the grid and focus the whole aperture. */
resonarsat_status_t rs_focus_full(const rs_cphd_t *cphd, const rs_grid_t *grid, rs_slc_t *img)
{
    if (!cphd || !grid || !img) return RS_ERR_ARG;
    resonarsat_status_t st = rs_slc_alloc(img, grid->n_x, grid->n_y);
    if (st != RS_OK) return st;
    st = rs_focus_backproject(cphd, grid, 0, cphd->n_pulse, img);
    if (st != RS_OK) rs_slc_free(img);
    return st;
}

/* On-disk header for the interchange format. Fixed-width fields so the file is
 * portable between builds; the magic and version make a stale file fail loudly
 * rather than being misinterpreted as valid geometry. */
typedef struct {
    uint32_t magic, version;
    uint64_t n_pulse, n_rbin;
    double   r_near, dr, fc, lambda, prf;
} rs_cphd_header_t;

/* Serialise a phase-history container. */
resonarsat_status_t rs_cphd_write(const rs_cphd_t *cphd, const char *path)
{
    if (!cphd || !path || !cphd->signal) return RS_ERR_ARG;

    FILE *f = fopen(path, "wb");
    if (!f) {
        rs_set_error("cphd write: cannot open %s", path);
        return RS_ERR_IO;
    }

    rs_cphd_header_t h = {
        .magic = RS_CPHD_MAGIC, .version = RS_CPHD_VERSION,
        .n_pulse = cphd->n_pulse, .n_rbin = cphd->n_rbin,
        .r_near = cphd->r_near, .dr = cphd->dr, .fc = cphd->fc,
        .lambda = cphd->lambda, .prf = cphd->prf
    };

    int ok = fwrite(&h, sizeof h, 1, f) == 1
          && fwrite(cphd->t, sizeof *cphd->t, cphd->n_pulse, f) == cphd->n_pulse
          && fwrite(cphd->r_ref, sizeof *cphd->r_ref, cphd->n_pulse, f) == cphd->n_pulse
          && fwrite(cphd->pos, sizeof *cphd->pos, 3 * cphd->n_pulse, f) == 3 * cphd->n_pulse
          && fwrite(cphd->signal, sizeof *cphd->signal,
                    cphd->n_pulse * cphd->n_rbin, f) == cphd->n_pulse * cphd->n_rbin;
    fclose(f);

    if (!ok) {
        rs_set_error("cphd write: short write to %s", path);
        return RS_ERR_IO;
    }
    return RS_OK;
}

/* Read a phase-history container back.
 *
 * Every field that determines an allocation size is validated before it is
 * used. A truncated or corrupt file must produce a described failure, never a
 * short buffer that the focusing loop then reads past. */
resonarsat_status_t rs_cphd_read(rs_cphd_t *cphd, const char *path)
{
    if (!cphd || !path) return RS_ERR_ARG;

    /* Zero the output before anything can fail.
     *
     * A caller that checks the status and frees the struct on the error path is
     * doing the ordinary thing, and would otherwise be freeing whatever was on
     * its stack. This is the same defect that rs_tomo_focus() carried (bug 19);
     * it survived here because no test had ever made this function fail. */
    memset(cphd, 0, sizeof *cphd);

    FILE *f = fopen(path, "rb");
    if (!f) {
        rs_set_error("cphd read: cannot open %s", path);
        return RS_ERR_IO;
    }

    rs_cphd_header_t h;
    if (fread(&h, sizeof h, 1, f) != 1) {
        fclose(f);
        rs_set_error("cphd read: %s is too short to contain a header", path);
        return RS_ERR_FORMAT;
    }
    if (h.magic != RS_CPHD_MAGIC) {
        fclose(f);
        rs_set_error("cphd read: %s is not a phase-history file (magic 0x%08x)",
                     path, h.magic);
        return RS_ERR_FORMAT;
    }
    if (h.version != RS_CPHD_VERSION) {
        fclose(f);
        rs_set_error("cphd read: %s has version %u, this build expects %u",
                     path, h.version, RS_CPHD_VERSION);
        return RS_ERR_UNSUPPORTED;
    }
    /* Bound the dimensions before they reach an allocator. */
    if (h.n_pulse == 0 || h.n_rbin == 0 ||
        h.n_pulse > (1u << 24) || h.n_rbin > (1u << 24)) {
        fclose(f);
        rs_set_error("cphd read: %s declares implausible dimensions %llu x %llu",
                     path, (unsigned long long)h.n_pulse, (unsigned long long)h.n_rbin);
        return RS_ERR_FORMAT;
    }

    resonarsat_status_t st = rs_cphd_alloc(cphd, (size_t)h.n_pulse, (size_t)h.n_rbin);
    if (st != RS_OK) { fclose(f); return st; }

    const size_t np = cphd->n_pulse, nb = cphd->n_rbin;
    int ok = fread(cphd->t, sizeof *cphd->t, np, f) == np
          && fread(cphd->r_ref, sizeof *cphd->r_ref, np, f) == np
          && fread(cphd->pos, sizeof *cphd->pos, 3 * np, f) == 3 * np
          && fread(cphd->signal, sizeof *cphd->signal, np * nb, f) == np * nb;
    fclose(f);

    if (!ok) {
        rs_cphd_free(cphd);
        rs_set_error("cphd read: %s is truncated", path);
        return RS_ERR_FORMAT;
    }

    cphd->r_near = h.r_near;
    cphd->dr     = h.dr;
    cphd->fc     = h.fc;
    cphd->lambda = h.lambda;
    cphd->prf    = h.prf;
    snprintf(cphd->source, sizeof cphd->source, "%s", "cphd-file");
    return RS_OK;
}

/* Azimuth resolution achievable from a given number of pulses.
 *
 * Platform speed is measured from the recorded positions rather than assumed,
 * so this stays correct for any collect geometry the container describes. The
 * aperture is the distance actually flown during the window. */
double rs_focus_azimuth_resolution(const rs_cphd_t *cphd, size_t pulse_count)
{
    if (!cphd || !cphd->pos || cphd->n_pulse < 2 || pulse_count < 2) return 0.0;
    if (pulse_count > cphd->n_pulse) pulse_count = cphd->n_pulse;
    if (cphd->lambda <= 0.0) return 0.0;

    /* Distance flown over the window, from first to last recorded position. */
    const size_t last = pulse_count - 1;
    const double dx = cphd->pos[3 * last + 0] - cphd->pos[0];
    const double dy = cphd->pos[3 * last + 1] - cphd->pos[1];
    const double dz = cphd->pos[3 * last + 2] - cphd->pos[2];
    const double aperture = sqrt(dx * dx + dy * dy + dz * dz);
    if (aperture <= 0.0) return 0.0;

    const double r_centre = (cphd->r_ref && cphd->r_ref[pulse_count / 2] > 0.0)
                          ? cphd->r_ref[pulse_count / 2]
                          : cphd->r_near;
    if (r_centre <= 0.0) return 0.0;

    return cphd->lambda * r_centre / (2.0 * aperture);
}
