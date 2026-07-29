/* Synthetic phase-history generator -- the only place ground truth exists.
 *
 * Emits range-compressed phase history for a straight-line spotlight collect
 * over a set of point scatterers, some of which vibrate sinusoidally. Because
 * the injected positions, frequencies and amplitudes are known exactly, this is
 * what every stage downstream is checked against.
 *
 * The generated geometry is deliberately simple: a platform flying along +x at
 * constant height and speed, staring at the scene centre. That is enough to
 * exercise focusing, sub-aperture decomposition, tracking and spectral
 * estimation, and it keeps the ground truth analytic. */

#include "resonarsat/focus.h"
#include "resonarsat/geom.h"

#include <math.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One simulated point scatterer. A stationary target has amplitude zero on its
 * vibration term; a vibrating one displaces sinusoidally along the vertical,
 * which projects onto the line of sight through the incidence angle. */
typedef struct {
    double x, y, z;      /* position in the local frame, m */
    double rcs;          /* reflectivity, arbitrary linear units */
    double vib_freq;     /* Hz, 0 for a static target */
    double vib_amp;      /* m, vertical displacement amplitude */
    double vib_phase;    /* rad, phase at t = 0 */
} rs_sim_target_t;

/* Simulation parameters, with defaults matching a free X-band spotlight
 * collect of the kind this project targets. */
typedef struct {
    double fc;            /* carrier, Hz */
    double prf;           /* pulses per second */
    double t_dwell;       /* dwell, s */
    double v_platform;    /* m/s */
    double height;        /* platform height above the scene plane, m */
    double range_offset;  /* cross-track offset of the track from the scene, m */
    size_t n_rbin;        /* range bins per pulse */
    double dr;            /* range bin spacing, m */
    double range_res;     /* range resolution, m (sets the compressed pulse width) */
} rs_sim_params_t;

/* Fill simulation parameters with a plausible X-band spotlight geometry. */
static void rs_sim_defaults(rs_sim_params_t *p)
{
    /* These match the configuration the project's published measurements were
     * taken at (see docs/FINDINGS.md), so that the documented reproduction
     * commands reproduce them. The dwell in particular is not arbitrary: the
     * sub-look count needed to satisfy the phase-ambiguity condition scales with
     * it, and a shorter dwell at the same look count gives coarser sub-looks and
     * a different answer. */
    p->fc = 9.6e9;
    p->prf = 400.0;
    p->t_dwell = 20.0;
    p->v_platform = 7500.0;
    p->height = 500000.0;
    p->range_offset = 350000.0;
    p->n_rbin = 256;
    p->dr = 0.5;
    p->range_res = 1.0;
}

/* Generate range-compressed phase history for a target list.
 *
 * For each pulse the platform position is advanced along the track, each
 * target's instantaneous position is evaluated (including its vibration), and
 * a compressed pulse response is deposited at the corresponding range bin.
 *
 * The compressed response is modelled as a Gaussian envelope of width
 * range_res carrying the propagation phase exp(-j*4*pi*R/lambda). A real system
 * would give a sinc; a Gaussian avoids sidelobes that would complicate the
 * ground-truth checks without changing anything the tests measure. The phase,
 * which is what the whole pipeline actually reads, is exact.
 *
 * The near range is chosen so the scene centre falls in the middle of the
 * collected swath. Returns RS_ERR_ALLOC on failure. */
static resonarsat_status_t rs_sim_generate(const rs_sim_params_t *p,
                                           const rs_sim_target_t *targets,
                                           size_t n_target,
                                           rs_cphd_t *cphd)
{
    const size_t n_pulse = (size_t)(p->prf * p->t_dwell);
    resonarsat_status_t st = rs_cphd_alloc(cphd, n_pulse, p->n_rbin);
    if (st != RS_OK) return st;

    cphd->fc = p->fc;
    cphd->lambda = RS_C_LIGHT / p->fc;
    cphd->prf = p->prf;
    cphd->dr = p->dr;
    snprintf(cphd->source, sizeof cphd->source, "simulated");

    /* Range to the scene centre at closest approach, used to centre the swath. */
    const double r_centre = sqrt(p->height * p->height + p->range_offset * p->range_offset);
    cphd->r_near = r_centre - 0.5 * (double)p->n_rbin * p->dr;

    const double k_phase = 4.0 * M_PI / cphd->lambda;
    const double sigma = p->range_res / 2.355;   /* FWHM to Gaussian sigma */

    for (size_t i = 0; i < n_pulse; i++) {
        const double t = (double)i / p->prf - 0.5 * p->t_dwell;
        cphd->t[i] = t + 0.5 * p->t_dwell;

        /* Platform flies along +x, offset cross-track in y, at fixed height. */
        cphd->pos[3 * i + 0] = p->v_platform * t;
        cphd->pos[3 * i + 1] = p->range_offset;
        cphd->pos[3 * i + 2] = p->height;

        /* Motion compensation: the receive window follows the scene reference
         * point (the origin), so bin n_rbin/2 always sits on it. Without this
         * a long dwell would walk the target clean out of a narrow swath. */
        cphd->r_ref[i] = sqrt(cphd->pos[3 * i + 0] * cphd->pos[3 * i + 0]
                            + cphd->pos[3 * i + 1] * cphd->pos[3 * i + 1]
                            + cphd->pos[3 * i + 2] * cphd->pos[3 * i + 2]);

        float complex *row = cphd->signal + i * p->n_rbin;

        for (size_t g = 0; g < n_target; g++) {
            const rs_sim_target_t *tg = &targets[g];

            /* Vibration displaces the target vertically. */
            double dz = 0.0;
            if (tg->vib_freq > 0.0 && tg->vib_amp != 0.0) {
                dz = tg->vib_amp * sin(2.0 * M_PI * tg->vib_freq * cphd->t[i] + tg->vib_phase);
            }

            const double dx = cphd->pos[3 * i + 0] - tg->x;
            const double dy = cphd->pos[3 * i + 1] - tg->y;
            const double dzz = cphd->pos[3 * i + 2] - (tg->z + dz);
            const double R = sqrt(dx * dx + dy * dy + dzz * dzz);

            const double fbin = (R - cphd->r_ref[i]) / p->dr + 0.5 * (double)p->n_rbin;
            if (fbin < 0.0 || fbin >= (double)p->n_rbin) continue;

            /* Deposit the compressed response over the bins the envelope
             * covers, carrying the exact propagation phase. */
            const long lo = (long)floor(fbin - 4.0 * sigma / p->dr);
            const long hi = (long)ceil(fbin + 4.0 * sigma / p->dr);
            const double ph = -k_phase * R;
            const double cr = cos(ph), ci = sin(ph);

            for (long b = lo; b <= hi; b++) {
                if (b < 0 || b >= (long)p->n_rbin) continue;
                const double d = ((double)b - fbin) * p->dr;
                const double env = tg->rcs * exp(-0.5 * (d * d) / (sigma * sigma));
                row[b] += (float)(env * cr) + (float)(env * ci) * I;
            }
        }
    }

    return RS_OK;
}

/* Emit the scene description alongside the data, so a later reader knows what
 * was injected without having to re-derive it from the command line. */
static void rs_sim_write_truth(const char *path, const rs_sim_params_t *p,
                               const rs_sim_target_t *t, size_t n)
{
    FILE *f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# ResonarSat synthetic scene ground truth\n");
    fprintf(f, "fc_hz %g\nprf_hz %g\nt_dwell_s %g\nv_platform_ms %g\n",
            p->fc, p->prf, p->t_dwell, p->v_platform);
    fprintf(f, "height_m %g\nrange_offset_m %g\nn_rbin %zu\ndr_m %g\n",
            p->height, p->range_offset, p->n_rbin, p->dr);
    fprintf(f, "# x_m y_m z_m rcs vib_freq_hz vib_amp_m vib_phase_rad\n");
    for (size_t i = 0; i < n; i++) {
        fprintf(f, "target %g %g %g %g %g %g %g\n",
                t[i].x, t[i].y, t[i].z, t[i].rcs,
                t[i].vib_freq, t[i].vib_amp, t[i].vib_phase);
    }
    fclose(f);
}

/* Deterministic uniform generator, so a seed reproduces a scene exactly.
 *
 * A 64-bit LCG rather than rand(): rand() differs between C libraries, and a
 * fixture whose speckle depends on the platform is not a fixture. */
static double rs_sim_uniform(uint64_t *state)
{
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((*state >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
}

/* Fill 'out' with 'n' static scatterers spread over a square of side
 * 'extent_m' centred on (cx, cy), with Rayleigh-distributed reflectivity.
 *
 * WHY THE PIPELINE NEEDS THIS. Correlation-based pixel tracking measures where
 * a PATCH sits, and a patch containing one bright point on empty background is
 * the worst case for it: almost every pixel in the window carries no
 * information, the correlation surface is dominated by the point's own
 * response, and its peak position is biased by the slowly varying shape of that
 * response rather than by the scene. Real structures sit in distributed
 * clutter, which is what makes the published coherence figures achievable and
 * what this generates.
 *
 * It matters most for the master-slave pair. That observable is a small
 * difference across a fixed lag, so it is far more easily buried by correlator
 * bias than a displacement measured against a fixed reference is -- which is
 * exactly the confound a textured scene exists to remove.
 *
 * Rayleigh amplitude is the standard model for the coherent sum of many
 * sub-resolution scatterers, and it is what rs_simulate_static_like() uses, so
 * the two agree on what "distributed" means. */
static void rs_sim_fill_clutter(rs_sim_target_t *out, size_t n,
                                double cx, double cy, double extent_m,
                                double mean_rcs, unsigned seed)
{
    uint64_t state = 0x9e3779b97f4a7c15ULL ^ (uint64_t)seed;
    /* Warm the generator: the first outputs of an LCG seeded with a small
     * integer are strongly correlated across seeds, which would make "distinct
     * realisations" share their first few scatterer positions. */
    for (int i = 0; i < 16; i++) (void)rs_sim_uniform(&state);

    for (size_t i = 0; i < n; i++) {
        out[i] = (rs_sim_target_t){ 0 };
        out[i].x = cx + (rs_sim_uniform(&state) - 0.5) * extent_m;
        out[i].y = cy + (rs_sim_uniform(&state) - 0.5) * extent_m;
        out[i].z = 0.0;

        /* Rayleigh amplitude by inverse transform: sqrt(-2 ln U) scaled so the
         * mean amplitude is 'mean_rcs'. Guarded against U = 0. */
        double u = rs_sim_uniform(&state);
        if (u < 1e-12) u = 1e-12;
        out[i].rcs = mean_rcs * sqrt(-2.0 * log(u)) / 1.2533141;  /* /sqrt(pi/2) */
    }
}

/* Parse one finite numeric argument without atof()'s dangerous "invalid means
 * zero" behaviour. The simulator is used to build long-running experiments, so
 * accepting a misspelled option as a physical zero is much worse than stopping
 * immediately with a precise diagnostic. */
static int rs_sim_number(const char *text, const char *name, double *out)
{
    char *end = NULL;
    errno = 0;
    const double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(value)) {
        fprintf(stderr, "sim_cphd: %s needs a finite number, got '%s'\n",
                name, text);
        return 0;
    }
    *out = value;
    return 1;
}

/* Fetch and parse the value following an option, advancing the argv cursor. */
static int rs_sim_option_number(int argc, char **argv, int *i,
                                const char *name, double *out)
{
    if (*i + 1 >= argc || strncmp(argv[*i + 1], "--", 2) == 0) {
        fprintf(stderr, "sim_cphd: %s needs a value\n", name);
        return 0;
    }
    (*i)++;
    return rs_sim_number(argv[*i], name, out);
}

/* Command-line entry point: build a default scene and write it out. */
int main(int argc, char **argv)
{
    if (argc > 1 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
        printf("usage: sim_cphd OUT.cphd [freq_hz] [amp_m] [offset_x_m] [offset_y_m]\n"
               "                [--clutter N] [--clutter-extent M] [--clutter-rcs F]\n"
               "                [--clutter-vib] [--seed S]\n"
               "                [--offset-x M] [--offset-y M]\n"
               "\n"
               "The offsets translate the whole scene, which is how independent\n"
               "realisations are produced: the geometry and the motion are identical,\n"
               "only the scatterer placement relative to the processing grid changes.\n"
               "Several such scenes are what the sweep needs to separate a real trend\n"
               "from the depth grid's quantisation.\n"
               "\n"
               "--clutter N adds N static Rayleigh scatterers around the vibrating\n"
               "target, giving a DISTRIBUTED-TEXTURE scene instead of an isolated\n"
               "point on empty background. Correlation tracking is biased by its own\n"
               "point response when the window is otherwise empty, and that bias is\n"
               "large next to a differential observable such as --reference pair. Use\n"
               "it whenever the question is about the tracker rather than about\n"
               "focusing. --seed selects the realisation.\n"
               "\n"
               "--clutter-vib makes every clutter scatterer vibrate coherently with\n"
               "the central target, so the patch moves AS A WHOLE -- which is what a\n"
               "structure's surface does, and what patch-level offset tracking is\n"
               "built to measure. Without it the clutter is static and the vibrating\n"
               "point is a minority of every correlation window's energy, so the\n"
               "patch observable is dominated by the static majority however large\n"
               "the injection: measured on the distributed fixture, the tracked\n"
               "series is a function of the clutter seed alone, unchanged by the\n"
               "injected frequency. A lone mover inside static clutter is a\n"
               "detection problem, not the tracking problem this models.\n");
        return 0;
    }

    const char *out_path = (argc > 1) ? argv[1] : "sim.cphd";
    double vib_freq = 2.0;
    double vib_amp = 0.005;
    double off_x = 0.0, off_y = 0.0;
    double clutter_count = 0.0;
    double clutter_extent = 60.0;
    double clutter_rcs = 0.30;
    double seed_value = 1.0;
    int clutter_vib = 0;
    int positional = 0;

    for (int i = 2; i < argc; i++) {
        const char *arg = argv[i];
        if (strcmp(arg, "--clutter-vib") == 0) {
            clutter_vib = 1;
        } else if (strcmp(arg, "--clutter") == 0) {
            if (!rs_sim_option_number(argc, argv, &i, arg, &clutter_count)) return 2;
        } else if (strcmp(arg, "--clutter-extent") == 0) {
            if (!rs_sim_option_number(argc, argv, &i, arg, &clutter_extent)) return 2;
        } else if (strcmp(arg, "--clutter-rcs") == 0) {
            if (!rs_sim_option_number(argc, argv, &i, arg, &clutter_rcs)) return 2;
        } else if (strcmp(arg, "--seed") == 0) {
            if (!rs_sim_option_number(argc, argv, &i, arg, &seed_value)) return 2;
        } else if (strcmp(arg, "--offset-x") == 0) {
            if (!rs_sim_option_number(argc, argv, &i, arg, &off_x)) return 2;
        } else if (strcmp(arg, "--offset-y") == 0) {
            if (!rs_sim_option_number(argc, argv, &i, arg, &off_y)) return 2;
        } else if (strncmp(arg, "--", 2) == 0) {
            fprintf(stderr, "sim_cphd: unknown option '%s'\n", arg);
            return 2;
        } else {
            double *slot = NULL;
            const char *name = NULL;
            switch (positional++) {
            case 0: slot = &vib_freq; name = "freq_hz"; break;
            case 1: slot = &vib_amp;  name = "amp_m"; break;
            case 2: slot = &off_x;    name = "offset_x_m"; break;
            case 3: slot = &off_y;    name = "offset_y_m"; break;
            default:
                fprintf(stderr, "sim_cphd: unexpected positional argument '%s'\n", arg);
                return 2;
            }
            if (!rs_sim_number(arg, name, slot)) return 2;
        }
    }

    if (clutter_count < 0.0 || clutter_count > (double)SIZE_MAX ||
        floor(clutter_count) != clutter_count) {
        fprintf(stderr, "sim_cphd: --clutter needs a non-negative integer\n");
        return 2;
    }
    if (seed_value < 0.0 || seed_value > (double)UINT_MAX ||
        floor(seed_value) != seed_value) {
        fprintf(stderr, "sim_cphd: --seed needs an integer from 0 to %u\n", UINT_MAX);
        return 2;
    }

    rs_sim_params_t p;
    rs_sim_defaults(&p);

    /* A static bright point at the scene centre, a vibrating point offset from
     * it, and two more static points to give the tracker something to lock
     * onto elsewhere in the patch. */
    /* The vibrating target sits at the scene origin so that a processing grid
     * centred on the origin contains it whatever its size. Putting it off-centre
     * -- as an earlier version did -- means a modest grid focuses everything
     * except the one target the run is about, and the results describe the
     * static scatterers instead. */
    const size_t n_clutter = (size_t)clutter_count;
    const unsigned seed = (unsigned)seed_value;

    const size_t n_fixed = 2;
    const size_t n_target = n_fixed + n_clutter;
    rs_sim_target_t *targets = calloc(n_target, sizeof *targets);
    if (!targets) {
        rs_report_error("sim_cphd", RS_ERR_ALLOC);
        return 1;
    }

    targets[0] = (rs_sim_target_t){ .x = 0.0, .y = 0.0, .z = 0.0, .rcs = 1.0,
                                    .vib_freq = vib_freq, .vib_amp = vib_amp,
                                    .vib_phase = 0.0 };
    targets[1] = (rs_sim_target_t){ .x = 9.0, .y = 6.0, .z = 0.0, .rcs = 0.8 };

    /* Clutter is centred on the vibrating target rather than on the scene
     * origin, so that the correlation window around it is filled whatever the
     * scene offset. Centring on the origin would leave the target sitting on
     * empty background again for any non-zero offset -- which is the exact
     * condition this option exists to remove. */
    if (n_clutter > 0) {
        rs_sim_fill_clutter(targets + n_fixed, n_clutter,
                            targets[0].x, targets[0].y,
                            clutter_extent, clutter_rcs, seed);
        if (clutter_vib) {
            for (size_t i = n_fixed; i < n_target; i++) {
                targets[i].vib_freq  = targets[0].vib_freq;
                targets[i].vib_amp   = targets[0].vib_amp;
                targets[i].vib_phase = targets[0].vib_phase;
            }
        }
    }

    for (size_t i = 0; i < n_target; i++) {
        targets[i].x += off_x;
        targets[i].y += off_y;
    }

    rs_cphd_t cphd;
    resonarsat_status_t st = rs_sim_generate(&p, targets, n_target, &cphd);
    if (st != RS_OK) {
        rs_report_error("sim_cphd", st);
        free(targets);
        return 1;
    }

    st = rs_cphd_write(&cphd, out_path);
    if (st != RS_OK) {
        rs_report_error("sim_cphd", st);
        rs_cphd_free(&cphd);
        free(targets);
        return 1;
    }

    char truth_path[512];
    snprintf(truth_path, sizeof truth_path, "%s.truth", out_path);
    rs_sim_write_truth(truth_path, &p, targets, n_target);

    printf("wrote %s: %zu pulses x %zu range bins\n",
           out_path, cphd.n_pulse, cphd.n_rbin);
    printf("  vibrating target at (%g, %g) m: %g Hz, %g m amplitude\n",
           off_x, off_y, vib_freq, vib_amp);
    if (n_clutter > 0) {
        printf("  %zu Rayleigh clutter scatterers over %g m, mean rcs %g, seed %u\n",
               n_clutter, clutter_extent, clutter_rcs, seed);
    }
    printf("  ground truth in %s\n", truth_path);

    rs_cphd_free(&cphd);
    free(targets);
    return 0;
}
