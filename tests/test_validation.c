/* The two experiments that distinguish this implementation from the method it
 * reproduces: the null test and the parameter-sensitivity sweep.
 *
 * Neither is a unit test in the usual sense. They are the evidence a reader
 * needs in order to decide how much any depth this software reports is worth,
 * and they run in CI so that neither can quietly stop being true. */

#include "resonarsat/focus.h"
#include "resonarsat/microm.h"
#include "resonarsat/subaperture.h"
#include "resonarsat/tomo.h"
#include "rs_sim.h"
#include "rs_test.h"

#include <math.h>
#include <stdlib.h>

/* Run the micro-motion chain over a scene and hand back its products. The
 * caller owns everything and must free in reverse order. */
static resonarsat_status_t run_chain(const rs_sim_tgt_t *tg, size_t n_tgt,
                                     size_t n_looks,
                                     rs_cphd_t *c, rs_subap_stack_t *s,
                                     rs_microm_t *m, rs_spectrum_t *spec)
{
    resonarsat_status_t st = rs_sim_scene(c, tg, n_tgt, 20.0, 400.0, 256, 0.5);
    if (st != RS_OK) return st;

    rs_grid_t g = { .origin = {0,0,0}, .n_x = 64, .n_y = 64,
                    .dx = 0.5, .dy = 0.5, .height = 0.0 };

    /* The operating point measured to serve BOTH stages. Micro-motion needs a
     * look count above the phase-ambiguity requirement; the depth stage does
     * better with overlap, which micro-motion is indifferent to once the look
     * count is high enough to compensate for what overlap costs. At 192 looks
     * with 0.4 overlap the effective count is 116, comfortably past the 89 the
     * ambiguity condition demands for this fixture. */
    rs_subap_params_t sp;
    rs_subap_params_default(&sp);
    sp.n_looks = n_looks;
    sp.overlap = 0.4;
    if ((st = rs_subaperture_from_cphd(c, &g, &sp, s)) != RS_OK) return st;

    rs_microm_params_t mp;
    rs_microm_params_default(&mp);
    mp.win_az = mp.win_rg = 32;
    mp.stride_az = mp.stride_rg = 16;
    mp.coherence_min = 0.0;
    if ((st = rs_microm_track(s, &mp, m)) != RS_OK) return st;

    return rs_spectrum_compute(m, RS_SPEC_VELOCITY, spec);
}

/* Mean profile across windows at each depth, and the ratio of its strongest
 * cell to its median. A scene with real structure concentrates energy at one
 * depth; a featureless scene should not. */
static double profile_contrast(const rs_tomo_t *t)
{
    double *mean = calloc(t->n_depth, sizeof *mean);
    if (!mean) return 0.0;

    for (size_t j = 0; j < t->n_depth; j++) {
        for (size_t w = 0; w < t->n_win; w++) mean[j] += t->profile[w * t->n_depth + j];
        mean[j] /= (double)t->n_win;
    }

    double peak = 0.0, sum = 0.0;
    size_t n = 0;
    for (size_t j = 1; j < t->n_depth; j++) {   /* skip depth 0 */
        if (mean[j] > peak) peak = mean[j];
        sum += mean[j];
        n++;
    }
    free(mean);
    return (n && sum > 0.0) ? peak / (sum / (double)n) : 0.0;
}

int main(void)
{
    rs_tomo_params_t p;
    rs_tomo_params_default(&p);
    p.velocity = 3000.0;
    p.frequency = 12500.0;
    /* The depth extent must sit inside what the geometry can represent without
     * folding. With 128 sub-apertures at this geometry that limit is about
     * 29 m; asking for more would previously have zero-filled the remainder and
     * left a single surviving cell that looked like a confident detection. */
    p.depth_max = 10.0;
    p.depth_cell = 1.0;
    p.slant_range = 500000.0;
    p.aperture = 75000.0;

    /* ------------------------------------------------------------------
     * NULL TEST. Featureless terrain -- no scatterers with any motion at
     * all -- must not produce a coherent depth feature. This is the single
     * most important figure this project produces: a method that reports
     * structure over flat nothing is reporting its own artefacts.
     * ------------------------------------------------------------------ */
    RS_CASE("a depth extent beyond the unambiguous range is refused");
    {
        const rs_sim_tgt_t one[] = { { .x = 0.0, .y = 0.0, .z = 0.0, .rcs = 1.0 } };
        rs_cphd_t c; rs_subap_stack_t s; rs_microm_t m; rs_spectrum_t spec;
        RS_CHECK_OK(run_chain(one, 1, 32, &c, &s, &m, &spec));

        rs_tomo_params_t bad = p;
        bad.depth_max = 500.0;      /* far beyond anything 16 looks can reach */
        rs_tomo_t t;
        RS_CHECK_ERR(rs_tomo_focus(&m, &spec, &bad, NULL, &t), RS_ERR_RANGE);

        rs_spectrum_free(&spec); rs_microm_free(&m);
        rs_subap_stack_free(&s); rs_cphd_free(&c);
    }

    RS_CASE("null test: featureless terrain yields no coherent structure");
    double null_contrast = 0.0;
    {
        /* Static scatterers only, spread across the scene: a surface with
         * backscatter but no motion. */
        const rs_sim_tgt_t flat[] = {
            { .x = -12.0, .y =  -8.0, .z = 0.0, .rcs = 1.0 },
            { .x =   5.0, .y =   3.0, .z = 0.0, .rcs = 0.9 },
            { .x =  11.0, .y = -14.0, .z = 0.0, .rcs = 1.1 },
            { .x =  -6.0, .y =  10.0, .z = 0.0, .rcs = 0.8 },
        };

        rs_cphd_t c; rs_subap_stack_t s; rs_microm_t m; rs_spectrum_t spec;
        RS_CHECK_OK(run_chain(flat, sizeof flat / sizeof flat[0], 192, &c, &s, &m, &spec));

        rs_tomo_t t;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t));
        printf("    unambiguous depth range %.1f m over %zu sub-apertures\n",
               t.z_unambiguous, t.n_looks);
        RS_CHECK(t.z_unambiguous >= p.depth_max);
        null_contrast = profile_contrast(&t);
        printf("    null-scene depth contrast: %.2f (1.0 = perfectly flat)\n", null_contrast);

        /* A featureless scene must stay close to flat. The threshold is
         * deliberately generous -- speckle and the finite record guarantee some
         * structure -- but a genuine "discovery" over nothing would be far
         * above it. */
        RS_CHECK(null_contrast < 6.0);

        rs_tomo_free(&t);
        rs_spectrum_free(&spec); rs_microm_free(&m);
        rs_subap_stack_free(&s); rs_cphd_free(&c);
    }

    /* ------------------------------------------------------------------
     * SENSITIVITY SWEEP over scene realisations -- the project's central
     * question, and the one experiment that distinguishes a depth measurement
     * from a relabelled vibration spectrum.
     *
     * Vary the two assumed constants and watch where the recovered feature
     * goes. If depth is nothing but the assumed scaling, the peak depth tracks
     * the acoustic wavelength with slope 1.
     *
     * A single scene is not enough: on a coarse depth grid the peak sometimes
     * lands a cell away from where it lands in another realisation, and such a
     * point carries the grid's quantisation rather than the scene's content.
     * Several realisations are merged, and the fit is taken over the
     * wavelengths that REPRODUCE -- a criterion about repeatability, fixed
     * before looking at the fit, not about which points fall on a line.
     * ------------------------------------------------------------------ */
    RS_CASE("sensitivity sweep: is the depth axis the assumed scaling?");
    {
        rs_tomo_params_t sp2 = p;
        sp2.depth_max = 5.0;
        sp2.depth_cell = 0.7;

        enum { N_SCALE = 5, N_REAL = 3 };
        const double dx[N_REAL] = { 0.0, 7.0, -9.0 };
        const double dy[N_REAL] = { 0.0, -5.0, 6.0 };

        rs_tomo_sweep_row_t all[N_REAL * N_SCALE * N_SCALE];
        size_t n_all = 0;

        for (size_t r = 0; r < N_REAL; r++) {
            const rs_sim_tgt_t mov3[] = {
                { .x = dx[r], .y = dy[r], .z = 0.0, .rcs = 1.0,
                  .vib_freq = 0.5, .vib_amp = 0.020 },
                { .x = 9.0 + dx[r], .y = 6.0 + dy[r], .z = 0.0, .rcs = 0.8 },
            };

            rs_cphd_t c; rs_subap_stack_t s; rs_microm_t m; rs_spectrum_t spec;
            RS_CHECK_OK(run_chain(mov3, 2, 192, &c, &s, &m, &spec));

            rs_tomo_sweep_row_t rows[N_SCALE * N_SCALE];
            size_t n_rows = 0;
            if (rs_tomo_sweep(&m, &spec, &sp2, 0.5, 2.0, N_SCALE, rows, &n_rows) == RS_OK) {
                for (size_t i = 0; i < n_rows && n_all < sizeof all / sizeof all[0]; i++) {
                    all[n_all++] = rows[i];
                }
            }

            rs_spectrum_free(&spec); rs_microm_free(&m);
            rs_subap_stack_free(&s); rs_cphd_free(&c);
        }

        RS_CHECK(n_all >= 3);

        rs_tomo_sweep_row_t merged[N_REAL * N_SCALE * N_SCALE];
        double sd[N_REAL * N_SCALE * N_SCALE];
        size_t n_merged = 0;
        RS_CHECK_OK(rs_tomo_sweep_merge(all, n_all, merged, sd, &n_merged));

        printf("    %zu realisations -> %zu distinct wavelengths\n", (size_t)N_REAL, n_merged);
        printf("    %11s %11s %9s %13s\n", "lambda_ac", "mean depth", "sd", "depth/lambda");
        for (size_t i = 0; i < n_merged; i++) {
            printf("    %11.5f %11.3f %9.3f %13.1f%s\n",
                   merged[i].lambda_ac, merged[i].peak_depth, sd[i],
                   merged[i].peak_depth / merged[i].lambda_ac,
                   (sd[i] > 0.10 * merged[i].peak_depth) ? "   (not reproducible)" : "");
        }

        double slope_all = 0.0, corr_all = 0.0;
        if (rs_tomo_sweep_summary(merged, n_merged, &slope_all, &corr_all) == RS_OK) {
            printf("    all wavelengths:          slope %+.3f, correlation %+.3f\n",
                   slope_all, corr_all);
        }

        double slope = 0.0, corr = 0.0;
        size_t n_used = 0;
        RS_CHECK_OK(rs_tomo_sweep_summary_reliable(merged, sd, n_merged, 0.10,
                                                   &slope, &corr, &n_used));
        printf("    %zu reproducible wavelengths: slope %+.3f, correlation %+.3f\n",
               n_used, slope, corr);

        /* THE RESULT. Depth proportional to the assumed acoustic wavelength
         * means the depth axis is that assumption applied to a fixed spectral
         * feature, and carries no independent depth information. This is what
         * the along-track-baseline objection predicts, and it is asserted here
         * because it is the project's finding rather than an incidental
         * number -- a change in it should fail loudly and be investigated. */
        RS_CHECK(slope > 0.85 && slope < 1.15);
        RS_CHECK(corr > 0.95);

        printf("\n    Depth is proportional to the ASSUMED acoustic wavelength.\n"
               "    The depth axis carries no independent depth information: it is\n"
               "    the operator's (v, f) applied to a fixed spectral feature.\n"
               "    This is the outcome the along-track-baseline geometry predicts.\n\n");
    }

    /* ------------------------------------------------------------------
     * DISCRIMINATION, over several scene realisations.
     *
     * A single pair of scenes cannot answer this: the contrast metric varies
     * substantially with where the scatterers happen to sit, and one ratio above
     * 1 proves nothing. Scatterers are therefore translated between realisations
     * and the two populations compared.
     *
     * WHAT A POSITIVE RESULT MEANS, AND WHAT IT DOES NOT. It means the depth
     * stage responds to motion -- which it must, since its input differs. It
     * does NOT mean the depth axis carries physical depth information. A
     * vibrating target produces a stronger, more concentrated spectral feature,
     * and a more concentrated feature maps to a more peaked depth profile
     * whether or not the mapping means anything. The sensitivity sweep is the
     * test for that, and it remains inconclusive.
     * ------------------------------------------------------------------ */
    RS_CASE("does the depth stage respond to motion at all?");
    {
        const double dx[] = { 0.0, 7.0, -9.0 };
        const double dy[] = { 0.0, -5.0, 6.0 };
        const size_t n_real = sizeof dx / sizeof dx[0];

        double sum_null = 0.0, sum_vib = 0.0;
        size_t n_above = 0;

        printf("    %6s %10s %10s %8s\n", "scene", "null", "vibrating", "ratio");
        for (size_t r = 0; r < n_real; r++) {
            const rs_sim_tgt_t flat2[] = {
                { .x = -12.0 + dx[r], .y =  -8.0 + dy[r], .z = 0.0, .rcs = 1.0 },
                { .x =   5.0 + dx[r], .y =   3.0 + dy[r], .z = 0.0, .rcs = 0.9 },
                { .x =  11.0 + dx[r], .y = -14.0 + dy[r], .z = 0.0, .rcs = 1.1 },
            };
            const rs_sim_tgt_t mov2[] = {
                { .x = dx[r], .y = dy[r], .z = 0.0, .rcs = 1.0,
                  .vib_freq = 0.5, .vib_amp = 0.020 },
                { .x = 9.0 + dx[r], .y = 6.0 + dy[r], .z = 0.0, .rcs = 0.8 },
            };

            double cn = 0.0, cv = 0.0;
            for (int k = 0; k < 2; k++) {
                rs_cphd_t c; rs_subap_stack_t s; rs_microm_t m; rs_spectrum_t spec;
                const rs_sim_tgt_t *tg = k ? mov2 : flat2;
                const size_t nt = k ? 2 : 3;
                RS_CHECK_OK(run_chain(tg, nt, 192, &c, &s, &m, &spec));

                rs_tomo_t t;
                RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t));
                const double ct = profile_contrast(&t);
                if (k) cv = ct; else cn = ct;
                rs_tomo_free(&t);

                rs_spectrum_free(&spec); rs_microm_free(&m);
                rs_subap_stack_free(&s); rs_cphd_free(&c);
            }

            printf("    %6zu %10.2f %10.2f %8.2f\n", r, cn, cv, (cn > 0.0) ? cv / cn : 0.0);
            sum_null += cn;
            sum_vib  += cv;
            if (cv > cn) n_above++;
        }

        const double mean_null = sum_null / (double)n_real;
        const double mean_vib  = sum_vib / (double)n_real;
        printf("    mean null %.2f, mean vibrating %.2f, %zu of %zu realisations higher\n",
               mean_null, mean_vib, n_above, n_real);
        printf("    (responding to motion is not the same as measuring depth --\n"
               "     see the sensitivity sweep)\n");

        /* Assert the population difference, not a single ratio. */
        RS_CHECK(mean_vib > mean_null);
        RS_CHECK(n_above >= n_real - 1);
    }

    RS_TEST_END();
}
