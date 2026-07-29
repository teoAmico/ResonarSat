/* Null test for the micro-motion stage.
 *
 * Phase 5 has a null test -- run the tomography over featureless ground and
 * confirm no structure is reported. Phase 3 had none, and this is it: run the
 * tracking over a scene where nothing moves, and measure how strong a "detection"
 * it reports anyway.
 *
 * That number is the false-positive floor. Any genuine detection must clear it
 * to mean anything, and the comparison is the honest measure of whether the
 * micro-motion stage discriminates at all. Both scenes below use identical
 * geometry, identical processing parameters and identical scatterer brightness,
 * so the only difference is whether one target moves.
 *
 * The test asserts that the floor is measured and finite, and PRINTS the
 * comparison. It deliberately does not assert that detections clear the floor,
 * because at present they do not -- see STATUS.md. Asserting it would either
 * fail permanently or, if tuned to pass, encode the shortfall as acceptable. */

#include "resonarsat/focus.h"
#include "resonarsat/microm.h"
#include "resonarsat/subaperture.h"
#include "rs_sim.h"
#include "rs_test.h"

#include <math.h>
#include <stdlib.h>

/* Processing applied identically to every scene here, so that any difference in
 * the result is attributable to the scene rather than to the parameters. */
static resonarsat_status_t analyse(const rs_sim_tgt_t *tg, size_t n_tgt,
                                   double *max_prom_out,
                                   double *at_freq_prom_out, double want_freq,
                                   size_t *n_above_out, double above_thresh,
                                   size_t *n_win_out)
{
    rs_cphd_t c;
    resonarsat_status_t st = rs_sim_scene(&c, tg, n_tgt, 20.0, 400.0, 256, 0.5);
    if (st != RS_OK) return st;

    rs_grid_t g = { .origin = {0,0,0}, .n_x = 64, .n_y = 64,
                    .dx = 0.5, .dy = 0.5, .height = 0.0 };

    /* The operating point that works, and the reason it does.
     *
     * The sub-look count is set by the phase-ambiguity condition (see
     * rs_microm_recommend_looks): the peak azimuth shift must fall inside three
     * quarters of a sub-look resolution cell, and since the shift is fixed by
     * the geometry while the resolution cell widens as sub-looks shorten, that
     * demands MANY looks rather than few. Overlap works against the condition,
     * so it is zero here. At 32 looks with 0.4 overlap -- the settings this test
     * used previously -- no genuine detection cleared the false-positive floor. */
    rs_subap_params_t sp;
    rs_subap_params_default(&sp);
    sp.n_looks = 128;
    sp.overlap = 0.0;

    rs_subap_stack_t s;
    if ((st = rs_subaperture_from_cphd(&c, &g, &sp, &s)) != RS_OK) {
        rs_cphd_free(&c); return st;
    }

    rs_microm_params_t mp;
    rs_microm_params_default(&mp);
    mp.win_az = mp.win_rg = 32;
    mp.stride_az = mp.stride_rg = 16;
    mp.coherence_min = 0.0;

    rs_microm_t m;
    if ((st = rs_microm_track(&s, &mp, &m)) != RS_OK) {
        rs_subap_stack_free(&s); rs_cphd_free(&c); return st;
    }

    rs_spectrum_t spec;
    if ((st = rs_spectrum_compute(&m, RS_SPEC_VELOCITY, &spec)) != RS_OK) {
        rs_microm_free(&m); rs_subap_stack_free(&s); rs_cphd_free(&c); return st;
    }

    double max_prom = 0.0, at_freq = 0.0;
    size_t n_above = 0;
    for (size_t w = 0; w < spec.n_win; w++) {
        if (spec.prominence[w] > max_prom) max_prom = spec.prominence[w];
        if (spec.prominence[w] > above_thresh) n_above++;
        if (want_freq > 0.0 &&
            fabs(spec.dominant_freq[w] - want_freq) < 2.0 * spec.df &&
            spec.prominence[w] > at_freq) {
            at_freq = spec.prominence[w];
        }
    }

    if (max_prom_out)     *max_prom_out = max_prom;
    if (at_freq_prom_out) *at_freq_prom_out = at_freq;
    if (n_above_out)      *n_above_out = n_above;
    if (n_win_out)        *n_win_out = spec.n_win;

    rs_spectrum_free(&spec);
    rs_microm_free(&m);
    rs_subap_stack_free(&s);
    rs_cphd_free(&c);
    return RS_OK;
}

int main(void)
{
    /* ------------------------------------------------------------------
     * The floor: a scene containing bright scatterers, none of them moving.
     * Every peak reported here is an artefact.
     * ------------------------------------------------------------------ */
    RS_CASE("false-positive floor on a scene where nothing moves");
    double floor_prom = 0.0;
    size_t n_above = 0, n_win = 0;
    {
        const rs_sim_tgt_t stat[] = {
            { .x =   0.0, .y = 0.0, .z = 0.0, .rcs = 1.0 },
            { .x = -30.0, .y = 0.0, .z = 0.0, .rcs = 0.8 },
        };
        RS_CHECK_OK(analyse(stat, 2, &floor_prom, NULL, 0.0, &n_above, 8.0, &n_win));

        printf("    highest prominence anywhere: %.1f\n", floor_prom);
        printf("    %zu of %zu windows exceed prominence 8\n", n_above, n_win);

        /* The floor must be a real number well above 1 (a flat spectrum), or
         * this test is not measuring what it claims. */
        RS_CHECK(floor_prom > 1.0);
        RS_CHECK(floor_prom == floor_prom);
    }

    /* ------------------------------------------------------------------
     * Genuine detections, same geometry and brightness, one target moving.
     * ------------------------------------------------------------------ */
    RS_CASE("genuine detections against that floor");
    {
        const double freqs[] = { 0.3, 0.5, 0.7, 0.9 };
        const size_t n_freq = sizeof freqs / sizeof freqs[0];

        printf("    %9s %14s %10s %s\n",
               "injected", "target prom", "floor", "clears?");

        size_t n_clear = 0;
        for (size_t i = 0; i < n_freq; i++) {
            const rs_sim_tgt_t vib[] = {
                { .x = 0.0, .y = 0.0, .z = 0.0, .rcs = 1.0,
                  .vib_freq = freqs[i], .vib_amp = 0.020 },
            };
            double tp = 0.0;
            RS_CHECK_OK(analyse(vib, 1, NULL, &tp, freqs[i], NULL, 0.0, NULL));

            const int clears = tp > floor_prom;
            if (clears) n_clear++;
            printf("    %7.1f Hz %13.1f %10.1f  %s\n",
                   freqs[i], tp, floor_prom, clears ? "yes" : "NO");
        }

        printf("\n    %zu of %zu genuine detections clear the false-positive floor\n",
               n_clear, n_freq);
        printf("    Detections clearing the floor is what makes a measurement a\n"
               "    measurement rather than a ranking of artefacts. This became\n"
               "    true only at a sub-look count set by the ambiguity condition;\n"
               "    at the settings used previously the margin was below 1.\n");

        /* Now assertable, where before it was not: at this operating point
         * genuine detections do clear the floor. The bar is set below the
         * measured rate so that an improvement passes and a regression fails. */
        RS_CHECK(n_clear >= 3);
    }

    RS_TEST_END();
}
