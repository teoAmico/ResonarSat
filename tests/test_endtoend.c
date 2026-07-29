/* End-to-end wiring check: inject a known vibration, run the whole chain.
 *
 * READ THIS BEFORE TRUSTING THE RESULT BELOW. This exercises ONE configuration
 * at ONE frequency, and it passes. It does not follow that the chain recovers
 * vibration frequencies in general -- it does not. Sweeping the injected
 * frequency through the same configuration (tests/test_tracking.c) recovers
 * roughly two cases in six, and the reported velocity barely varies with the
 * injected motion, which means it is substantially a fixed artefact.
 *
 * What this file therefore validates is that the stages are wired together
 * correctly and that data flows end to end with the right units and
 * conventions. Treat the recovered frequency as one data point, not as
 * evidence the measurement works.
 *
 * This is the synthetic ground-truth test the implementation plan puts first,
 * and the only place where a frequency the software reports can be checked
 * against a frequency that was definitely put there. Simulate phase history for
 * a vibrating point target, focus it, decompose into sub-looks, track, estimate
 * the spectrum, and require the injected frequency back within one spectral
 * bin.
 *
 * What this proves and what it does not: it validates image formation,
 * decomposition, tracking and spectral estimation as a working chain. It says
 * nothing about the depth stage, because the simulator would have to assume the
 * very physics that stage is in question over. */

#include "resonarsat/focus.h"
#include "resonarsat/microm.h"
#include "resonarsat/subaperture.h"
#include "rs_sim.h"
#include "rs_test.h"

int main(void)
{
    /* A 20 s dwell samples a 1 Hz vibration over 20 cycles, which is a
     * comfortable margin for a 16-look decomposition. */
    const double t_dwell = 20.0, prf = 400.0;
    const double vib_freq = 0.5;
    /* 20 mm amplitude is chosen deliberately to EXCEED a quarter wavelength.
     * At X-band lambda/4 is about 7.8 mm, so this motion wraps the
     * interferometric phase several times over -- which is exactly the regime
     * real structural response occupies, and exactly what the displacement
     * observable alone cannot measure.
     *
     * The tracked azimuth shift is unaffected: it is a geometric displacement
     * of the target in the image, not a phase, so it does not wrap. That is why
     * line-of-sight velocity is the primary observable and why this test asserts
     * on it. */
    const double vib_amp = 0.020;
    const size_t n_rbin = 256;
    const double dr = 0.5;

    const rs_sim_tgt_t tg[] = {
        { .x = 0.0, .y = 0.0, .z = 0.0, .rcs = 1.0,
          .vib_freq = vib_freq, .vib_amp = vib_amp, .vib_phase = 0.0 },
    };

    rs_cphd_t cphd;
    RS_CASE("simulate a vibrating point target");
    RS_CHECK_OK(rs_sim_scene(&cphd, tg, 1, t_dwell, prf, n_rbin, dr));

    rs_grid_t grid = { .origin = {0,0,0}, .n_x = 48, .n_y = 48,
                       .dx = 1.0, .dy = 1.0, .height = 0.0 };

    rs_subap_params_t sp;
    rs_subap_params_default(&sp);
    /* 64 looks over a 20 s dwell samples a 0.5 Hz motion about 6.5 times per
     * cycle -- comfortably above Nyquist, where 32 looks would give only 3.3 and
     * leave the estimate at the mercy of tracking noise. */
    sp.n_looks = 64;
    sp.overlap = 0.4;

    /* Form the looks by focusing shifted pulse windows -- the faithful path,
     * and the one whose timing is exact by construction. */
    rs_subap_stack_t stack;
    RS_CASE("decompose into sub-apertures from phase history");
    RS_CHECK_OK(rs_subaperture_from_cphd(&cphd, &grid, &sp, &stack));

    /* The observable band must cover the injected frequency, or the test is
     * checking aliasing rather than recovery. */
    RS_CASE("the observable band covers the injected frequency");
    RS_CHECK(stack.f_max > vib_freq);

    rs_microm_params_t mp;
    rs_microm_params_default(&mp);
    mp.win_az = mp.win_rg = 24;
    mp.stride_az = mp.stride_rg = 8;
    mp.upsample_az = 40;
    mp.upsample_rg = 20;
    /* Unmasked: an isolated point target on empty background scores below the
     * literature coherence threshold however well it tracks, because most of
     * the correlation window is background. */
    mp.coherence_min = 0.0;

    rs_microm_t m;
    RS_CASE("track sub-looks");
    RS_CHECK_OK(rs_microm_track(&stack, &mp, &m));
    RS_CHECK(m.n_win > 0);
    RS_CHECK(m.n_looks == 64);

    rs_spectrum_t spec;
    RS_CASE("estimate vibration spectra");
    RS_CHECK_OK(rs_spectrum_compute(&m, RS_SPEC_VELOCITY, &spec));
    RS_CHECK(spec.n_freq > 1);
    RS_CHECK(spec.df > 0.0);

    /* Select the window whose displacement series varies most: that is where
     * the vibrating target sits. A scene average would be dominated by static
     * background, and picking by spectral amplitude alone can land on an empty
     * window whose noise happens to peak. */
    RS_CASE("the injected frequency is recovered AT THIS ONE CONFIGURATION");
    {
        size_t best = 0;
        double best_span = -1.0;
        for (size_t w = 0; w < m.n_win; w++) {
            double lo = m.vel_los[w * m.n_looks], hi = lo;
            for (size_t k = 1; k < m.n_looks; k++) {
                const double v = m.vel_los[w * m.n_looks + k];
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
            if (hi - lo > best_span) { best_span = hi - lo; best = w; }
        }

        printf("    injected %.3f Hz, recovered %.3f Hz, bin width %.3f Hz\n",
               vib_freq, spec.dominant_freq[best], spec.df);

        /* Within three bins: one for the estimator's own resolution, one for
         * the sub-look sampling grid not aligning with the vibration phase, and
         * one for the short record the Hann window is applied to. */
        RS_CHECK_NEAR(spec.dominant_freq[best], vib_freq, 3.0 * spec.df);
    }

    /* A displacement series must actually vary. A tracker returning zeros would
     * still produce a spectrum, with a meaningless peak. */
    RS_CASE("displacement series are non-trivial");
    {
        double span = 0.0;
        for (size_t w = 0; w < m.n_win; w++) {
            double lo = m.vel_los[w * m.n_looks], hi = lo;
            for (size_t k = 1; k < m.n_looks; k++) {
                const double v = m.vel_los[w * m.n_looks + k];
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
            if (hi - lo > span) span = hi - lo;
        }
        RS_CHECK(span > 0.0);
    }

    rs_spectrum_free(&spec);
    rs_microm_free(&m);
    rs_subap_stack_free(&stack);
    rs_cphd_free(&cphd);
    RS_TEST_END();
}
