/* Sub-aperture decomposition: band layout, timing, and the resolution trade. */

#include "resonarsat/subaperture.h"
#include "resonarsat/geom.h"
#include "rs_test.h"

#include <math.h>
#include <stdlib.h>

/* Build a focused image containing a single bright point plus a weak
 * background, with self-consistent azimuth timing. */
static resonarsat_status_t make_image(rs_slc_t *img, size_t n, double t_dwell)
{
    resonarsat_status_t st = rs_slc_alloc(img, n, n);
    if (st != RS_OK) return st;

    img->fc = 9.6e9;
    img->azimuth_time_interval = t_dwell / (double)n;
    img->t_dwell = t_dwell;
    img->v_platform = 7500.0;
    img->r0 = 500000.0;
    img->az_spacing_m = img->rg_spacing_m = 1.0;
    st = rs_slc_finalise_metadata(img);
    if (st != RS_OK) return st;

    for (size_t a = 0; a < n; a++) {
        for (size_t r = 0; r < n; r++) {
            const double da = ((double)a - (double)n / 2) / 2.0;
            const double dr = ((double)r - (double)n / 2) / 2.0;
            img->data[a * n + r] = (float)(exp(-0.5 * (da * da + dr * dr))
                                  + 0.001 * ((double)((a * 7 + r * 13) % 17) / 17.0));
        }
    }
    return RS_OK;
}

int main(void)
{
    const size_t n = 128;
    const double t_dwell = 10.0;

    rs_slc_t img;
    RS_CASE("build a test image");
    RS_CHECK_OK(make_image(&img, n, t_dwell));

    RS_CASE("uniform mode produces the requested looks");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.n_looks = 16;
        p.overlap = 0.4;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        RS_CHECK(s.n_looks == 16);
        RS_CHECK(s.dt > 0.0);
        RS_CHECK(s.f_max > 0.0);

        /* f_max and dt are two views of one quantity and must agree. */
        RS_CHECK_NEAR(s.f_max, 1.0 / (2.0 * s.dt), 1e-9);

        /* Every look must carry energy: a band filter that zeroed a look
         * would leave a silently dead entry in the stack. */
        for (size_t k = 0; k < s.n_looks; k++) {
            double e = 0.0;
            for (size_t i = 0; i < n * n; i++) e += (double)cabsf(s.look[k].data[i]);
            RS_CHECK(e > 0.0);
        }

        /* Sub-look resolution must be coarser than the full aperture's. */
        const double full = rs_azimuth_resolution(img.lambda, img.r0,
                                                  img.v_platform, t_dwell);
        RS_CHECK(s.az_resolution > full);

        rs_subap_stack_free(&s);
    }

    /* More looks buy a wider observable band and cost resolution. This is the
     * project's central trade and it must show up in the reported numbers, not
     * only in the documentation. */
    RS_CASE("more looks widen the band and coarsen the resolution");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.overlap = 0.0;

        rs_subap_stack_t s8, s32;
        p.n_looks = 8;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s8));
        p.n_looks = 32;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s32));

        RS_CHECK(s32.f_max > s8.f_max);
        RS_CHECK(s32.az_resolution > s8.az_resolution);
        RS_CHECK(s32.t_sap < s8.t_sap);

        rs_subap_stack_free(&s8);
        rs_subap_stack_free(&s32);
    }

    RS_CASE("paper mode runs and holds out bandwidth");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 16;
        p.left_out_frac = 0.5;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        RS_CHECK(s.n_looks == 16);
        /* Each look processes half the total band, so its aperture time is
         * half the dwell -- much longer than a uniform look, which is why the
         * paper's looks overlap so heavily. */
        RS_CHECK(s.t_sap > 0.4 * t_dwell);
        rs_subap_stack_free(&s);
    }

    RS_CASE("paper mode reports derived B_shift");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 16;
        p.left_out_frac = 0.5;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        const double expected = (s.doppler_bandwidth -
                                 s.doppler_bandwidth * p.left_out_frac) /
                                (double)p.n_looks;
        RS_CHECK_NEAR(s.b_shift_hz, expected, 1e-9 * fmax(1.0, expected));
        rs_subap_stack_free(&s);
    }

    RS_CASE("B_shift is free when set, derived when not");
    {
        /* WO2024008365A1 [0004] makes B_shift an operator choice independent of
         * the sweep step. Setting it must be honoured, not quietly overridden by
         * the step, or the parameter the patent calls "the precise vibrational
         * frequency one wishes to observe" would have no effect. */
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 16;
        p.left_out_frac = 0.5;
        p.b_shift_hz = 3.0;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        RS_CHECK_NEAR(s.b_shift_hz, 3.0, 1e-12);
        /* The gap is a time lag: dt = B_shift * t_dwell / B_CD. */
        const double expect_lag = 3.0 * t_dwell / s.doppler_bandwidth;
        RS_CHECK_NEAR(s.pair_lag_s, expect_lag, 1e-9 * fmax(1.0, expect_lag));
        rs_subap_stack_free(&s);
    }

    RS_CASE("pair mode builds slave bands offset by B_shift");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 8;
        p.left_out_frac = 0.5;
        p.pair = 1;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        RS_CHECK(s.slave != NULL);
        RS_CHECK(s.b_shift_hz > 0.0);
        for (size_t i = 0; i < s.n_looks; i++) {
            RS_CHECK(s.slave[i].data != NULL);
            RS_CHECK(s.slave[i].data != s.look[i].data);
            RS_CHECK(s.slave[i].n_az == s.look[i].n_az);
            RS_CHECK(s.slave[i].n_rg == s.look[i].n_rg);
        }
        rs_subap_stack_free(&s);
    }

    RS_CASE("without pair the slave array stays absent");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 8;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        RS_CHECK(s.slave == NULL);
        rs_subap_stack_free(&s);
    }

    RS_CASE("pair mode is refused outside the paper band layout");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_UNIFORM;
        p.n_looks = 8;
        p.pair = 1;

        rs_subap_stack_t s;
        RS_CHECK_ERR(rs_subaperture_split(&img, &p, &s), RS_ERR_ARG);
    }

    /* The pair is "rigidly held at a distance B_shift" for the WHOLE sweep
     * (WO2024008365A1 [0004]). If the last slave runs past the edge of the
     * measured Doppler support, the band that filter passes is clipped, the
     * separation stops being rigid at the end of the sweep, and the resulting
     * master-slave offset drifts with sweep position -- a systematic
     * low-frequency term added to exactly the series the method reads a
     * vibration out of. That must be refused rather than produced quietly. */
    RS_CASE("the rigid pair is kept inside the Doppler support");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 16;
        p.left_out_frac = 0.5;
        p.pair = 1;

        /* The default B_shift is the sweep step, and the layout is sized so
         * that the last slave lands exactly on the band edge: it fits. */
        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        const double step = s.b_shift_hz;
        RS_CHECK(step > 0.0);
        rs_subap_stack_free(&s);

        /* One step of headroom is all the layout has, so twice the step cannot
         * fit and the split must say so instead of clipping the last slave. */
        p.b_shift_hz = 2.0 * step;
        RS_CHECK_ERR(rs_subaperture_split(&img, &p, &s), RS_ERR_ARG);

        /* The same oversized B_shift without a slave to place is fine: with no
         * pair there is nothing that can fall off the edge. */
        p.pair = 0;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        RS_CHECK(s.slave == NULL);
        rs_subap_stack_free(&s);
    }

    /* THE ONE PUBLISHED B_SHIFT, AND WHERE IT FITS.
     *
     * Biondi's power-line paper (Preprints 2023, doi
     * 10.20944/preprints202308.0926.v1) states B_shift = B_CD/100 -- the only
     * value any source gives. The sweep geometry allows at most one step of
     * headroom, and with B_DL = B_CD/2 that step is B_CD/(2*N_D), so the
     * published value is admissible only for N_D <= 50 and coincides with this
     * code's derived default exactly at N_D = 50.
     *
     * That arithmetic is argued in rs_subap_params_t.b_shift_hz and is checked
     * here, because a claim about which published configurations are
     * representable should not rest on a comment. */
    RS_CASE("the published B_shift fits only up to fifty sub-apertures");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.left_out_frac = 0.5;
        p.pair = 1;

        /* The Doppler band the split will measure, taken from the stack itself
         * so the test uses the same number the layout does. */
        p.n_looks = 50;
        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        const double bw = s.doppler_bandwidth;
        const double step_50 = s.b_shift_hz;
        rs_subap_stack_free(&s);

        const double published = bw / 100.0;
        printf("    B_CD %.4g Hz; published B_CD/100 = %.4g Hz; "
               "sweep step at 50 looks = %.4g Hz\n", bw, published, step_50);

        /* Equal at fifty, to floating-point noise. */
        RS_CHECK_REL(published, step_50, 1e-9);
        p.b_shift_hz = published;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        rs_subap_stack_free(&s);

        /* Below fifty the step is wider than the published value, so it fits. */
        p.n_looks = 32;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        rs_subap_stack_free(&s);

        /* Above fifty the step has shrunk past it and the layout refuses it.
         * At the 128 looks the Giza runs used, the published configuration is
         * not representable at all. */
        p.n_looks = 128;
        printf("    at 128 looks the published value is %.3f of the sweep step\n",
               published / (bw / (2.0 * 128.0)));
        RS_CHECK_ERR(rs_subaperture_split(&img, &p, &s), RS_ERR_ARG);
    }

    /* THE DEFAULT B_shift COLLAPSES THE PAIR, and this measures it rather than
     * arguing it. B_shift defaults to the sweep step, the slave of look k is
     * therefore centred exactly where the master of look k+1 is, and both are
     * cut from the same spectrum with the same width and the same taper. So
     * slave[k] is not merely similar to master[k+1] -- it is the same image,
     * sample for sample.
     *
     * What that costs: RS_MICROM_REF_PAIR then measures the offset between
     * consecutive masters, which is ADJACENT's differential without ADJACENT's
     * accumulation, not an independent slave observation. The patent's
     * construction is only distinct from adjacent-look differencing when
     * B_shift differs from the step, and the patent gives no value for it.
     *
     * See also the case above: the layout holds exactly one step of headroom,
     * so every B_shift that escapes this collapse is SMALLER than the step, and
     * a smaller lag attenuates the differential by |2 sin(pi f dt)|. There is no
     * way out of the degenerate case in the direction [0004] describes as the
     * useful one ("the higher this parameter, the lower the mechanical
     * frequency observed"). */
    RS_CASE("the default B_shift makes each slave its neighbour's master");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        p.mode = RS_SUBAP_PAPER;
        p.n_looks = 8;
        p.left_out_frac = 0.5;
        p.pair = 1;

        rs_subap_stack_t s;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));

        double worst = 0.0, scale = 0.0;
        for (size_t k = 0; k + 1 < s.n_looks; k++) {
            const rs_slc_t *sl = &s.slave[k], *ma = &s.look[k + 1];
            for (size_t i = 0; i < sl->n_az * sl->n_rg; i++) {
                const double d = cabs((double complex)sl->data[i]
                                    - (double complex)ma->data[i]);
                const double m = cabs((double complex)ma->data[i]);
                if (d > worst)  worst = d;
                if (m > scale)  scale = m;
            }
        }
        printf("    max|slave[k] - master[k+1]| = %.3e against a peak of %.3e\n",
               worst, scale);
        RS_CHECK(scale > 0.0);
        RS_CHECK(worst == 0.0);
        const double step = s.b_shift_hz;
        rs_subap_stack_free(&s);

        /* Half a step puts the slave between two masters, so it becomes a band
         * neither master covers and the collapse must disappear. Asserting this
         * too keeps the case above from passing for a trivial reason -- an
         * all-zero stack would satisfy an equality test on its own. */
        p.b_shift_hz = 0.5 * step;
        RS_CHECK_OK(rs_subaperture_split(&img, &p, &s));
        worst = 0.0;
        for (size_t k = 0; k + 1 < s.n_looks; k++) {
            const rs_slc_t *sl = &s.slave[k], *ma = &s.look[k + 1];
            for (size_t i = 0; i < sl->n_az * sl->n_rg; i++) {
                const double d = cabs((double complex)sl->data[i]
                                    - (double complex)ma->data[i]);
                if (d > worst) worst = d;
            }
        }
        printf("    at half a step: %.3e (%.1f%% of peak)\n",
               worst, 100.0 * worst / scale);
        RS_CHECK(worst > 0.01 * scale);
        rs_subap_stack_free(&s);
    }

    RS_CASE("degenerate parameters are refused");
    {
        rs_subap_params_t p;
        rs_subap_params_default(&p);
        rs_subap_stack_t s;

        p.n_looks = 0;
        RS_CHECK_ERR(rs_subaperture_split(&img, &p, &s), RS_ERR_ARG);

        p.n_looks = n;   /* needs 2*n azimuth lines, image has only n */
        RS_CHECK_ERR(rs_subaperture_split(&img, &p, &s), RS_ERR_ARG);

        p.n_looks = 8;
        p.mode = RS_SUBAP_PAPER;
        p.left_out_frac = 1.5;
        RS_CHECK_ERR(rs_subaperture_split(&img, &p, &s), RS_ERR_ARG);
    }

    RS_CASE("doppler estimators run on a real image");
    {
        double fdc = 0.0, bw = 0.0;
        RS_CHECK_OK(rs_estimate_doppler_centroid(&img, &fdc));
        RS_CHECK_OK(rs_estimate_doppler_bandwidth(&img, 0.5, &bw));
        RS_CHECK(bw > 0.0);
        RS_CHECK(bw <= img.fs_az);
        RS_CHECK_ERR(rs_estimate_doppler_bandwidth(&img, 1.5, &bw), RS_ERR_ARG);
    }

    rs_slc_free(&img);
    RS_TEST_END();
}
