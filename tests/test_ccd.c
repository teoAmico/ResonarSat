/* Tests for the scale-invariant CCD micro-motion locator.
 *
 * The first case is the one that matters, and it needs no scene: the statistic
 * is defined to be invariant to a positive scaling between the two covariances,
 * so multiplying a whole sub-look by any constant must leave the map unchanged.
 * That invariance is the entire reason to prefer this detector over an
 * amplitude-based one -- it is what keeps a bright but stationary target quiet
 * -- and it is checkable in closed form, so it is asserted tightly.
 *
 * The second case builds a scene with one vibrating and one static target of
 * equal brightness and REPORTS the contrast between them rather than asserting
 * a threshold. That follows tests/test_nullmotion.c: the source paper
 * implements no detection threshold and offers no false-alarm rate, so there is
 * no published figure to assert against, and inventing one here would encode a
 * number nobody has justified. The contrast is printed so a regression is
 * visible to a reader even though it is not gated. */

#include "resonarsat/ccd.h"
#include "resonarsat/focus.h"
#include "resonarsat/subaperture.h"
#include "rs_sim.h"
#include "rs_test.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Mean of the map over a square of side 'side' centred on (row, col).
 *
 * Used to compare two targets without depending on exactly which pixel a peak
 * lands in, which shifts with the grid origin. */
static double rs_ccd_patch_mean(const rs_ccd_t *ccd, size_t row, size_t col, size_t side)
{
    const size_t half = side / 2;
    double sum = 0.0;
    size_t n = 0;
    for (size_t r = row - half; r <= row + half; r++) {
        for (size_t c = col - half; c <= col + half; c++) {
            if (r < ccd->n_row && c < ccd->n_col) {
                sum += ccd->map[r * ccd->n_col + c];
                n++;
            }
        }
    }
    return n ? sum / (double)n : 0.0;
}

int main(void)
{
    /* A stack of pure noise is enough for the invariance case: the identity
     * holds per window regardless of what the windows contain. */
    RS_CASE("the map is invariant to a positive scaling of the whole stack");
    {
        const size_t n_looks = 5, n = 24;
        rs_subap_stack_t stack;
        memset(&stack, 0, sizeof stack);
        stack.n_looks = n_looks;
        stack.look = calloc(n_looks, sizeof *stack.look);
        RS_CHECK(stack.look != NULL);

        unsigned seed = 12345u;
        for (size_t i = 0; i < n_looks; i++) {
            RS_CHECK_OK(rs_slc_alloc(&stack.look[i], n, n));
            for (size_t p = 0; p < n * n; p++) {
                /* xorshift, so the fixture is identical on every platform. */
                seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5;
                const double re = (double)(seed % 2000) / 1000.0 - 1.0;
                seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5;
                const double im = (double)(seed % 2000) / 1000.0 - 1.0;
                stack.look[i].data[p] = (float)re + (float)im * I;
            }
        }

        rs_ccd_params_t cp;
        rs_ccd_params_default(&cp);

        rs_ccd_t before;
        RS_CHECK_OK(rs_ccd_locate(&stack, &cp, &before));

        /* Scale the WHOLE stack. Both covariances scale together, so
         * S_X S_Y^-1 is untouched and the map must not move at all -- including
         * the diagonal loading, which is derived from the stack's own mean power
         * and therefore scales with it. Scaling a single look would NOT test
         * this; see rs_ccd_stat_t. */
        const double gamma = 37.5;
        for (size_t i = 0; i < n_looks; i++) {
            for (size_t p = 0; p < n * n; p++) {
                stack.look[i].data[p] *= (float)gamma;
            }
        }

        rs_ccd_t after;
        RS_CHECK_OK(rs_ccd_locate(&stack, &cp, &after));

        double worst = 0.0;
        for (size_t p = 0; p < before.n_row * before.n_col; p++) {
            const double d = fabs(after.map[p] - before.map[p]);
            const double rel = before.map[p] > 0.0 ? d / before.map[p] : d;
            if (rel > worst) worst = rel;
        }
        printf("    worst relative change under a %.1fx scaling: %.3g\n", gamma, worst);
        RS_CHECK(worst < 1e-6);

        rs_ccd_free(&before);
        rs_ccd_free(&after);
        rs_subap_stack_free(&stack);
    }

    /* A degenerate stack must be refused rather than producing an empty map,
     * which would read downstream as a scene with nothing moving in it. */
    RS_CASE("a stack of fewer than three looks is refused");
    {
        rs_subap_stack_t stack;
        memset(&stack, 0, sizeof stack);
        stack.n_looks = 2;
        stack.look = calloc(2, sizeof *stack.look);
        RS_CHECK(stack.look != NULL);
        for (size_t i = 0; i < 2; i++) RS_CHECK_OK(rs_slc_alloc(&stack.look[i], 16, 16));

        rs_ccd_params_t cp;
        rs_ccd_params_default(&cp);
        rs_ccd_t ccd;
        RS_CHECK_ERR(rs_ccd_locate(&stack, &cp, &ccd), RS_ERR_ARG);
        rs_subap_stack_free(&stack);
    }

    /* One vibrating target and one static target of identical brightness, on a
     * grid wide enough to hold both with clear separation. */
    RS_CASE("a vibrating target against a static one of equal brightness");
    {
        const rs_sim_tgt_t tg[2] = {
            { -20.0, 0.0, 0.0, 1.0, 2.0, 0.005, 0.0 },   /* vibrating */
            {  20.0, 0.0, 0.0, 1.0, 0.0, 0.0,   0.0 },   /* static, same RCS */
        };

        rs_cphd_t cphd;
        RS_CHECK_OK(rs_sim_scene(&cphd, tg, 2, 20.0, 400.0, 256, 0.5));

        rs_grid_t grid = { .n_x = 96, .n_y = 96, .dx = 1.0, .dy = 1.0, .height = 0.0 };
        grid.origin[0] = grid.origin[1] = grid.origin[2] = 0.0;

        /* The source paper's operating point, and its stated design rule: the
         * interval between successive sub-apertures must be well under the
         * vibration period. These give t_sap 0.655 s (aperture fraction 3.3%,
         * matching both their trials) and dt 0.1225 s against a 0.5 s period,
         * so dt/T = 0.24 where theirs is 0.25.
         *
         * The first version of this test used 24 looks at 0.5 overlap, which is
         * dt/T = 1.60 -- the sub-apertures were sampled more slowly than the
         * target vibrated, aliasing the motion away entirely. It is recorded
         * here because the resulting map looked plausible: the targets were
         * present and the numbers finite, and only the comparison between them
         * gave it away. */
        rs_subap_params_t sp;
        rs_subap_params_default(&sp);
        sp.n_looks = 156;
        sp.overlap = 0.81;
        sp.mode = RS_SUBAP_UNIFORM;

        rs_subap_stack_t stack;
        RS_CHECK_OK(rs_subaperture_from_cphd(&cphd, &grid, &sp, &stack));

        rs_ccd_params_t cp;
        rs_ccd_params_default(&cp);
        rs_ccd_t ccd;
        RS_CHECK_OK(rs_ccd_locate(&stack, &cp, &ccd));

        /* Targets sit at x = -20 and +20 m on a 1 m grid centred on the origin,
         * so they fall 20 cells either side of the middle row. Image row is the
         * grid's x. */
        const size_t mid_r = ccd.n_row / 2, mid_c = ccd.n_col / 2;
        const double vib = rs_ccd_patch_mean(&ccd, mid_r - 20, mid_c, 5);
        const double sta = rs_ccd_patch_mean(&ccd, mid_r + 20, mid_c, 5);

        double bg = 0.0;
        size_t nbg = 0;
        for (size_t r = 5; r < ccd.n_row - 5; r++) {
            for (size_t c = 5; c < ccd.n_col - 5; c++) {
                const size_t dr = r > mid_r ? r - mid_r : mid_r - r;
                if (dr < 30) continue;          /* skip both target rows */
                bg += ccd.map[r * ccd.n_col + c];
                nbg++;
            }
        }
        if (nbg) bg /= (double)nbg;

        printf("    vibrating %.4f   static %.4f   background %.4f   ratio v/s %.3f\n",
               vib, sta, bg, sta > 0.0 ? vib / sta : 0.0);

        /* The behaviour the detector exists for, asserted as a comparison on a
         * controlled fixture rather than as a threshold: two targets of equal
         * brightness, one vibrating, and the statistic must separate them.
         *
         * The margin is deliberately loose against the measured 2.76x. What is
         * NOT asserted is any absolute value -- the source paper implements no
         * detection threshold and offers no false-alarm rate, so there is no
         * published figure to pin, and inventing one here would encode a number
         * nobody has justified. */
        RS_CHECK(isfinite(vib) && isfinite(sta));
        RS_CHECK(vib > 1.5 * sta);

        /* The static target must also be close to the no-change value of 1,
         * which is the scale-invariance claim showing up on a scene: a bright
         * stationary scatterer is not merely quieter than the vibrating one, it
         * is quiet in absolute terms. */
        RS_CHECK(sta < 1.3);

        rs_ccd_free(&ccd);
        rs_subap_stack_free(&stack);
        rs_cphd_free(&cphd);
    }

    RS_TEST_END();
}
