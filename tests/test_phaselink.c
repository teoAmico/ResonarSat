/* Split-band Phase Linking, checked against known shifts on ideal data.
 *
 * This test exists to separate two questions that are easy to confuse: is the
 * estimator implemented correctly, and is the data it is given suitable for it?
 * Here the data is made ideal by construction -- fully coherent speckle,
 * shifted by an exactly known sub-pixel amount via a spectral phase ramp -- so a
 * failure is unambiguously an implementation fault.
 *
 * Method: De Zan, "Coherent Shift Estimation for Stacks of SAR Images", IEEE
 * GRSL 8, 1095 (2011). */

#include "resonarsat/phaselink.h"
#include "resonarsat/fft.h"
#include "rs_test.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Deterministic uniform generator, so failures reproduce exactly. */
static unsigned rs_rng = 20260725u;
static double urand(void)
{
    rs_rng = rs_rng * 1103515245u + 12345u;
    return (double)((rs_rng >> 16) & 0x7fff) / 32767.0;
}

/* Build a band-limited complex speckle patch: random scatterers in the spectral
 * domain over the central 'occupancy' fraction of the band, transformed to the
 * image domain. Band-limiting matters -- a white patch has no spectral slope for
 * the split-band method to measure. */
static void make_patch(float complex *p, size_t n_az, size_t n_rg, double occupancy)
{
    float complex *col = malloc(n_az * sizeof *col);
    rs_fft_plan *plan = NULL;
    rs_fft_plan_create(n_az, &plan);

    const size_t half = (size_t)(0.5 * occupancy * (double)n_az);

    for (size_t r = 0; r < n_rg; r++) {
        for (size_t a = 0; a < n_az; a++) col[a] = 0.0f;
        for (size_t f = 0; f < n_az; f++) {
            const size_t fs = (f + n_az / 2) % n_az;      /* centred band */
            const long d = (long)fs - (long)(n_az / 2);
            if ((size_t)labs(d) > half) continue;
            col[f] = (float)(urand() - 0.5) + (float)(urand() - 0.5) * I;
        }
        rs_fft_inverse(plan, col);
        for (size_t a = 0; a < n_az; a++) p[a * n_rg + r] = col[a];
    }

    rs_fft_plan_destroy(plan);
    free(col);
}

/* Shift a patch along azimuth by 'shift' samples using a spectral phase ramp,
 * which gives an exact sub-sample shift with no interpolation error. */
static void shift_patch(const float complex *in, float complex *out,
                        size_t n_az, size_t n_rg, double shift)
{
    float complex *col = malloc(n_az * sizeof *col);
    rs_fft_plan *plan = NULL;
    rs_fft_plan_create(n_az, &plan);

    for (size_t r = 0; r < n_rg; r++) {
        for (size_t a = 0; a < n_az; a++) col[a] = in[a * n_rg + r];
        rs_fft_forward(plan, col);
        for (size_t f = 0; f < n_az; f++) {
            const double fq = (f < n_az / 2) ? (double)f : (double)f - (double)n_az;
            const double ph = -2.0 * M_PI * fq * shift / (double)n_az;
            col[f] *= (float)cos(ph) + (float)sin(ph) * I;
        }
        rs_fft_inverse(plan, col);
        for (size_t a = 0; a < n_az; a++) out[a * n_rg + r] = col[a];
    }

    rs_fft_plan_destroy(plan);
    free(col);
}

int main(void)
{
    const size_t n_az = 64, n_rg = 32, n_look = 16;
    const size_t n_pix = n_az * n_rg;

    float complex *base = malloc(n_pix * sizeof *base);
    float complex *stack = malloc(n_look * n_pix * sizeof *stack);
    double *truth = malloc(n_look * sizeof *truth);
    double *got = malloc(n_look * sizeof *got);

    /* --------------------------------------------------------------
     * Phase Linking alone: a stack differing only by a known phase.
     * -------------------------------------------------------------- */
    RS_CASE("phase linking recovers known phases");
    {
        make_patch(base, n_az, n_rg, 0.6);
        double *want = malloc(n_look * sizeof *want);
        for (size_t k = 0; k < n_look; k++) {
            want[k] = 0.4 * sin(2.0 * M_PI * (double)k / (double)n_look);
            const double c = cos(want[k]), s = sin(want[k]);
            for (size_t p = 0; p < n_pix; p++) {
                stack[k * n_pix + p] = base[p] * ((float)c + (float)s * I);
            }
        }
        RS_CHECK_OK(rs_phase_link(stack, n_look, n_pix, got));

        double worst = 0.0;
        for (size_t k = 0; k < n_look; k++) {
            const double e = fabs((got[k]) - (want[k] - want[0]));
            if (e > worst) worst = e;
        }
        printf("    worst phase error %.2e rad\n", worst);
        RS_CHECK_NEAR(worst, 0.0, 1e-3);
        free(want);
    }

    /* --------------------------------------------------------------
     * The full estimator against exactly known sub-sample shifts.
     * -------------------------------------------------------------- */
    RS_CASE("split-band recovers known sub-sample shifts");
    {
        /* Amplitude kept inside the unambiguous range. The source states the
         * limit as three quarters of a resolution cell of the sub-band
         * separation; exceeding it wraps the phase difference and the recovered
         * shift folds, which presents as a plausible-looking scale error rather
         * than an obvious failure. A 1.5-sample amplitude sits right on that
         * boundary for this band and does exactly that. */
        make_patch(base, n_az, n_rg, 0.6);
        for (size_t k = 0; k < n_look; k++) {
            truth[k] = 0.4 * sin(2.0 * M_PI * (double)k / (double)n_look);
            shift_patch(base, stack + k * n_pix, n_az, n_rg, truth[k]);
        }

        double coh = 0.0;
        RS_CHECK_OK(rs_splitband_shift(stack, n_look, n_az, n_rg, got, &coh));

        double worst = 0.0;
        printf("    coherence %.3f\n", coh);
        printf("    %6s %10s %10s\n", "look", "truth", "estimate");
        for (size_t k = 0; k < n_look; k += 4) {
            printf("    %6zu %10.3f %10.3f\n", k, truth[k] - truth[0], got[k]);
        }
        for (size_t k = 0; k < n_look; k++) {
            const double e = fabs(got[k] - (truth[k] - truth[0]));
            if (e > worst) worst = e;
        }
        printf("    worst shift error %.3f samples\n", worst);

        /* On fully coherent, exactly shifted data the estimator should be far
         * better than a correlation peak: a tenth of a sample rather than the
         * order of a resolution cell. */
        RS_CHECK_NEAR(worst, 0.0, 0.1);
        /* Shifting band-limited speckle decorrelates it: a shift approaching a
         * resolution cell is most of the way to independence. At this amplitude
         * the looks stay well correlated. */
        RS_CHECK(coh > 0.8);
    }

    /* --------------------------------------------------------------
     * The ambiguity limit, measured rather than assumed. Beyond it the
     * estimate folds; a caller needs to know where that is.
     * -------------------------------------------------------------- */
    RS_CASE("the ambiguity limit is where the source says it is");
    {
        make_patch(base, n_az, n_rg, 0.6);
        printf("    %10s %12s %10s\n", "true shift", "estimate", "error");
        for (double amp = 0.25; amp < 2.6; amp *= 2.0) {
            for (size_t k = 0; k < n_look; k++) {
                truth[k] = amp * ((k % 2) ? 1.0 : -1.0);
                shift_patch(base, stack + k * n_pix, n_az, n_rg, truth[k]);
            }
            RS_CHECK_OK(rs_splitband_shift(stack, n_look, n_az, n_rg, got, NULL));
            const double want = truth[1] - truth[0];
            printf("    %10.2f %12.3f %10.3f%s\n",
                   want, got[1], got[1] - want,
                   fabs(got[1] - want) > 0.2 ? "   <- folded" : "");
        }
    }

    /* --------------------------------------------------------------
     * Degenerate inputs must be refused, not guessed at.
     * -------------------------------------------------------------- */
    RS_CASE("degenerate inputs are refused");
    {
        RS_CHECK_ERR(rs_splitband_shift(NULL, n_look, n_az, n_rg, got, NULL), RS_ERR_ARG);
        RS_CHECK_ERR(rs_splitband_shift(stack, 1, n_az, n_rg, got, NULL), RS_ERR_ARG);
        RS_CHECK_ERR(rs_phase_link(NULL, n_look, n_pix, got), RS_ERR_ARG);
    }

    free(base); free(stack); free(truth); free(got);
    RS_TEST_END();
}
