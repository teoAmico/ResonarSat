/* Sub-pixel offset tracking: can it recover a known shift?
 *
 * Everything the micro-motion stage reports rests on this primitive locating a
 * correlation peak to a fraction of a pixel, so the tolerances here are the
 * tightest in the suite. */

#include "resonarsat/coreg.h"
#include "rs_test.h"

#include <math.h>
#include <stdlib.h>

/* Fill a patch with a Gaussian blob centred at (cy, cx). Sub-pixel centres are
 * representable because the Gaussian is evaluated analytically. */
static void make_patch(float complex *p, size_t n, double cy, double cx, double sigma)
{
    for (size_t a = 0; a < n; a++) {
        for (size_t r = 0; r < n; r++) {
            const double dy = ((double)a - cy) / sigma;
            const double dx = ((double)r - cx) / sigma;
            p[a * n + r] = (float)exp(-0.5 * (dy * dy + dx * dx));
        }
    }
}

/* Remove the mean, as rs_coreg_extract() does for real patches. */
static void demean(float complex *p, size_t n)
{
    float complex s = 0.0f;
    for (size_t i = 0; i < n * n; i++) s += p[i];
    const float complex m = s / (float)(n * n);
    for (size_t i = 0; i < n * n; i++) p[i] -= m;
}

int main(void)
{
    const size_t n = 32;
    const double sigma = 2.0;
    float complex *ref = malloc(n * n * sizeof *ref);
    float complex *cur = malloc(n * n * sizeof *cur);

    /* Integer shifts first: if these fail nothing else is meaningful. */
    RS_CASE("integer shifts recovered exactly");
    {
        const double shifts[][2] = { {0,0}, {1,0}, {0,1}, {2,-3}, {-2,2} };
        for (size_t i = 0; i < sizeof shifts / sizeof shifts[0]; i++) {
            make_patch(ref, n, 16.0, 16.0, sigma);
            make_patch(cur, n, 16.0 + shifts[i][0], 16.0 + shifts[i][1], sigma);
            demean(ref, n); demean(cur, n);

            double sa = 0, sr = 0, pk = 0;
            RS_CHECK_OK(rs_coreg_shift(ref, cur, n, n, 10, 20, &sa, &sr, &pk));
            RS_CHECK_NEAR(sa, shifts[i][0], 0.2);
            RS_CHECK_NEAR(sr, shifts[i][1], 0.2);
            RS_CHECK(pk > 0.5);
        }
    }

    /* Sub-pixel shifts are the real requirement: structural vibration moves a
     * scatterer by a small fraction of a resolution cell, so a tracker that
     * only resolves whole pixels measures nothing. */
    RS_CASE("sub-pixel shifts recovered to better than 0.15 px");
    {
        const double shifts[] = { 0.25, 0.5, -0.4, 1.3, -1.75 };
        for (size_t i = 0; i < sizeof shifts / sizeof shifts[0]; i++) {
            make_patch(ref, n, 16.0, 16.0, sigma);
            make_patch(cur, n, 16.0 + shifts[i], 16.0, sigma);
            demean(ref, n); demean(cur, n);

            double sa = 0, sr = 0, pk = 0;
            RS_CHECK_OK(rs_coreg_shift(ref, cur, n, n, 20, 20, &sa, &sr, &pk));
            RS_CHECK_NEAR(sa, shifts[i], 0.15);
            RS_CHECK_NEAR(sr, 0.0, 0.15);
        }
    }

    RS_CASE("identical patches correlate at unity");
    {
        make_patch(ref, n, 16.0, 16.0, sigma);
        demean(ref, n);
        double sa = 0, sr = 0, pk = 0;
        RS_CHECK_OK(rs_coreg_shift(ref, ref, n, n, 10, 10, &sa, &sr, &pk));
        RS_CHECK_NEAR(pk, 1.0, 0.05);
    }

    /* A blank patch must mask itself out rather than divide by zero. */
    RS_CASE("degenerate patches yield zero peak, not a crash");
    {
        for (size_t i = 0; i < n * n; i++) { ref[i] = 0.0f; cur[i] = 0.0f; }
        double sa = 1.0, sr = 1.0, pk = 1.0;
        RS_CHECK_OK(rs_coreg_shift(ref, cur, n, n, 10, 10, &sa, &sr, &pk));
        RS_CHECK_NEAR(pk, 0.0, 1e-9);
        RS_CHECK_NEAR(sa, 0.0, 1e-9);
    }

    RS_CASE("null arguments are refused");
    {
        double sa, sr, pk;
        RS_CHECK_ERR(rs_coreg_shift(NULL, cur, n, n, 1, 1, &sa, &sr, &pk), RS_ERR_ARG);
        RS_CHECK_ERR(rs_coreg_shift(ref, cur, 0, n, 1, 1, &sa, &sr, &pk), RS_ERR_ARG);
    }

    free(ref);
    free(cur);
    RS_TEST_END();
}
