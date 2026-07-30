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

    /* --no-optimize's audit path. The whole value of a baseline is that it agrees
     * with the fast path on cases where the fast path is known to be right, so
     * that a disagreement elsewhere means something. Tested to one lattice step
     * rather than exactly: the two search the same lattice but arrive by
     * different arithmetic -- a padded transform against a direct summation --
     * and single-precision rounding can put the crest one step either way when
     * two neighbouring points are within an ulp of each other. */
    RS_CASE("exhaustive search agrees with local refinement on clean shifts");
    {
        const size_t up = 20;
        const double tol = 1.5 / (double)up;
        const double shifts[][2] = { {0,0}, {1,0}, {0,1}, {2,-3}, {0.25,0},
                                     {-0.4,0.5}, {1.3,-1.75} };
        for (size_t i = 0; i < sizeof shifts / sizeof shifts[0]; i++) {
            make_patch(ref, n, 16.0, 16.0, sigma);
            make_patch(cur, n, 16.0 + shifts[i][0], 16.0 + shifts[i][1], sigma);
            demean(ref, n); demean(cur, n);

            double la = 0, lr = 0, lpk = 0, xa = 0, xr = 0, xpk = 0;
            RS_CHECK_OK(rs_coreg_shift_ex(ref, cur, n, n, up, up,
                                          RS_COREG_REFINE_LOCAL, &la, &lr, &lpk));
            RS_CHECK_OK(rs_coreg_shift_ex(ref, cur, n, n, up, up,
                                          RS_COREG_REFINE_EXHAUSTIVE, &xa, &xr, &xpk));

            /* Against the truth, so a bug that moved both paths equally cannot
             * pass by having them agree with each other. */
            RS_CHECK_NEAR(xa, shifts[i][0], 0.15);
            RS_CHECK_NEAR(xr, shifts[i][1], 0.15);
            RS_CHECK_NEAR(xa, la, tol);
            RS_CHECK_NEAR(xr, lr, tol);
            /* The normalisations must match too, or the coherence mask would
             * threshold two different quantities depending on the mode. */
            RS_CHECK_NEAR(xpk, lpk, 0.02);
        }
    }

    /* The exhaustive path must land on the same 1/upsample lattice as the local
     * one. A padding or index-unwrap error would still produce plausible shifts,
     * just off the lattice -- which is invisible in a tolerance check. */
    RS_CASE("exhaustive shifts land on the 1/upsample lattice");
    {
        const size_t up = 8;
        make_patch(ref, n, 16.0, 16.0, sigma);
        make_patch(cur, n, 16.0 + 0.375, 16.0 - 1.25, sigma);
        demean(ref, n); demean(cur, n);

        double sa = 0, sr = 0, pk = 0;
        RS_CHECK_OK(rs_coreg_shift_ex(ref, cur, n, n, up, up,
                                      RS_COREG_REFINE_EXHAUSTIVE, &sa, &sr, &pk));
        RS_CHECK_NEAR(sa * (double)up, round(sa * (double)up), 1e-9);
        RS_CHECK_NEAR(sr * (double)up, round(sr * (double)up), 1e-9);
    }

    /* A degenerate patch must mask itself out on this path too. The exhaustive
     * branch returns before the search, so this checks the early exit is shared
     * rather than duplicated and forgotten. */
    RS_CASE("exhaustive search masks out a blank patch");
    {
        for (size_t i = 0; i < n * n; i++) { ref[i] = 0.0f; cur[i] = 0.0f; }
        double sa = 1.0, sr = 1.0, pk = 1.0;
        RS_CHECK_OK(rs_coreg_shift_ex(ref, cur, n, n, 10, 10,
                                      RS_COREG_REFINE_EXHAUSTIVE, &sa, &sr, &pk));
        RS_CHECK_NEAR(pk, 0.0, 1e-9);
        RS_CHECK_NEAR(sa, 0.0, 1e-9);
        RS_CHECK_NEAR(sr, 0.0, 1e-9);
    }

    /* An impossible surface must be refused with a described error rather than
     * attempted. The failure mode this guards is a tracker that returns a
     * complete result with every window zero. */
    RS_CASE("an oversized exhaustive surface is refused");
    {
        RS_CHECK_ERR(rs_coreg_surface_check(4096, 4096, 64, 64), RS_ERR_RANGE);
        RS_CHECK_OK(rs_coreg_surface_check(32, 32, 20, 20));

        double sa, sr, pk;
        RS_CHECK_ERR(rs_coreg_shift_ex(ref, cur, n, n, 4096, 4096,
                                       RS_COREG_REFINE_EXHAUSTIVE, &sa, &sr, &pk),
                     RS_ERR_RANGE);
    }

    free(ref);
    free(cur);
    RS_TEST_END();
}
