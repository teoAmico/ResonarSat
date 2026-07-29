/* FFT wrapper: round trips, normalisation, and shift correctness. */

#include "resonarsat/fft.h"
#include "rs_test.h"

#include <complex.h>
#include <stdlib.h>

/* A forward transform followed by an inverse must return the input, for both
 * power-of-two and awkward composite lengths. */
static void test_round_trip(size_t n)
{
    float complex *a = malloc(n * sizeof *a);
    float complex *orig = malloc(n * sizeof *orig);
    rs_fft_plan *plan = NULL;

    RS_CHECK_OK(rs_fft_plan_create(n, &plan));
    for (size_t i = 0; i < n; i++) {
        a[i] = (float)(i % 7) - 3.0f + (float)((i * 3) % 5) * I;
        orig[i] = a[i];
    }

    RS_CHECK_OK(rs_fft_forward(plan, a));
    RS_CHECK_OK(rs_fft_inverse(plan, a));

    double worst = 0.0;
    for (size_t i = 0; i < n; i++) {
        const double e = (double)cabsf(a[i] - orig[i]);
        if (e > worst) worst = e;
    }
    RS_CHECK_NEAR(worst, 0.0, 1e-4);

    rs_fft_plan_destroy(plan);
    free(a);
    free(orig);
}

int main(void)
{
    RS_CASE("round trip, power of two");
    test_round_trip(64);
    test_round_trip(1024);

    RS_CASE("round trip, non power of two");
    test_round_trip(60);
    test_round_trip(97);   /* prime */

    /* A pure tone must land in exactly one bin, at the expected index and with
     * the expected magnitude. This pins the sign and normalisation convention;
     * a transform that is correct up to a conjugate would pass a round-trip
     * test and fail here. */
    RS_CASE("single tone lands in one bin");
    {
        const size_t n = 64;
        const size_t k0 = 5;
        float complex *a = malloc(n * sizeof *a);
        rs_fft_plan *plan = NULL;
        RS_CHECK_OK(rs_fft_plan_create(n, &plan));

        for (size_t i = 0; i < n; i++) {
            const double ph = 2.0 * M_PI * (double)k0 * (double)i / (double)n;
            a[i] = (float)cos(ph) + (float)sin(ph) * I;
        }
        RS_CHECK_OK(rs_fft_forward(plan, a));

        size_t peak = 0;
        double peak_mag = 0.0;
        for (size_t i = 0; i < n; i++) {
            const double m = (double)cabsf(a[i]);
            if (m > peak_mag) { peak_mag = m; peak = i; }
        }
        RS_CHECK(peak == k0);
        RS_CHECK_REL(peak_mag, (double)n, 1e-3);

        rs_fft_plan_destroy(plan);
        free(a);
    }

    /* fftshift then ifftshift must be the identity for both parities. Odd
     * lengths are the case where a naive halves-swap implementation breaks. */
    RS_CASE("shift and inverse shift round trip");
    for (size_t n = 2; n <= 9; n++) {
        float complex *a = malloc(n * sizeof *a);
        for (size_t i = 0; i < n; i++) a[i] = (float)i;
        rs_fft_shift(a, n);
        rs_ifft_shift(a, n);
        for (size_t i = 0; i < n; i++) RS_CHECK_NEAR(crealf(a[i]), (double)i, 1e-6);
        free(a);
    }

    RS_CASE("shift moves DC to the centre");
    {
        const size_t n = 8;
        float complex a[8] = {0};
        a[0] = 1.0f;  /* DC */
        rs_fft_shift(a, n);
        RS_CHECK_NEAR(crealf(a[n / 2]), 1.0, 1e-6);
    }

    RS_CASE("2-D round trip");
    {
        const size_t nr = 16, nc = 24;
        float complex *a = malloc(nr * nc * sizeof *a);
        float complex *o = malloc(nr * nc * sizeof *o);
        for (size_t i = 0; i < nr * nc; i++) { a[i] = (float)(i % 11); o[i] = a[i]; }
        RS_CHECK_OK(rs_fft2(a, nr, nc, 0));
        RS_CHECK_OK(rs_fft2(a, nr, nc, 1));
        double worst = 0.0;
        for (size_t i = 0; i < nr * nc; i++) {
            const double e = (double)cabsf(a[i] - o[i]);
            if (e > worst) worst = e;
        }
        RS_CHECK_NEAR(worst, 0.0, 1e-3);
        free(a); free(o);
    }

    RS_CASE("invalid arguments are refused");
    {
        rs_fft_plan *p = NULL;
        RS_CHECK_ERR(rs_fft_plan_create(0, &p), RS_ERR_ARG);
        RS_CHECK(p == NULL);
        RS_CHECK_ERR(rs_fft_forward(NULL, NULL), RS_ERR_ARG);
        rs_fft_plan_destroy(NULL);   /* must not crash */
    }

    RS_TEST_END();
}
