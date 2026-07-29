/* FFT wrapper implementation -- the only file in ResonarSat that names the
 * FFT backend. See include/resonarsat/fft.h for the rationale.
 *
 * Backend: PocketFFT (BSD-3-Clause), vendored under third_party/pocketfft.
 *
 * One backend property leaks into the implementation and is worth stating
 * plainly: the C build of PocketFFT transforms double-precision interleaved
 * data only, while this pipeline carries single-precision complex samples.
 * Every transform therefore widens to double on entry and narrows on exit,
 * through a scratch buffer owned by the plan. That costs memory bandwidth and
 * an allocation per plan, and it is a deliberate trade: correctness and a
 * clean licence now, with the option of a float-native backend later behind an
 * unchanged interface. */

#include "resonarsat/fft.h"

#include <stdlib.h>
#include <string.h>

#include "pocketfft.h"

struct rs_fft_plan {
    cfft_plan backend;  /* PocketFFT handle */
    size_t    n;        /* transform length */
    double   *scratch;  /* 2*n doubles, interleaved re/im */
};

/* Create a reusable plan for transforms of length n. */
resonarsat_status_t rs_fft_plan_create(size_t n, rs_fft_plan **out)
{
    if (!out) return RS_ERR_ARG;
    *out = NULL;
    if (n == 0) {
        rs_set_error("fft: transform length must be positive");
        return RS_ERR_ARG;
    }

    rs_fft_plan *p = calloc(1, sizeof *p);
    if (!p) return RS_ERR_ALLOC;

    p->n = n;
    p->scratch = malloc(2 * n * sizeof *p->scratch);
    if (!p->scratch) {
        free(p);
        return RS_ERR_ALLOC;
    }

    p->backend = make_cfft_plan(n);
    if (!p->backend) {
        free(p->scratch);
        free(p);
        rs_set_error("fft: backend rejected length %zu", n);
        return RS_ERR_ALLOC;
    }

    *out = p;
    return RS_OK;
}

/* Release a plan; accepts NULL. */
void rs_fft_plan_destroy(rs_fft_plan *plan)
{
    if (!plan) return;
    if (plan->backend) destroy_cfft_plan(plan->backend);
    free(plan->scratch);
    free(plan);
}

/* Transform length this plan was built for. */
size_t rs_fft_plan_length(const rs_fft_plan *plan)
{
    return plan ? plan->n : 0;
}

/* Widen caller samples into the plan's double scratch buffer. Shared by the
 * forward and inverse paths so the two cannot drift apart. */
static void rs_fft_widen(const rs_fft_plan *plan, const float complex *src)
{
    for (size_t i = 0; i < plan->n; i++) {
        plan->scratch[2 * i]     = (double)crealf(src[i]);
        plan->scratch[2 * i + 1] = (double)cimagf(src[i]);
    }
}

/* Narrow the plan's scratch buffer back into caller samples, applying a scale
 * factor as it goes (1.0 forward, 1/n inverse). */
static void rs_fft_narrow(const rs_fft_plan *plan, float complex *dst, double scale)
{
    for (size_t i = 0; i < plan->n; i++) {
        dst[i] = (float)(plan->scratch[2 * i] * scale)
               + (float)(plan->scratch[2 * i + 1] * scale) * I;
    }
}

/* Forward transform, unnormalised. */
resonarsat_status_t rs_fft_forward(rs_fft_plan *plan, float complex *data)
{
    if (!plan || !data) return RS_ERR_ARG;
    rs_fft_widen(plan, data);
    if (cfft_forward(plan->backend, plan->scratch, 1.0) != 0) {
        rs_set_error("fft: forward transform failed at length %zu", plan->n);
        return RS_ERR_SINGULAR;
    }
    rs_fft_narrow(plan, data, 1.0);
    return RS_OK;
}

/* Inverse transform, normalised by 1/n so that forward-then-inverse is the
 * identity. The backend applies the scale itself, so the narrowing pass uses
 * unit scale. */
resonarsat_status_t rs_fft_inverse(rs_fft_plan *plan, float complex *data)
{
    if (!plan || !data) return RS_ERR_ARG;
    rs_fft_widen(plan, data);
    if (cfft_backward(plan->backend, plan->scratch, 1.0 / (double)plan->n) != 0) {
        rs_set_error("fft: inverse transform failed at length %zu", plan->n);
        return RS_ERR_SINGULAR;
    }
    rs_fft_narrow(plan, data, 1.0);
    return RS_OK;
}

/* Move the zero-frequency component from index 0 to the centre.
 *
 * Implemented as a rotation left by n/2 via three reversals, which needs no
 * scratch buffer and handles odd n correctly. For even n the operation is an
 * involution; for odd n it is not, which is why rs_ifft_shift() exists. */
static void rs_reverse(float complex *a, size_t lo, size_t hi)
{
    while (lo < hi) {
        float complex t = a[lo];
        a[lo++] = a[hi];
        a[hi--] = t;
    }
}

/* Rotate an array left by k positions using the three-reversal identity
 * reverse(0,k-1), reverse(k,n-1), reverse(0,n-1). No scratch buffer, and it
 * handles any k including k >= n. */
static void rs_rotate_left(float complex *data, size_t n, size_t k)
{
    if (n == 0) return;
    k %= n;
    if (k == 0) return;
    rs_reverse(data, 0, k - 1);
    rs_reverse(data, k, n - 1);
    rs_reverse(data, 0, n - 1);
}

/* Zero frequency to the centre. */
void rs_fft_shift(float complex *data, size_t n)
{
    if (!data || n < 2) return;
    rs_rotate_left(data, n, (n + 1) / 2);
}

/* Centre back to index 0; exact inverse of rs_fft_shift() for all n. */
void rs_ifft_shift(float complex *data, size_t n)
{
    if (!data || n < 2) return;
    rs_rotate_left(data, n, n - (n + 1) / 2);
}

/* Two-dimensional transform, rows then columns.
 *
 * Columns are gathered into a contiguous scratch buffer before transforming
 * because the backend needs contiguous interleaved input; a strided column
 * transform would require a backend-specific API this wrapper deliberately does
 * not expose. The gather costs one extra pass over the array and keeps the
 * interface clean.
 *
 * The inverse direction is normalised once per axis, giving 1/(n_row*n_col)
 * overall, so a round trip through this function is the identity. */
resonarsat_status_t rs_fft2(float complex *data, size_t n_row, size_t n_col, int inverse)
{
    if (!data) return RS_ERR_ARG;
    if (n_row == 0 || n_col == 0) {
        rs_set_error("fft2: dimensions must be positive (got %zux%zu)", n_row, n_col);
        return RS_ERR_ARG;
    }

    resonarsat_status_t st;
    rs_fft_plan *row_plan = NULL, *col_plan = NULL;
    float complex *col = NULL;

    if ((st = rs_fft_plan_create(n_col, &row_plan)) != RS_OK) goto done;
    if ((st = rs_fft_plan_create(n_row, &col_plan)) != RS_OK) goto done;

    col = malloc(n_row * sizeof *col);
    if (!col) { st = RS_ERR_ALLOC; goto done; }

    for (size_t r = 0; r < n_row; r++) {
        st = inverse ? rs_fft_inverse(row_plan, data + r * n_col)
                     : rs_fft_forward(row_plan, data + r * n_col);
        if (st != RS_OK) goto done;
    }

    for (size_t c = 0; c < n_col; c++) {
        for (size_t r = 0; r < n_row; r++) col[r] = data[r * n_col + c];
        st = inverse ? rs_fft_inverse(col_plan, col)
                     : rs_fft_forward(col_plan, col);
        if (st != RS_OK) goto done;
        for (size_t r = 0; r < n_row; r++) data[r * n_col + c] = col[r];
    }

done:
    free(col);
    rs_fft_plan_destroy(row_plan);
    rs_fft_plan_destroy(col_plan);
    return st;
}
