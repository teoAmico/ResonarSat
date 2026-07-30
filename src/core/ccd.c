/* Scale-invariant coherent change detection over a sub-aperture stack.
 *
 * See include/resonarsat/ccd.h for what this stage is for and what a map from
 * it does and does not mean. */

#include "resonarsat/ccd.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* A 2x2 Hermitian covariance, stored by its independent entries.
 *
 * The diagonal of R R^H is real by construction; keeping it as a double rather
 * than taking creal() of a complex accumulator avoids carrying an imaginary
 * part that rounding makes non-zero and that every later expression would then
 * have to ignore. */
typedef struct {
    double a;              /* S[0][0], real */
    double d;              /* S[1][1], real */
    double complex b;      /* S[0][1]; S[1][0] is its conjugate */
} rs_ccd_herm2_t;

/* Accumulate S = R R^H over one window, for the channel pair (c0, c1).
 *
 * 'c0' and 'c1' are whole sub-look images on a shared grid; the window is the
 * square of side 'win' whose top-left corner is (row0, col0). Each pixel of the
 * window contributes one column of R, so S is the sum over the window of the
 * outer product of the per-pixel channel vector with itself. */
static rs_ccd_herm2_t rs_ccd_cov(const float complex *c0, const float complex *c1,
                                 size_t n_col, size_t row0, size_t col0, size_t win)
{
    rs_ccd_herm2_t s = { 0.0, 0.0, 0.0 };
    for (size_t r = 0; r < win; r++) {
        const size_t base = (row0 + r) * n_col + col0;
        for (size_t c = 0; c < win; c++) {
            const double complex x0 = c0[base + c];
            const double complex x1 = c1[base + c];
            s.a += creal(x0) * creal(x0) + cimag(x0) * cimag(x0);
            s.d += creal(x1) * creal(x1) + cimag(x1) * cimag(x1);
            s.b += x0 * conj(x1);
        }
    }
    return s;
}

/* The ratio of the eigenvalues of S_X S_Y^-1, largest over smallest.
 *
 * Both matrices are Hermitian positive semi-definite, so the eigenvalues of the
 * product are real and non-negative and the quadratic below has a non-negative
 * discriminant in exact arithmetic. It is clamped anyway, because a window of
 * nearly collinear samples can drive the discriminant a few ulp below zero.
 *
 * 'load' is an ABSOLUTE diagonal term, in the units of S's diagonal, added to
 * both matrices before the ratio is formed. It must not be a fraction of each
 * window's own trace: that scales away exactly where it is needed. A window
 * containing no scatterers has a trace near zero, its covariance is noise, and
 * a trace-relative load leaves the ratio free to take any value at all -- so the
 * detector lights up brightest over empty ground, which is the precise opposite
 * of its purpose. Measured on the synthetic two-target scene before this was
 * fixed: background 75612 against 1.5 at the targets. An absolute floor derived
 * from the scene's own mean power makes an empty window resolve to load*I on
 * both sides, hence a ratio of 1, which is the null value.
 *
 * Adding it to both matrices rather than to S_Y alone keeps the scaling
 * identity intact: if the whole stack is scaled, the floor scales with it and
 * the ratio is unchanged.
 *
 * Returns 1.0 -- the no-change value -- when S_Y is degenerate even after
 * loading. */
static double rs_ccd_eig_ratio(rs_ccd_herm2_t sx, rs_ccd_herm2_t sy, double load)
{
    if (load > 0.0) {
        sx.a += load;
        sx.d += load;
        sy.a += load;
        sy.d += load;
    }

    const double det_y = sy.a * sy.d - (creal(sy.b) * creal(sy.b) + cimag(sy.b) * cimag(sy.b));
    const double det_x = sx.a * sx.d - (creal(sx.b) * creal(sx.b) + cimag(sx.b) * cimag(sx.b));
    if (!(det_y > 0.0) || !isfinite(det_y)) return 1.0;

    /* A = S_X S_Y^-1 with S_Y^-1 = adj(S_Y)/det(S_Y). Only the trace and the
     * determinant of A are needed, and det(A) = det(S_X)/det(S_Y) directly. */
    const double complex adj01 = -sy.b;
    const double complex adj10 = -conj(sy.b);
    const double complex a00 = sx.a * sy.d + sx.b * adj10;
    const double complex a11 = conj(sx.b) * adj01 + sx.d * sy.a;
    const double tr = (creal(a00) + creal(a11)) / det_y;
    const double det = det_x / det_y;

    double disc = tr * tr - 4.0 * det;
    if (disc < 0.0) disc = 0.0;
    const double root = sqrt(disc);
    const double l1 = 0.5 * (tr + root);
    const double l2 = 0.5 * (tr - root);

    if (!(l2 > 0.0) || !isfinite(l1) || !isfinite(l2)) return 1.0;
    return l1 / l2;
}

/* Fill 'params' with the source paper's window and a conservative loading. */
void rs_ccd_params_default(rs_ccd_params_t *params)
{
    if (!params) return;
    params->win = 5;              /* the paper's 5x5 sliding window */
    params->stat = RS_CCD_G12;
    params->loading = 1e-3;
}

/* Release a map and everything it owns. */
void rs_ccd_free(rs_ccd_t *ccd)
{
    if (!ccd) return;
    free(ccd->map);
    memset(ccd, 0, sizeof *ccd);
}

/* Accumulate the change statistic over every sub-aperture triple in 'stack'. */
resonarsat_status_t rs_ccd_locate(const rs_subap_stack_t *stack,
                                  const rs_ccd_params_t *params,
                                  rs_ccd_t *out)
{
    if (!out) return RS_ERR_ARG;
    memset(out, 0, sizeof *out);
    if (!stack || !params || !stack->look) return RS_ERR_ARG;

    if (stack->n_looks < 3) {
        rs_set_error("ccd: %zu looks is too few; a triple needs three",
                     stack->n_looks);
        return RS_ERR_ARG;
    }
    if (params->win == 0) {
        rs_set_error("ccd: window must be positive");
        return RS_ERR_ARG;
    }
    if (params->stat != RS_CCD_G12) {
        rs_set_error("ccd: only the two-channel statistic is implemented");
        return RS_ERR_ARG;
    }

    const size_t n_row = stack->look[0].n_az;
    const size_t n_col = stack->look[0].n_rg;
    const size_t win = params->win;
    if (win > n_row || win > n_col) {
        rs_set_error("ccd: %zux%zu window does not fit a %zux%zu image",
                     win, win, n_row, n_col);
        return RS_ERR_ARG;
    }

    double *map = calloc(n_row * n_col, sizeof *map);
    if (!map) return RS_ERR_ALLOC;

    const size_t half = win / 2;
    const size_t n_triples = stack->n_looks - 2;

    /* The absolute noise floor for the diagonal loading, derived from the
     * scene's own mean power so that it scales with the data rather than
     * assuming a calibration. Multiplied by win*win because S accumulates over
     * the window, so its diagonal is that many pixel powers. See
     * rs_ccd_eig_ratio() for why a window-relative load does not work. */
    double ref_power = 0.0;
    {
        double sum = 0.0;
        const size_t npx = n_row * n_col;
        for (size_t i = 0; i < stack->n_looks; i++) {
            const float complex *d = stack->look[i].data;
            for (size_t p = 0; p < npx; p++) {
                sum += creal(d[p]) * creal(d[p]) + cimag(d[p]) * cimag(d[p]);
            }
        }
        ref_power = sum / ((double)npx * (double)stack->n_looks);
    }
    const double load = params->loading * ref_power * (double)(win * win);

    /* Every window is independent, and the accumulator is private per pixel, so
     * the row loop parallelises without a reduction. The triple loop is inside
     * so that each thread walks one row of the image across all looks rather
     * than the whole image once per look, which keeps the three sub-look rows it
     * touches resident. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t r = half; r < n_row - half; r++) {
        for (size_t c = half; c < n_col - half; c++) {
            const size_t row0 = r - half, col0 = c - half;
            double acc = 0.0;
            for (size_t m = 1; m + 1 < stack->n_looks; m++) {
                const float complex *lm1 = stack->look[m - 1].data;
                const float complex *lm  = stack->look[m].data;
                const float complex *lp1 = stack->look[m + 1].data;
                const rs_ccd_herm2_t sx = rs_ccd_cov(lm1, lm, n_col, row0, col0, win);
                const rs_ccd_herm2_t sy = rs_ccd_cov(lm, lp1, n_col, row0, col0, win);
                acc += rs_ccd_eig_ratio(sx, sy, load);
            }
            map[r * n_col + c] = acc / (double)n_triples;
        }
    }

    out->map = map;
    out->n_row = n_row;
    out->n_col = n_col;
    out->n_triples = n_triples;
    out->params = *params;
    return RS_OK;
}
