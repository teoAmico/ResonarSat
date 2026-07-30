/* Tomographic focusing. See include/resonarsat/tomo.h for the contract and for
 * the two caveats that govern how any output of this file may be presented. */

#include "resonarsat/tomo.h"
#include "resonarsat/fft.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Defaults for everything except the two required assumptions. */
void rs_tomo_params_default(rs_tomo_params_t *params)
{
    if (!params) return;
    memset(params, 0, sizeof *params);
    params->model = RS_TOMO_MODEL_A;
    params->solver = RS_TOMO_SOLVER_DFT;
    params->y_source = RS_TOMO_Y_SHIFTS;
    params->convention = RS_WAVELEN_PAPER;
    params->velocity = 0.0;    /* required from the caller */
    params->frequency = 0.0;   /* required from the caller */
    params->depth_max = 60.0;
    params->depth_cell = 1.0;
    params->slant_range = 500000.0;
    params->incidence = 35.0 * M_PI / 180.0;
    params->aperture = 75000.0;
    params->radar_wavelength = 0.031;   /* X-band */
    params->regularisation = 1e-3;
    params->eq22_literal_t = 0.0;       /* conventional exp(j*Kz*z) */
    params->window = 1;
    params->remove_y_mean = 1;
}

/* Scale applied to Kz by the experimental rendered-Eq.-22 mode. */
static double rs_tomo_eq22_scale(const rs_tomo_params_t *params)
{
    return (params && params->eq22_literal_t > 0.0)
         ? 2.0 * M_PI * params->eq22_literal_t : 1.0;
}

/* Validate a parameter set. */
resonarsat_status_t rs_tomo_params_check(const rs_tomo_params_t *params)
{
    if (!params) return RS_ERR_ARG;

    if (!(params->velocity > 0.0)) {
        rs_set_error("tomo: --velocity is required and has no default; it is an "
                     "assumed material wave speed that scales the depth axis");
        return RS_ERR_ARG;
    }
    /* Model D reads its frequency from the measured spectrum rather than taking
     * one on faith -- that is the whole point of it, and the one respect in
     * which it assumes less than the published method. Demanding --frequency
     * anyway forced callers to supply a number the focusing never used, and
     * then reported it in the sidecar as though it had. */
    if (params->model != RS_TOMO_MODEL_D && !(params->frequency > 0.0)) {
        rs_set_error("tomo: --frequency is required and has no default; it is the "
                     "assumed investigation frequency, NOT a measured vibration "
                     "frequency, and it scales the depth axis");
        return RS_ERR_ARG;
    }
    if (params->velocity < 100.0 || params->velocity > 20000.0) {
        rs_set_error("tomo: assumed wave speed %g m/s outside [100, 20000]",
                     params->velocity);
        return RS_ERR_RANGE;
    }
    if (!(params->depth_max > 0.0) || !(params->depth_cell > 0.0)) {
        rs_set_error("tomo: depth extent %g m and cell %g m must both be positive",
                     params->depth_max, params->depth_cell);
        return RS_ERR_ARG;
    }
    if (params->depth_cell > params->depth_max) {
        rs_set_error("tomo: depth cell %g m exceeds the depth extent %g m",
                     params->depth_cell, params->depth_max);
        return RS_ERR_ARG;
    }
    if (!(params->slant_range > 0.0) || !(params->aperture > 0.0)) {
        rs_set_error("tomo: slant range %g m and aperture %g m must be positive",
                     params->slant_range, params->aperture);
        return RS_ERR_ARG;
    }
    if (!(params->incidence > 0.0 && params->incidence < M_PI / 2.0)) {
        rs_set_error("tomo: incidence angle %g rad outside (0, pi/2)", params->incidence);
        return RS_ERR_RANGE;
    }
    if (!isfinite(params->eq22_literal_t) || params->eq22_literal_t < 0.0) {
        rs_set_error("tomo: --eq22-literal-t must be finite and positive when "
                     "provided (got %g)", params->eq22_literal_t);
        return RS_ERR_ARG;
    }
    if (params->eq22_literal_t > 0.0 &&
        (params->model != RS_TOMO_MODEL_A ||
         params->solver != RS_TOMO_SOLVER_LSTSQ)) {
        rs_set_error("tomo: --eq22-literal-t is experimental and requires "
                     "--model A --solver lstsq; the DFT shortcut implements "
                     "only conventional exp(j*Kz*z)");
        return RS_ERR_ARG;
    }

    /* The check that matters: refuse a grid finer than the geometry supports.
     *
     * Model C inverts with the radar wavelength, not the acoustic substitution,
     * so applying the acoustic figure to it would reject grids it can perfectly
     * well support -- and would do so in the direction that makes the
     * uncontested reference model look worse than it is. */
    const double lam = (params->model == RS_TOMO_MODEL_C)
        ? ((params->radar_wavelength > 0.0) ? params->radar_wavelength : 0.031)
        : rs_acoustic_wavelength(params->velocity, params->frequency,
                                 params->convention);
    /* A literal 2*pi*t multiplier scales the sampled Kz span by the same
     * factor. This derived resolution intentionally differs from the patent's
     * separately printed lambda*R/(2*A) whenever the experimental factor is
     * enabled -- one of the source's internal inconsistencies. */
    const double dT = rs_tomo_resolution(lam, params->slant_range, params->aperture)
                    / rs_tomo_eq22_scale(params);
    if (dT > 0.0 && params->depth_cell < dT) {
        rs_set_error("tomo: depth cell %g m is finer than the resolution %g m that "
                     "v=%g m/s, f=%g Hz, R=%g m, A=%g m support; a finer grid "
                     "interpolates, it does not resolve",
                     params->depth_cell, dT, params->velocity, params->frequency,
                     params->slant_range, params->aperture);
        return RS_ERR_RANGE;
    }

    return RS_OK;
}

/* Release everything a tomogram owns. */
void rs_tomo_free(rs_tomo_t *t)
{
    if (!t) return;
    free(t->profile);
    free(t->depth);
    free(t->quality);
    memset(t, 0, sizeof *t);
}

/* Greatest unambiguously representable depth for a given look count. */
double rs_tomo_max_depth(const rs_tomo_params_t *params, size_t n_looks)
{
    if (!params || n_looks == 0) return 0.0;
    /* Each model's own wavelength; see rs_tomo_params_check(). */
    const double lam = (params->model == RS_TOMO_MODEL_C)
        ? ((params->radar_wavelength > 0.0) ? params->radar_wavelength : 0.031)
        : rs_acoustic_wavelength(params->velocity, params->frequency,
                                 params->convention);
    const double kz_max = rs_tomo_wavenumber(params->aperture, lam,
                                             params->slant_range, params->incidence)
                        * rs_tomo_eq22_scale(params);
    if (!(kz_max > 0.0)) return 0.0;
    return 2.0 * M_PI * (double)n_looks / kz_max;
}

/* Assemble the data vector Y of Eq. 21 for one window.
 *
 * Under RS_TOMO_Y_SHIFTS the two coregistrator shift components become the real
 * and imaginary parts, converted from pixels to metres -- the paper's own
 * construction, see rs_tomo_y_source_t. Under RS_TOMO_Y_LOS the line-of-sight
 * displacement is used as a real series.
 *
 * The mean is removed either way: a constant offset across looks is a
 * registration bias, not a scatterer, and leaving it in deposits a spurious
 * feature at zero depth in every window. */
static void rs_tomo_build_y(const rs_microm_t *m, size_t w,
                            rs_tomo_y_source_t source, int remove_mean,
                            double complex *y)
{
    const size_t k = m->n_looks;
    const double az_sp = (m->az_spacing_m > 0.0) ? m->az_spacing_m : 1.0;
    const double rg_sp = (m->rg_spacing_m > 0.0) ? m->rg_spacing_m : 1.0;

    double complex mean = 0.0;
    for (size_t i = 0; i < k; i++) {
        const size_t idx = w * k + i;
        y[i] = (source == RS_TOMO_Y_LOS)
             ? (double complex)m->disp_los[idx]
             : m->disp_az[idx] * az_sp + I * (m->disp_rg[idx] * rg_sp);
        mean += y[i];
    }
    mean /= (double)k;
    if (remove_mean) for (size_t i = 0; i < k; i++) y[i] -= mean;
}

/* Model A via the Fourier route.
 *
 * Both primary sources state that the steering matrix is a DFT operator, so
 * this evaluates it as one: the data vector Y is transformed along the
 * sub-aperture index and the resulting frequency axis is relabelled to depth
 * through the Kz mapping.
 *
 * A Hann window is applied before the transform. Without it the depth profile
 * carries the sidelobes of a rectangular k-point transform, which appear as
 * evenly spaced "layers" beneath any bright scatterer -- an artefact that is
 * very easy to mistake for structure, and precisely the kind of artefact this
 * project exists to rule out.
 *
 * The depth grid is resampled from the transform's natural frequency grid by
 * linear interpolation, so that the caller gets the depth axis it asked for
 * rather than one dictated by the look count. */
static resonarsat_status_t rs_tomo_model_a_dft(const rs_microm_t *m,
                                               const rs_tomo_params_t *params,
                                               double lambda_ac,
                                               rs_tomo_t *out)
{
    const size_t k = m->n_looks;
    const size_t n_depth = out->n_depth;

    rs_fft_plan *plan = NULL;
    resonarsat_status_t st = rs_fft_plan_create(k, &plan);
    if (st != RS_OK) return st;

    float complex *buf = malloc(k * sizeof *buf);
    double complex *y = malloc(k * sizeof *y);
    double *win = malloc(k * sizeof *win);
    double *mag = malloc(k * sizeof *mag);
    if (!buf || !y || !win || !mag) {
        free(buf); free(y); free(win); free(mag);
        rs_fft_plan_destroy(plan);
        return RS_ERR_ALLOC;
    }

    for (size_t i = 0; i < k; i++) {
        win[i] = params->window
               ? 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)(k - 1)))
               : 1.0;
    }

    /* Map spectral index to depth. The along-orbit baseline of sub-aperture i
     * relative to the reference is (i/k) * A, so Kz grows linearly with index
     * and the transform's conjugate variable is proportional to depth. The
     * span the transform can represent is set by the largest Kz step. */
    const double b_max = params->aperture;
    const double kz_max = rs_tomo_wavenumber(b_max, lambda_ac, params->slant_range,
                                             params->incidence);
    /* Depth corresponding to one full turn across the aperture. */
    const double depth_span = (kz_max > 0.0) ? (2.0 * M_PI * (double)k / kz_max) : 0.0;
    const double depth_step = (k > 0 && depth_span > 0.0) ? depth_span / (double)k : 0.0;

    for (size_t w = 0; w < m->n_win; w++) {
        rs_tomo_build_y(m, w, params->y_source, params->remove_y_mean, y);

        for (size_t i = 0; i < k; i++) {
            buf[i] = (float)(creal(y[i]) * win[i]) + (float)(cimag(y[i]) * win[i]) * I;
        }

        if (rs_fft_forward(plan, buf) != RS_OK) continue;
        for (size_t i = 0; i < k; i++) mag[i] = (double)cabsf(buf[i]);

        /* Resample the transform's natural grid onto the requested depths. */
        double *prof = out->profile + w * n_depth;
        for (size_t d = 0; d < n_depth; d++) {
            const double z = out->depth[d];
            if (depth_step <= 0.0) { prof[d] = 0.0; continue; }
            const double fidx = z / depth_step;
            if (fidx < 0.0 || fidx >= (double)(k - 1)) { prof[d] = 0.0; continue; }
            const size_t i0 = (size_t)fidx;
            const double fr = fidx - (double)i0;
            prof[d] = mag[i0] * (1.0 - fr) + mag[i0 + 1] * fr;
        }
    }

    free(buf); free(y); free(win); free(mag);
    rs_fft_plan_destroy(plan);
    return RS_OK;
}

/* LU factorise an n-by-n matrix in place with partial pivoting. */
static int rs_lu_factor(double complex *M, size_t n, size_t *piv)
{
    for (size_t i = 0; i < n; i++) piv[i] = i;
    for (size_t c = 0; c < n; c++) {
        size_t best = c;
        double best_mag = cabs(M[c * n + c]);
        for (size_t r = c + 1; r < n; r++) {
            const double mg = cabs(M[r * n + c]);
            if (mg > best_mag) { best_mag = mg; best = r; }
        }
        if (best_mag < 1e-300) return 0;
        if (best != c) {
            for (size_t q = 0; q < n; q++) {
                double complex t = M[c * n + q];
                M[c * n + q] = M[best * n + q];
                M[best * n + q] = t;
            }
            size_t tp = piv[c]; piv[c] = piv[best]; piv[best] = tp;
        }
        for (size_t r = c + 1; r < n; r++) {
            const double complex f = M[r * n + c] / M[c * n + c];
            M[r * n + c] = f;
            for (size_t q = c + 1; q < n; q++) M[r * n + q] -= f * M[c * n + q];
        }
    }
    return 1;
}

/* Solve using a factorisation from rs_lu_factor(). 'b' is overwritten. */
static void rs_lu_solve(const double complex *M, size_t n, const size_t *piv,
                        const double complex *rhs, double complex *b)
{
    for (size_t i = 0; i < n; i++) b[i] = rhs[piv[i]];
    for (size_t i = 1; i < n; i++)
        for (size_t q = 0; q < i; q++) b[i] -= M[i * n + q] * b[q];
    for (size_t i = n; i-- > 0; ) {
        for (size_t q = i + 1; q < n; q++) b[i] -= M[i * n + q] * b[q];
        b[i] /= M[i * n + i];
    }
}

/* Model A via explicit least squares, kept as a cross-check on the DFT route.
 *
 * Builds the dimensionally conventional steering matrix
 * A[i][j] = exp(j*Kz_i*z_j) and solves the Tikhonov-regularised normal
 * equations (A^H A + mu I) h = A^H Y by Gaussian elimination with partial
 * pivoting.
 *
 * This is NOT rendered Eq. 22 literally. The PDF prints
 * exp(j*2*pi*Kz_i*t*z_j), while defining Kz_i with 4*pi already and never
 * defining 't' in this equation. Its prose says A is k-by-F, its displayed
 * matrix is F-by-k, and Eq. 23 only multiplies consistently under the former.
 * We use k-by-F and, by default, omit the dimensionally inconsistent 2*pi*t.
 * The experimental eq22_literal_t option applies that printed scale explicitly.
 *
 * Regularisation is not optional here. The along-orbit baseline spread is
 * small, so the steering matrix is severely ill-conditioned and an unregularised
 * solve amplifies noise into exactly the kind of sharp spurious peaks that would
 * be reported as discovered structure. The condition of the system is what the
 * open question about this geometry looks like numerically.
 *
 * Solved in double precision throughout; the matrices are k-by-F with k of order
 * tens, so the cost is irrelevant next to the tracking stage. */
static resonarsat_status_t rs_tomo_model_a_lstsq(const rs_microm_t *m,
                                                 const rs_tomo_params_t *params,
                                                 double lambda_ac,
                                                 rs_tomo_t *out)
{
    const size_t k = m->n_looks;
    const size_t F = out->n_depth;
    const double eq22_scale = rs_tomo_eq22_scale(params);

    double complex *A = malloc(k * F * sizeof *A);
    double complex *rhs = malloc(F * sizeof *rhs);
    double complex *y = malloc(k * sizeof *y);
    if (!A || !rhs || !y) {
        free(A); free(rhs); free(y);
        return RS_ERR_ALLOC;
    }

    /* Steering matrix: baseline of look i is (i/k)*A along the orbit. */
    for (size_t i = 0; i < k; i++) {
        const double b_perp = params->aperture * ((double)i / (double)k);
        const double kz = rs_tomo_wavenumber(b_perp, lambda_ac, params->slant_range,
                                             params->incidence) * eq22_scale;
        for (size_t j = 0; j < F; j++) {
            A[i * F + j] = cexp(I * kz * out->depth[j]);
        }
    }

    /* The pseudoinverse of Eq. 24 has two forms and the right one depends on
     * the shape of the steering matrix. Using the wrong one is not a numerical
     * detail -- it is a singular matrix.
     *
     *   k >= F, more sub-looks than depth cells: A has full column rank and
     *           A_dagger = (A^H A)^-1 A^H. The normal matrix is F-by-F.
     *
     *   F > k,  more depth cells than sub-looks: A^H A is F-by-F with rank at
     *           most k, so it is SINGULAR and no amount of care in the solve
     *           recovers it. The minimum-norm pseudoinverse is instead
     *           A_dagger = A^H (A A^H)^-1, a k-by-k solve.
     *
     * With regularisation zero, each branch computes the exact Moore-Penrose
     * pseudoinverse the patent states. A non-zero factor adds a Tikhonov ridge,
     * which stabilises an ill-conditioned steering matrix and is a departure
     * from Eq. 24 -- see rs_tomo_params_t.regularisation. */
    const int overdet = (k >= F);
    const size_t n = overdet ? F : k;

    double complex *M = malloc(n * n * sizeof *M);
    size_t *piv = malloc(n * sizeof *piv);
    double complex *b = malloc(n * sizeof *b);
    if (!M || !piv || !b) {
        free(A); free(rhs); free(y); free(M); free(piv); free(b);
        return RS_ERR_ALLOC;
    }

    double trace = 0.0;
    if (overdet) {
        for (size_t a = 0; a < F; a++) {
            for (size_t bb = 0; bb < F; bb++) {
                double complex sum = 0.0;
                for (size_t i = 0; i < k; i++) sum += conj(A[i * F + a]) * A[i * F + bb];
                M[a * F + bb] = sum;
            }
            trace += creal(M[a * F + a]);
        }
    } else {
        for (size_t a = 0; a < k; a++) {
            for (size_t bb = 0; bb < k; bb++) {
                double complex sum = 0.0;
                for (size_t j = 0; j < F; j++) sum += A[a * F + j] * conj(A[bb * F + j]);
                M[a * k + bb] = sum;
            }
            trace += creal(M[a * k + a]);
        }
    }
    const double mu = params->regularisation * (trace / (double)n);
    for (size_t a = 0; a < n; a++) M[a * n + a] += mu;

    if (!rs_lu_factor(M, n, piv)) {
        rs_set_error("tomo: the %s matrix is singular at regularisation %g; "
                     "raise --regularisation, or change the depth grid so that "
                     "the steering matrix is better conditioned",
                     overdet ? "normal" : "Gram", params->regularisation);
        free(A); free(rhs); free(y); free(M); free(piv); free(b);
        return RS_ERR_SINGULAR;
    }

    for (size_t w = 0; w < m->n_win; w++) {
        rs_tomo_build_y(m, w, params->y_source, params->remove_y_mean, y);
        double *prof = out->profile + w * F;

        if (overdet) {
            for (size_t j = 0; j < F; j++) {
                double complex sum = 0.0;
                for (size_t i = 0; i < k; i++) sum += conj(A[i * F + j]) * y[i];
                rhs[j] = sum;
            }
            rs_lu_solve(M, n, piv, rhs, b);
            for (size_t j = 0; j < F; j++) prof[j] = cabs(b[j]);
        } else {
            /* Solve (A A^H) w = Y, then h = A^H w. */
            rs_lu_solve(M, n, piv, y, b);
            for (size_t j = 0; j < F; j++) {
                double complex sum = 0.0;
                for (size_t i = 0; i < k; i++) sum += conj(A[i * F + j]) * b[i];
                prof[j] = cabs(sum);
            }
        }
    }

    free(A); free(rhs); free(y); free(M); free(piv); free(b);
    return RS_OK;
}

/* Model B: standing-wave back-projection.
 *
 * Each spectral peak is read as a standing wave of acoustic wavelength
 * lambda_i = v/(C*f_i), and energy is deposited at the depths where such a wave
 * has antinodes, d = n*lambda_i/2, with a Gaussian kernel one resolution cell
 * wide. This is a frankly heuristic model and is included because it fails in a
 * different way from Model A: agreement between two models with different
 * failure modes is weak evidence, but disagreement is informative.
 *
 * Note that this is the one place where a MEASURED vibration frequency enters a
 * depth mapping, which is a deliberate difference from Model A. It is what the
 * standing-wave interpretation means; it is also why this model is a
 * cross-check and not the primary. */
static resonarsat_status_t rs_tomo_model_b(const rs_spectrum_t *spec,
                                           const rs_tomo_params_t *params,
                                           rs_tomo_t *out)
{
    if (!spec) {
        rs_set_error("tomo: model B requires vibration spectra");
        return RS_ERR_ARG;
    }

    const double dT = out->dT > 0.0 ? out->dT : params->depth_cell;
    const double sigma = fmax(dT, params->depth_cell);

    for (size_t w = 0; w < out->n_win; w++) {
        const double *psd = spec->psd + w * spec->n_freq;
        double *prof = out->profile + w * out->n_depth;

        for (size_t k = 1; k < spec->n_freq; k++) {
            const double f_i = spec->freq[k];
            if (f_i <= 0.0 || psd[k] <= 0.0) continue;

            const double lam_i = rs_acoustic_wavelength(params->velocity, f_i,
                                                        params->convention);
            if (lam_i <= 0.0) continue;

            const double amp = sqrt(psd[k]);
            for (int n = 1; n * lam_i / 2.0 <= params->depth_max; n++) {
                const double d = (double)n * lam_i / 2.0;
                for (size_t j = 0; j < out->n_depth; j++) {
                    const double dz = (out->depth[j] - d) / sigma;
                    if (fabs(dz) > 3.0) continue;
                    prof[j] += amp * exp(-0.5 * dz * dz);
                }
            }
        }
    }
    return RS_OK;
}

/* Model D: depth from the measured vibration frequency, z = v / (2 f).
 *
 * Model A treats the sub-aperture index as a tomographic baseline, so its depth
 * axis is linear in the transform's bin index. This treats the same index as
 * time, transforms to frequency, and reads depth from a resonance condition: a
 * layer of thickness d rings at f = v/(2d), the mass-and-spring picture the
 * presentation material sets out.
 *
 * The two differ in kind, not degree. Model A's depth grows with bin index;
 * here it falls as 1/f, so the shallowest representable depth is set by the
 * HIGHEST frequency the sub-look series supports and the deepest by the lowest.
 * That inverts which end of the spectrum is well sampled, and it is what lets
 * this reach the hundreds of metres to kilometres the source reports.
 *
 * Only the velocity is assumed here. The frequency is read from the data, which
 * is the substantive difference from Model A and the reason this deserved
 * implementing rather than arguing about: it removes the second free constant.
 *
 * Energy is accumulated rather than interpolated. Each spectral bin maps to one
 * depth, the mapping is non-uniform, and several bins can fall in one depth
 * cell at the shallow end while deep cells may receive none. Interpolating
 * across that would invent a smoothness the measurement does not have. */
static resonarsat_status_t rs_tomo_model_d(const rs_microm_t *m,
                                           const rs_spectrum_t *spec,
                                           const rs_tomo_params_t *params,
                                           rs_tomo_t *out)
{
    if (!spec || !spec->psd || spec->n_freq < 2) {
        rs_set_error("tomo: model D needs a spectrum; pass one to rs_tomo_focus()");
        return RS_ERR_ARG;
    }

    const double v = params->velocity;
    const size_t n_depth = out->n_depth;

    /* The band this model can actually represent, so a caller asking for depths
     * outside it is told rather than handed empty cells. */
    const double f_hi = spec->freq[spec->n_freq - 1];
    const double f_lo = spec->freq[1];
    if (f_hi > 0.0 && f_lo > 0.0) {
        const double z_shallow = v / (2.0 * f_hi);
        const double z_deep    = v / (2.0 * f_lo);
        if (out->depth[n_depth - 1] < z_shallow || out->depth[0] > z_deep) {
            rs_set_error("tomo: model D represents %g to %g m at v=%g with a "
                         "%g to %g Hz spectrum; the requested %g to %g m lies "
                         "outside it entirely",
                         z_shallow, z_deep, v, f_lo, f_hi,
                         out->depth[0], out->depth[n_depth - 1]);
            return RS_ERR_RANGE;
        }
    }

    for (size_t w = 0; w < m->n_win && w < spec->n_win; w++) {
        double *prof = out->profile + w * n_depth;
        const double *psd = spec->psd + w * spec->n_freq;
        size_t *hits = calloc(n_depth, sizeof *hits);
        if (!hits) return RS_ERR_ALLOC;

        /* Skip bin 0: it is the DC term, and z = v/(2*0) is unbounded. */
        for (size_t k = 1; k < spec->n_freq; k++) {
            const double f = spec->freq[k];
            if (!(f > 0.0)) continue;
            const double z = v / (2.0 * f);
            if (z < out->depth[0] || z > out->depth[n_depth - 1]) continue;
            const size_t j = (size_t)((z - out->depth[0]) / params->depth_cell + 0.5);
            if (j >= n_depth) continue;
            prof[j] += psd[k];
            hits[j]++;
        }
        /* Average within a cell so that cells collecting many bins are not
         * favoured purely for being wide in frequency. */
        for (size_t j = 0; j < n_depth; j++) {
            if (hits[j] > 1) prof[j] /= (double)hits[j];
        }
        free(hits);
    }
    return RS_OK;
}

/* Model C: classic multi-baseline tomography over genuine perpendicular
 * baselines.
 *
 * Beamforming rather than inversion: for each depth, coherently sum the
 * observations with the phase the geometry predicts for a scatterer there. This
 * is the standard, uncontested formulation, and its resolution is bounded by
 * real baseline diversity. It is the reference against which Model A's apparent
 * depth resolution should be judged: if Model A reports finer structure than
 * Model C can support on the same scene, that gap is the disputed claim,
 * quantified.
 *
 * Uses the electromagnetic wavelength, not the acoustic one -- this model makes
 * no acoustic substitution. */
static resonarsat_status_t rs_tomo_model_c(const rs_microm_t *m,
                                           const rs_tomo_params_t *params,
                                           const double *baselines,
                                           double lambda_em,
                                           rs_tomo_t *out)
{
    if (!baselines) {
        rs_set_error("tomo: model C requires an array of %zu perpendicular "
                     "baselines in metres", m->n_looks);
        return RS_ERR_ARG;
    }

    const size_t k = m->n_looks;
    double *kz = malloc(k * sizeof *kz);
    if (!kz) return RS_ERR_ALLOC;

    for (size_t i = 0; i < k; i++) {
        kz[i] = rs_tomo_wavenumber(baselines[i], lambda_em, params->slant_range,
                                   params->incidence);
    }

    double complex *y = malloc(k * sizeof *y);
    if (!y) { free(kz); return RS_ERR_ALLOC; }

    for (size_t w = 0; w < m->n_win; w++) {
        rs_tomo_build_y(m, w, params->y_source, params->remove_y_mean, y);

        double *prof = out->profile + w * out->n_depth;
        for (size_t j = 0; j < out->n_depth; j++) {
            double complex acc = 0.0;
            for (size_t i = 0; i < k; i++) {
                acc += y[i] * cexp(-I * kz[i] * out->depth[j]);
            }
            prof[j] = cabs(acc) / (double)k;
        }
    }

    free(y);
    free(kz);
    return RS_OK;
}

/* Focus vibration observations into depth profiles. */
resonarsat_status_t rs_tomo_focus(const rs_microm_t *m,
                                  const rs_spectrum_t *spec,
                                  const rs_tomo_params_t *params,
                                  const double *baselines,
                                  rs_tomo_t *out)
{
    if (!m || !params || !out || !m->disp_los) return RS_ERR_ARG;

    /* Zero the output before anything can fail. A caller that checks the status
     * and then frees the struct on the error path -- which is the ordinary
     * pattern, and what the test suite does -- would otherwise be freeing
     * whatever happened to be on its stack. */
    memset(out, 0, sizeof *out);

    resonarsat_status_t st = rs_tomo_params_check(params);
    if (st != RS_OK) return st;

    const double lambda_ac = rs_acoustic_wavelength(params->velocity, params->frequency,
                                                    params->convention);

    /* Refuse a depth extent the geometry cannot represent. Beyond the
     * unambiguous range the transform folds, and the previous behaviour --
     * zero-filling everything past the limit -- presented a profile that was
     * mostly fabricated emptiness, with a single surviving cell that looked
     * like a confident detection. Models A and C both sample Kz; Model B
     * deposits energy directly at computed depths and is not periodic, so the
     * limit does not apply to it. */
    /* Model D does not sample Kz either, and applying this limit to it would
     * repeat the error of judging Model C by Model A's wavelength (bug 18). Its
     * depth axis comes from the measured frequency band, so its range is set by
     * the dwell and the sub-look count -- z from v/(2*f_max) to v/(2*f_min) --
     * and is validated inside rs_tomo_model_d() where the spectrum is in hand.
     * Refusing a 4 km grid here because a baseline formula says 54 m would
     * refuse the very thing the model exists to represent. */
    if (params->model != RS_TOMO_MODEL_B && params->model != RS_TOMO_MODEL_D) {
        const double z_max = rs_tomo_max_depth(params, m->n_looks);
        if (z_max > 0.0 && params->depth_max > z_max) {
            rs_set_error("tomo: depth extent %g m exceeds the %g m this geometry "
                         "represents without ambiguity at %zu sub-apertures "
                         "(v=%g, f=%g, A=%g, R=%g). Beyond it the profile folds. "
                         "Use --depth %g or fewer, or raise --n to about %zu",
                         params->depth_max, z_max, m->n_looks,
                         params->velocity, params->frequency,
                         params->aperture, params->slant_range,
                         z_max,
                         (size_t)ceil((double)m->n_looks * params->depth_max / z_max));
            return RS_ERR_RANGE;
        }
    }

    const size_t n_depth = (size_t)(params->depth_max / params->depth_cell) + 1;

    out->profile = calloc(m->n_win * n_depth, sizeof *out->profile);
    out->depth   = calloc(n_depth, sizeof *out->depth);
    out->quality = calloc(m->n_win, sizeof *out->quality);
    if (!out->profile || !out->depth || !out->quality) {
        rs_tomo_free(out);
        return RS_ERR_ALLOC;
    }

    out->n_win = m->n_win;
    out->n_win_az = m->n_win_az;
    out->n_win_rg = m->n_win_rg;
    out->n_depth = n_depth;
    out->lambda_ac = lambda_ac;
    /* Resolution is quoted in the wavelength each model actually inverts with.
     * Model C uses the radar wavelength, not the acoustic substitution, and
     * quoting the acoustic figure for it would make any comparison between the
     * two models meaningless in exactly the direction that flatters Model A. */
    out->dT = rs_tomo_resolution(
        (params->model == RS_TOMO_MODEL_C)
          ? ((params->radar_wavelength > 0.0) ? params->radar_wavelength : 0.031)
          : lambda_ac,
        params->slant_range, params->aperture);
    if (params->model == RS_TOMO_MODEL_A)
        out->dT /= rs_tomo_eq22_scale(params);
    out->n_looks = m->n_looks;
    out->z_unambiguous = rs_tomo_max_depth(params, m->n_looks);
    out->params = *params;

    for (size_t j = 0; j < n_depth; j++) out->depth[j] = (double)j * params->depth_cell;
    memcpy(out->quality, m->quality, m->n_win * sizeof *out->quality);

    switch (params->model) {
    case RS_TOMO_MODEL_A:
        st = (params->solver == RS_TOMO_SOLVER_LSTSQ)
           ? rs_tomo_model_a_lstsq(m, params, lambda_ac, out)
           : rs_tomo_model_a_dft(m, params, lambda_ac, out);
        break;
    case RS_TOMO_MODEL_B:
        st = rs_tomo_model_b(spec, params, out);
        break;
    case RS_TOMO_MODEL_C:
        /* The electromagnetic wavelength: Model C makes no acoustic
         * substitution. X-band is assumed when nothing better is known. */
        st = rs_tomo_model_c(m, params, baselines,
                             (params->radar_wavelength > 0.0)
                               ? params->radar_wavelength : 0.031,
                             out);
        break;
    case RS_TOMO_MODEL_D:
        st = rs_tomo_model_d(m, spec, params, out);
        break;
    default:
        rs_set_error("tomo: unknown model %d", (int)params->model);
        st = RS_ERR_ARG;
        break;
    }

    if (st != RS_OK) rs_tomo_free(out);
    return st;
}

/* Emit the parameters and derived constants that must travel with a tomogram. */
void rs_tomo_write_metadata(const rs_tomo_t *t, void *stream)
{
    if (!t || !stream) return;
    FILE *f = (FILE *)stream;

    /* Every model must be named here. A cube whose sidecar says "unknown" is
     * unattributable, and Model D was omitted from this table for as long as it
     * existed -- so the comparison in MODEL-A-VS-D.md was run on products that
     * did not record which model produced them. The bound below is derived from
     * the table rather than written as a literal so that adding a model and
     * forgetting to name it cannot silently reintroduce that. */
    static const char *model_name[] = { "A (steering-matrix, along-orbit)",
                                        "B (standing-wave back-projection)",
                                        "C (multi-baseline, uncontested)",
                                        "D (resonance mapping, z = v/2f)" };
    const int mi = (int)t->params.model;
    const int n_model = (int)(sizeof model_name / sizeof model_name[0]);

    fprintf(f, "model                 %s\n",
            (mi >= 0 && mi < n_model) ? model_name[mi] : "unknown");

    /* Which arithmetic produced the measurement, printed high in the sidecar and
     * tagged so it reads at a glance.
     *
     * The two modes can report different tracked shifts, so two cubes of the same
     * scene may differ for this reason and no other. A reader who does not know
     * which is which will reach for a physical explanation of the difference. */
    fprintf(f, "arithmetic_mode       %s\n",
            t->params.no_optimize
              ? "[UNOPTIMIZED] exhaustive correlation peak search, serial execution"
              : "optimised (local peak refinement, threaded)");
    fprintf(f, "measurement_chain     %s%s\n",
            t->params.no_optimize ? "[UNOPTIMIZED] " : "",
            t->params.provenance[0] ? t->params.provenance : "UNRECORDED");
    fprintf(f, "solver                %s\n",
            t->params.solver == RS_TOMO_SOLVER_LSTSQ ? "least squares" : "DFT");
    fprintf(f, "assumed_velocity_ms   %g\n", t->params.velocity);
    fprintf(f, "assumed_frequency_hz  %g\n", t->params.frequency);
    fprintf(f, "y_vector_source       %s\n",
            t->params.y_source == RS_TOMO_Y_SHIFTS
              ? "complex tracked shifts (paper Eq. 20-21)"
              : "real line-of-sight displacement");
    fprintf(f, "wavelength_convention %s\n",
            t->params.convention == RS_WAVELEN_PAPER ? "paper (v/2f)" : "patent (v/f)");
    fprintf(f, "acoustic_wavelength_m %g\n", t->lambda_ac);
    fprintf(f, "slant_range_m         %g\n", t->params.slant_range);
    fprintf(f, "aperture_m            %g\n", t->params.aperture);
    fprintf(f, "incidence_deg         %g\n", t->params.incidence * 180.0 / M_PI);
    fprintf(f, "depth_resolution_m    %g\n", t->dT);
    fprintf(f, "depth_unambiguous_m   %g\n", t->z_unambiguous);
    fprintf(f, "n_sub_apertures       %zu\n", t->n_looks);
    fprintf(f, "depth_cell_m          %g\n", t->params.depth_cell);
    fprintf(f, "depth_max_m           %g\n", t->params.depth_max);
    fprintf(f, "n_windows             %zu\n", t->n_win);
    fprintf(f, "n_depth_cells         %zu\n", t->n_depth);

    /* The conditioning steps, printed because a reader cannot otherwise tell a
     * product that follows the published equations from one that does not. The
     * three below are additions this project makes and the sources do not
     * describe; each is defensible and none is in the papers. */
    fprintf(f, "y_mean_removed        %s\n", t->params.remove_y_mean ? "yes" : "no");
    fprintf(f, "depth_taper           %s\n", t->params.window ? "hann" : "none");
    fprintf(f, "subaperture_taper     %s\n",
            t->params.subap_window ? "raised-cosine" : "none");
    fprintf(f, "coherence_min         %g\n", t->params.coherence_min);
    fprintf(f, "tracking_reference    %s\n",
            t->params.pair_reference ? "master-slave pair (B_shift apart)"
                                     : "not the master-slave pair");
    fprintf(f, "regularisation        %g%s\n", t->params.regularisation,
            (t->params.solver == RS_TOMO_SOLVER_LSTSQ &&
             t->params.regularisation > 0.0) ? " (Tikhonov, not a pseudoinverse)" : "");
    fprintf(f, "eq22_steering         %s\n",
            t->params.eq22_literal_t > 0.0
              ? "EXPERIMENTAL literal exp(j*2*pi*Kz*t*z)"
              : "conventional exp(j*Kz*z)");
    fprintf(f, "eq22_literal_t        %g\n", t->params.eq22_literal_t);
    fprintf(f, "eq22_kz_scale         %g\n",
            rs_tomo_eq22_scale(&t->params));

    if (t->params.model == RS_TOMO_MODEL_A) {
        /* 'exact' below identifies the unconditioned choices described by the
         * source. It deliberately does NOT certify Eqs. 21-24 "as written":
         * rendered Eq. 22 contains an undefined, dimensionally inconsistent
         * 2*pi*t factor and the Eq. 22-23 matrix dimensions conflict. */
        const int exact = t->params.patent_exact &&
                          t->params.y_source == RS_TOMO_Y_SHIFTS &&
                          t->params.convention == RS_WAVELEN_PATENT &&
                          t->params.solver == RS_TOMO_SOLVER_LSTSQ &&
                          !t->params.remove_y_mean &&
                          !t->params.window &&
                          !t->params.subap_window &&
                          t->params.pair_reference &&
                          t->params.coherence_min == 0.0 &&
                          !(t->params.solver == RS_TOMO_SOLVER_LSTSQ &&
                            t->params.regularisation > 0.0);
        if (exact && t->params.eq22_literal_t > 0.0) {
            fprintf(f,
                "\nEXPERIMENTAL LITERAL EQ. 22: the steering exponent includes\n"
                "the rendered 2*pi*t factor, with caller-supplied t=%g.\n"
                "The source does not define t here and its physical dimensions\n"
                "remain unresolved; this mode is not conventional TomoSAR.\n",
                t->params.eq22_literal_t);
        } else if (exact) {
            fprintf(f,
                "\nPATENT-CHAIN INTERPRETATION: raw Eq. 21 data, the source's\n"
                "unconditioned Figure 0.5 choices, and Eq. 24's pseudoinverse.\n"
                "NOT LITERAL EQ. 22: the rendered source prints an undefined\n"
                "2*pi*t factor and mutually inconsistent matrix dimensions;\n"
                "this code uses the conventional exp(j*Kz*z), k-by-F form.\n");
        } else {
            fprintf(f,
                "\nNOT THE UNCONDITIONED PATENT CHAIN. This run applies\n"
                "conditioning the sources do not describe -- see the three fields\n"
                "above. Each removes a measured artefact and each changes the\n"
                "profile. --patent-exact selects the patent-chain interpretation.\n"
                "Use --eq22-literal-t only to test the rendered Eq. 22 factor;\n"
                "its t is undefined by the source and is not conventional.\n");
        }
        fprintf(f,
            "\nCAVEAT: in Eqs. 21-24 the tomographic baseline is formed from\n"
            "along-orbit sub-aperture phase-centre separations. An along-track\n"
            "separation is not an elevation baseline. Whether such a geometry\n"
            "carries depth information is the method's principal open question.\n"
            "The depth axis above is scaled by the ASSUMED velocity and frequency\n"
            "recorded here, neither of which is measured by this software. Do not\n"
            "present these depths as measurements without the null test and the\n"
            "parameter-sensitivity sweep alongside.\n");
    }
}

/* Locate the strongest feature in a set of depth profiles.
 *
 * Averages the profiles across windows before searching, so the reported peak
 * describes the scene rather than whichever single window happened to be
 * noisiest. The zero-depth cell is excluded: every profile has energy at zero
 * by construction and it is not a feature.
 *
 * The peak position is refined below one cell by fitting a parabola through the
 * strongest cell and its two neighbours. Without that refinement a sweep across
 * parameters reports the same quantised depth for every combination -- the
 * feature moves, but by less than a cell each time -- which looks exactly like
 * a depth that is pinned by the data when in fact it is pinned by the grid.
 * That is a false negative on the one question the sweep exists to answer. */
static void rs_tomo_peak(const rs_tomo_t *t, double *depth_out, double *value_out)
{
    *depth_out = 0.0;
    *value_out = 0.0;
    if (t->n_depth < 2 || t->n_win == 0) return;

    double *mean = malloc(t->n_depth * sizeof *mean);
    if (!mean) return;

    for (size_t j = 0; j < t->n_depth; j++) {
        double acc = 0.0;
        for (size_t w = 0; w < t->n_win; w++) acc += t->profile[w * t->n_depth + j];
        mean[j] = acc / (double)t->n_win;
    }

    size_t best = 1;
    for (size_t j = 2; j < t->n_depth; j++) if (mean[j] > mean[best]) best = j;

    double offset = 0.0;
    if (best > 0 && best + 1 < t->n_depth) {
        const double y0 = mean[best - 1], y1 = mean[best], y2 = mean[best + 1];
        const double denom = y0 - 2.0 * y1 + y2;
        if (denom != 0.0) {
            offset = 0.5 * (y0 - y2) / denom;
            if (offset < -1.0) offset = -1.0;
            if (offset >  1.0) offset =  1.0;
        }
    }

    *depth_out = t->depth[best] + offset * t->params.depth_cell;
    *value_out = mean[best];
    free(mean);
}

/* Sweep the two assumed constants and record where the strongest feature lands. */
resonarsat_status_t rs_tomo_sweep(const rs_microm_t *m,
                                  const rs_spectrum_t *spec,
                                  const rs_tomo_params_t *params,
                                  double scale_min, double scale_max, size_t n_scale,
                                  rs_tomo_sweep_row_t *rows, size_t *n_rows)
{
    if (!m || !params || !rows || !n_rows) return RS_ERR_ARG;
    if (n_scale < 2 || !(scale_min > 0.0) || !(scale_max > scale_min)) {
        rs_set_error("tomo sweep: need n_scale >= 2 and 0 < scale_min < scale_max");
        return RS_ERR_ARG;
    }

    *n_rows = 0;
    const double log_lo = log(scale_min);
    const double log_hi = log(scale_max);

    for (size_t iv = 0; iv < n_scale; iv++) {
        const double fv = exp(log_lo + (log_hi - log_lo) * (double)iv / (double)(n_scale - 1));
        for (size_t jf = 0; jf < n_scale; jf++) {
            const double ff = exp(log_lo + (log_hi - log_lo) * (double)jf / (double)(n_scale - 1));

            rs_tomo_params_t p = *params;
            p.velocity  = params->velocity * fv;
            p.frequency = params->frequency * ff;

            /* Skip combinations the depth grid cannot support rather than
             * quietly accepting a grid finer than the resolution. */
            if (rs_tomo_params_check(&p) != RS_OK) continue;

            rs_tomo_t t;
            if (rs_tomo_focus(m, spec, &p, NULL, &t) != RS_OK) continue;

            rs_tomo_sweep_row_t *r = &rows[(*n_rows)++];
            r->velocity  = p.velocity;
            r->frequency = p.frequency;
            r->lambda_ac = t.lambda_ac;
            r->dT        = t.dT;
            rs_tomo_peak(&t, &r->peak_depth, &r->peak_value);

            rs_tomo_free(&t);
        }
    }

    if (*n_rows == 0) {
        rs_set_error("tomo sweep: no parameter combination produced a usable grid; "
                     "widen --cell or narrow the sweep range");
        return RS_ERR_RANGE;
    }
    return RS_OK;
}

/* Group sweep rows by acoustic wavelength and average within each group. */
resonarsat_status_t rs_tomo_sweep_merge(const rs_tomo_sweep_row_t *rows, size_t n_rows,
                                        rs_tomo_sweep_row_t *out, double *sd_out,
                                        size_t *n_out)
{
    if (!rows || !out || !n_out || n_rows == 0) return RS_ERR_ARG;

    size_t n_group = 0;

    for (size_t i = 0; i < n_rows; i++) {
        if (rows[i].peak_depth <= 0.0 || rows[i].lambda_ac <= 0.0) continue;

        /* Already collected under this wavelength? */
        int seen = 0;
        for (size_t g = 0; g < n_group; g++) {
            if (fabs(out[g].lambda_ac - rows[i].lambda_ac) <=
                1e-9 * rows[i].lambda_ac) { seen = 1; break; }
        }
        if (seen) continue;

        /* Accumulate every row sharing this wavelength. */
        double sum = 0.0, sumsq = 0.0, sum_val = 0.0;
        size_t n = 0;
        for (size_t j = 0; j < n_rows; j++) {
            if (rows[j].peak_depth <= 0.0 || rows[j].lambda_ac <= 0.0) continue;
            if (fabs(rows[j].lambda_ac - rows[i].lambda_ac) >
                1e-9 * rows[i].lambda_ac) continue;
            sum += rows[j].peak_depth;
            sumsq += rows[j].peak_depth * rows[j].peak_depth;
            sum_val += rows[j].peak_value;
            n++;
        }
        if (n == 0) continue;

        out[n_group] = rows[i];
        out[n_group].peak_depth = sum / (double)n;
        out[n_group].peak_value = sum_val / (double)n;

        if (sd_out) {
            const double mean = sum / (double)n;
            const double var = sumsq / (double)n - mean * mean;
            sd_out[n_group] = (var > 0.0) ? sqrt(var) : 0.0;
        }
        n_group++;
    }

    if (n_group == 0) {
        rs_set_error("tomo sweep merge: no usable rows");
        return RS_ERR_ARG;
    }
    *n_out = n_group;
    return RS_OK;
}

/* Regress log(peak depth) on log(acoustic wavelength). */
resonarsat_status_t rs_tomo_sweep_summary(const rs_tomo_sweep_row_t *rows, size_t n_rows,
                                          double *slope, double *correlation)
{
    if (!rows || !slope || !correlation) return RS_ERR_ARG;

    /* Only rows with a positive peak depth can be logged, and only DISTINCT
     * acoustic wavelengths carry information.
     *
     * The sweep varies velocity and frequency independently, but they enter the
     * depth mapping solely through lambda_ac = v/(C*f), so a 5x5 grid contains
     * just 9 distinct wavelengths and many exact duplicates. Regressing over the
     * duplicates does not change the fitted slope but it inflates the apparent
     * sample size and makes the correlation look far better determined than the
     * evidence supports -- which is precisely the wrong way for this particular
     * statistic to mislead, since it is the one deciding whether the depth axis
     * carries information. */
    size_t n = 0;
    double sx = 0.0, sy = 0.0;
    for (size_t i = 0; i < n_rows; i++) {
        if (rows[i].peak_depth <= 0.0 || rows[i].lambda_ac <= 0.0) continue;

        int dup = 0;
        for (size_t j = 0; j < i; j++) {
            if (rows[j].peak_depth <= 0.0 || rows[j].lambda_ac <= 0.0) continue;
            if (fabs(rows[j].lambda_ac - rows[i].lambda_ac) <=
                1e-9 * rows[i].lambda_ac) { dup = 1; break; }
        }
        if (dup) continue;

        sx += log(rows[i].lambda_ac);
        sy += log(rows[i].peak_depth);
        n++;
    }
    if (n < 3) {
        rs_set_error("tomo sweep: only %zu distinct acoustic wavelengths, "
                     "need at least 3", n);
        return RS_ERR_ARG;
    }

    const double mx = sx / (double)n, my = sy / (double)n;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (size_t i = 0; i < n_rows; i++) {
        if (rows[i].peak_depth <= 0.0 || rows[i].lambda_ac <= 0.0) continue;
        int dup = 0;
        for (size_t j = 0; j < i; j++) {
            if (rows[j].peak_depth <= 0.0 || rows[j].lambda_ac <= 0.0) continue;
            if (fabs(rows[j].lambda_ac - rows[i].lambda_ac) <=
                1e-9 * rows[i].lambda_ac) { dup = 1; break; }
        }
        if (dup) continue;
        const double dx = log(rows[i].lambda_ac) - mx;
        const double dy = log(rows[i].peak_depth) - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }

    *slope = (sxx > 0.0) ? sxy / sxx : 0.0;
    *correlation = (sxx > 0.0 && syy > 0.0) ? sxy / sqrt(sxx * syy) : 0.0;
    return RS_OK;
}


/* Regression restricted to reproducible wavelengths. */
resonarsat_status_t rs_tomo_sweep_summary_reliable(const rs_tomo_sweep_row_t *rows,
                                                   const double *sd, size_t n_rows,
                                                   double max_rel_sd,
                                                   double *slope, double *correlation,
                                                   size_t *n_used)
{
    if (!rows || !sd || !slope || !correlation) return RS_ERR_ARG;
    if (!(max_rel_sd >= 0.0)) return RS_ERR_ARG;

    size_t n = 0;
    double sx = 0.0, sy = 0.0;
    for (size_t i = 0; i < n_rows; i++) {
        if (rows[i].peak_depth <= 0.0 || rows[i].lambda_ac <= 0.0) continue;
        if (sd[i] > max_rel_sd * rows[i].peak_depth) continue;
        sx += log(rows[i].lambda_ac);
        sy += log(rows[i].peak_depth);
        n++;
    }

    if (n_used) *n_used = n;
    if (n < 3) {
        rs_set_error("tomo sweep: only %zu of %zu wavelengths reproduce within "
                     "%.0f%%; run more realisations or use a finer depth grid",
                     n, n_rows, max_rel_sd * 100.0);
        return RS_ERR_ARG;
    }

    const double mx = sx / (double)n, my = sy / (double)n;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (size_t i = 0; i < n_rows; i++) {
        if (rows[i].peak_depth <= 0.0 || rows[i].lambda_ac <= 0.0) continue;
        if (sd[i] > max_rel_sd * rows[i].peak_depth) continue;
        const double dx = log(rows[i].lambda_ac) - mx;
        const double dy = log(rows[i].peak_depth) - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }

    *slope = (sxx > 0.0) ? sxy / sxx : 0.0;
    *correlation = (sxx > 0.0 && syy > 0.0) ? sxy / sqrt(sxx * syy) : 0.0;
    return RS_OK;
}
