/* Split-band Phase Linking shift estimation.
 * See include/resonarsat/phaselink.h for the method and its provenance. */

#include "resonarsat/phaselink.h"
#include "resonarsat/fft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define RS_PL_MAX_ITER 32
#define RS_PL_TOL      1e-9

/* Phase Linking over a stack's coherence matrix. */
resonarsat_status_t rs_phase_link(const float complex *sig,
                                  size_t n_sig, size_t n_samp,
                                  double *phase_out)
{
    if (!sig || !phase_out || n_sig == 0 || n_samp == 0) return RS_ERR_ARG;

    double complex *gamma = calloc(n_sig * n_sig, sizeof *gamma);
    double *energy = calloc(n_sig, sizeof *energy);
    double complex *x = calloc(n_sig, sizeof *x);
    double complex *xn = calloc(n_sig, sizeof *xn);
    if (!gamma || !energy || !x || !xn) {
        free(gamma); free(energy); free(x); free(xn);
        return RS_ERR_ALLOC;
    }

    for (size_t i = 0; i < n_sig; i++) {
        const float complex *a = sig + i * n_samp;
        double e = 0.0;
        for (size_t p = 0; p < n_samp; p++) {
            const double m = (double)cabsf(a[p]);
            e += m * m;
        }
        energy[i] = e;
    }

    /* Sample coherence over every pair. Using all N^2 interferograms rather
     * than the N-1 against a single reference is the reason this estimator
     * beats pairwise correlation. */
    for (size_t i = 0; i < n_sig; i++) {
        for (size_t j = i; j < n_sig; j++) {
            const float complex *a = sig + i * n_samp;
            const float complex *b = sig + j * n_samp;
            double complex acc = 0.0;
            for (size_t p = 0; p < n_samp; p++) {
                acc += (double complex)a[p] * conj((double complex)b[p]);
            }
            const double norm = sqrt(energy[i] * energy[j]);
            const double complex g = (norm > 0.0) ? acc / norm : 0.0;
            gamma[i * n_sig + j] = g;
            gamma[j * n_sig + i] = conj(g);
        }
    }

    for (size_t i = 0; i < n_sig; i++) x[i] = 1.0;

    /* Fixed-point iteration toward the maximum-likelihood phase vector. Each
     * component is re-estimated from every other, then renormalised to unit
     * modulus, which is what confines the solution to phases only. */
    for (int iter = 0; iter < RS_PL_MAX_ITER; iter++) {
        double delta = 0.0;
        for (size_t i = 0; i < n_sig; i++) {
            double complex acc = 0.0;
            for (size_t j = 0; j < n_sig; j++) {
                if (i == j) continue;
                acc += gamma[i * n_sig + j] * x[j];
            }
            const double mag = cabs(acc);
            xn[i] = (mag > 0.0) ? acc / mag : x[i];
            delta += cabs(xn[i] - x[i]);
        }
        memcpy(x, xn, n_sig * sizeof *x);
        if (delta < RS_PL_TOL) break;
    }

    /* Reference to signal 0 so the returned vector is the relative phases the
     * problem actually determines. */
    const double ref = carg(x[0]);
    for (size_t i = 0; i < n_sig; i++) {
        double p = carg(x[i]) - ref;
        while (p >  M_PI) p -= 2.0 * M_PI;
        while (p < -M_PI) p += 2.0 * M_PI;
        phase_out[i] = (energy[i] > 0.0) ? p : 0.0;
    }

    free(gamma); free(energy); free(x); free(xn);
    return RS_OK;
}

/* Mean azimuth power spectrum of a patch stack, fftshifted so index 0 is the
 * most negative frequency. Caller frees. */
static double *rs_pl_mean_spectrum(const float complex *patch, size_t n_look,
                                   size_t n_az, size_t n_rg)
{
    double *psd = calloc(n_az, sizeof *psd);
    float complex *col = malloc(n_az * sizeof *col);
    rs_fft_plan *plan = NULL;

    if (!psd || !col || rs_fft_plan_create(n_az, &plan) != RS_OK) {
        free(psd); free(col); rs_fft_plan_destroy(plan);
        return NULL;
    }

    for (size_t k = 0; k < n_look; k++) {
        const float complex *p = patch + k * n_az * n_rg;
        for (size_t r = 0; r < n_rg; r++) {
            for (size_t a = 0; a < n_az; a++) col[a] = p[a * n_rg + r];
            if (rs_fft_forward(plan, col) != RS_OK) continue;
            rs_fft_shift(col, n_az);
            for (size_t f = 0; f < n_az; f++) {
                const double m = (double)cabsf(col[f]);
                psd[f] += m * m;
            }
        }
    }

    rs_fft_plan_destroy(plan);
    free(col);
    return psd;
}

/* Locate the band holding RS_PL_ENERGY_FRAC of the spectrum's energy, growing
 * outward from the peak.
 *
 * This exists because the source assumes signals scaled to unit bandwidth while
 * these sub-looks are oversampled several times over. Taking thirds of the
 * SAMPLED band would place the outer thirds in empty spectrum and return pure
 * noise; taking thirds of the OCCUPIED band is what the method intends. */
static void rs_pl_occupied_band(const double *psd, size_t n, size_t *lo, size_t *hi)
{
    double total = 0.0;
    size_t peak = 0;
    for (size_t i = 0; i < n; i++) {
        total += psd[i];
        if (psd[i] > psd[peak]) peak = i;
    }

    if (total <= 0.0) { *lo = 0; *hi = n - 1; return; }

    double acc = psd[peak];
    size_t l = peak, h = peak;
    while (acc < RS_PL_ENERGY_FRAC * total && (l > 0 || h + 1 < n)) {
        const double left  = (l > 0)     ? psd[l - 1] : -1.0;
        const double right = (h + 1 < n) ? psd[h + 1] : -1.0;
        if (right >= left) { h++; acc += psd[h]; }
        else               { l--; acc += psd[l]; }
    }
    *lo = l;
    *hi = h;
}

/* Band-filter one patch stack along azimuth into 'out', keeping spectral indices
 * [f0, f1) of the fftshifted spectrum and zeroing the rest. */
static resonarsat_status_t rs_pl_bandpass(const float complex *patch, size_t n_look,
                                          size_t n_az, size_t n_rg,
                                          size_t f0, size_t f1,
                                          float complex *out)
{
    float complex *col = malloc(n_az * sizeof *col);
    rs_fft_plan *plan = NULL;
    resonarsat_status_t st = rs_fft_plan_create(n_az, &plan);
    if (!col || st != RS_OK) {
        free(col); rs_fft_plan_destroy(plan);
        return (st == RS_OK) ? RS_ERR_ALLOC : st;
    }

    for (size_t k = 0; k < n_look; k++) {
        const float complex *p = patch + k * n_az * n_rg;
        float complex *o = out + k * n_az * n_rg;
        for (size_t r = 0; r < n_rg; r++) {
            for (size_t a = 0; a < n_az; a++) col[a] = p[a * n_rg + r];
            if ((st = rs_fft_forward(plan, col)) != RS_OK) goto done;
            rs_fft_shift(col, n_az);
            for (size_t f = 0; f < n_az; f++) {
                if (f < f0 || f >= f1) col[f] = 0.0f;
            }
            rs_ifft_shift(col, n_az);
            if ((st = rs_fft_inverse(plan, col)) != RS_OK) goto done;
            for (size_t a = 0; a < n_az; a++) o[a * n_rg + r] = col[a];
        }
    }

done:
    rs_fft_plan_destroy(plan);
    free(col);
    return st;
}

/* Split-band Phase Linking shift estimate. */
resonarsat_status_t rs_splitband_shift(const float complex *patch,
                                       size_t n_look, size_t n_az, size_t n_rg,
                                       double *shift_out,
                                       double *coherence_out)
{
    if (!patch || !shift_out || n_look < 2 || n_az < 6 || n_rg == 0)
        return RS_ERR_ARG;

    const size_t n_pix = n_az * n_rg;
    resonarsat_status_t st = RS_OK;

    double *psd = rs_pl_mean_spectrum(patch, n_look, n_az, n_rg);
    if (!psd) return RS_ERR_ALLOC;

    size_t lo = 0, hi = n_az - 1;
    rs_pl_occupied_band(psd, n_az, &lo, &hi);
    free(psd);

    const size_t width = hi - lo + 1;
    if (width < 6) {
        /* Too little occupied bandwidth to split three ways and retain anything
         * meaningful in each third. */
        rs_set_error("phaselink: occupied azimuth band is %zu bins, too narrow to split",
                     width);
        return RS_ERR_RANGE;
    }

    const size_t third = width / 3;
    const size_t lo0 = lo,              lo1 = lo + third;
    const size_t hi0 = hi + 1 - third,  hi1 = hi + 1;

    float complex *low  = malloc(n_look * n_pix * sizeof *low);
    float complex *high = malloc(n_look * n_pix * sizeof *high);
    double *pl = malloc(n_look * sizeof *pl);
    double *ph = malloc(n_look * sizeof *ph);
    if (!low || !high || !pl || !ph) { st = RS_ERR_ALLOC; goto done; }

    if ((st = rs_pl_bandpass(patch, n_look, n_az, n_rg, lo0, lo1, low)) != RS_OK) goto done;
    if ((st = rs_pl_bandpass(patch, n_look, n_az, n_rg, hi0, hi1, high)) != RS_OK) goto done;

    if ((st = rs_phase_link(low,  n_look, n_pix, pl)) != RS_OK) goto done;
    if ((st = rs_phase_link(high, n_look, n_pix, ph)) != RS_OK) goto done;

    /* Delay from the phase difference between sub-bands.
     *
     * The source's 3/(4*pi) assumes the split spans a signal of unit bandwidth,
     * so the sub-band centres sit 2/3 of the sampled band apart. Here the signal
     * occupies only 'width' of 'n_az' bins, so the centre separation is
     * (2/3)*(width/n_az) of the sampled band and the scaling grows by the
     * reciprocal of that occupancy. Omitting this factor would under-report
     * every shift by the oversampling ratio -- a factor of five on the
     * geometries used here, which would look like a working estimator with a
     * calibration error rather than an obvious failure. */
    const double occupancy = (double)width / (double)n_az;

    /* Note the sign. A signal delayed by +d samples acquires a spectral phase
     * of -2*pi*f*d, so the phase difference between the upper and lower
     * sub-bands runs OPPOSITE to the shift. Omitting the negation yields
     * magnitudes accurate to a couple of percent with every shift inverted --
     * which in a power spectrum is invisible, and which cost a round of
     * debugging to find here even with exact ground truth available. */
    const double scale = -3.0 / (4.0 * M_PI * occupancy);

    for (size_t k = 0; k < n_look; k++) {
        double d = ph[k] - pl[k];
        while (d >  M_PI) d -= 2.0 * M_PI;
        while (d < -M_PI) d += 2.0 * M_PI;
        shift_out[k] = d * scale;
    }
    shift_out[0] = 0.0;

    if (coherence_out) {
        /* Mean off-diagonal coherence of the full-band stack, as a quality
         * measure comparable with the correlation tracker's peak value. */
        double acc = 0.0;
        size_t n_pair = 0;
        for (size_t i = 0; i < n_look; i++) {
            for (size_t j = i + 1; j < n_look; j++) {
                const float complex *a = patch + i * n_pix;
                const float complex *b = patch + j * n_pix;
                double complex c = 0.0;
                double ea = 0.0, eb = 0.0;
                for (size_t p = 0; p < n_pix; p++) {
                    c += (double complex)a[p] * conj((double complex)b[p]);
                    const double ma = (double)cabsf(a[p]), mb = (double)cabsf(b[p]);
                    ea += ma * ma; eb += mb * mb;
                }
                const double nrm = sqrt(ea * eb);
                if (nrm > 0.0) { acc += cabs(c) / nrm; n_pair++; }
            }
        }
        *coherence_out = n_pair ? acc / (double)n_pair : 0.0;
    }

done:
    free(low); free(high); free(pl); free(ph);
    return st;
}
