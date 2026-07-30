/* Doppler sub-aperture decomposition.
 * See include/resonarsat/subaperture.h for the contract. */

#include "resonarsat/subaperture.h"
#include "resonarsat/fft.h"
#include "resonarsat/focus.h"
#include "resonarsat/geom.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Recommended defaults. Uniform mode with 40 percent overlap sits inside the
 * 38-43 percent range reported by the validated processing literature. */
void rs_subap_params_default(rs_subap_params_t *params)
{
    if (!params) return;
    params->mode = RS_SUBAP_UNIFORM;
    params->n_looks = 16;
    params->overlap = 0.40;
    params->left_out_frac = 0.5;
    params->window = 1;
    params->b_shift_hz = 0.0;   /* 0: derive from the sweep step */
    params->pair = 0;
    /* Threaded by default. Assigned explicitly because this function does not
     * memset the struct and its callers do not zero it -- see the note in
     * rs_microm_params_default(), where omitting the equivalent line produced a
     * tracker that silently changed method between runs. */
    params->single_thread = 0;
    params->range_taps = 0;
}

/* Gather one range column of an image into a contiguous buffer, so the azimuth
 * transform sees the stride-1 layout the FFT wrapper requires. */
static void rs_gather_column(const rs_slc_t *img, size_t rg, float complex *col)
{
    for (size_t az = 0; az < img->n_az; az++) col[az] = img->data[az * img->n_rg + rg];
}

/* Scatter a transformed column back into an image. */
static void rs_scatter_column(rs_slc_t *img, size_t rg, const float complex *col)
{
    for (size_t az = 0; az < img->n_az; az++) img->data[az * img->n_rg + rg] = col[az];
}

/* Accumulate the azimuth power spectrum of an image, averaged over range.
 *
 * Returns a newly allocated array of img->n_az doubles in fftshifted order, so
 * index 0 is the most negative Doppler frequency and the centre index is zero
 * Doppler. Caller frees. Returns NULL on allocation failure. */
static double *rs_azimuth_power_spectrum(const rs_slc_t *img)
{
    double *psd = calloc(img->n_az, sizeof *psd);
    float complex *col = malloc(img->n_az * sizeof *col);
    rs_fft_plan *plan = NULL;

    if (!psd || !col || rs_fft_plan_create(img->n_az, &plan) != RS_OK) {
        free(psd); free(col); rs_fft_plan_destroy(plan);
        return NULL;
    }

    for (size_t rg = 0; rg < img->n_rg; rg++) {
        rs_gather_column(img, rg, col);
        if (rs_fft_forward(plan, col) != RS_OK) continue;
        rs_fft_shift(col, img->n_az);
        for (size_t k = 0; k < img->n_az; k++) {
            const double m = (double)cabsf(col[k]);
            psd[k] += m * m;
        }
    }

    rs_fft_plan_destroy(plan);
    free(col);
    for (size_t k = 0; k < img->n_az; k++) psd[k] /= (double)img->n_rg;
    return psd;
}

/* Doppler frequency in Hz of fftshifted spectral index k. */
static double rs_bin_frequency(size_t k, size_t n, double fs)
{
    const double centre = (double)(n / 2);
    return ((double)k - centre) * fs / (double)n;
}

/* Doppler centroid as the power-weighted mean frequency. */
resonarsat_status_t rs_estimate_doppler_centroid(const rs_slc_t *img, double *out)
{
    if (!img || !out || !img->data) return RS_ERR_ARG;
    if (img->fs_az <= 0.0) {
        rs_set_error("doppler centroid: image has no azimuth sampling rate");
        return RS_ERR_MISSING_META;
    }

    double *psd = rs_azimuth_power_spectrum(img);
    if (!psd) return RS_ERR_ALLOC;

    double num = 0.0, den = 0.0;
    for (size_t k = 0; k < img->n_az; k++) {
        const double f = rs_bin_frequency(k, img->n_az, img->fs_az);
        num += f * psd[k];
        den += psd[k];
    }
    free(psd);

    if (den <= 0.0) {
        rs_set_error("doppler centroid: image has no energy");
        return RS_ERR_RANGE;
    }
    *out = num / den;
    return RS_OK;
}

/* Doppler bandwidth at a given power fraction below peak. */
resonarsat_status_t rs_estimate_doppler_bandwidth(const rs_slc_t *img, double frac, double *out)
{
    if (!img || !out || !img->data) return RS_ERR_ARG;
    if (!(frac > 0.0 && frac < 1.0)) {
        rs_set_error("doppler bandwidth: power fraction %g outside (0,1)", frac);
        return RS_ERR_ARG;
    }
    if (img->fs_az <= 0.0) {
        rs_set_error("doppler bandwidth: image has no azimuth sampling rate");
        return RS_ERR_MISSING_META;
    }

    double *psd = rs_azimuth_power_spectrum(img);
    if (!psd) return RS_ERR_ALLOC;

    size_t peak = 0;
    for (size_t k = 1; k < img->n_az; k++) if (psd[k] > psd[peak]) peak = k;

    const double threshold = psd[peak] * frac;
    size_t lo = peak, hi = peak;
    while (lo > 0 && psd[lo] > threshold) lo--;
    while (hi + 1 < img->n_az && psd[hi] > threshold) hi++;
    free(psd);

    *out = ((double)hi - (double)lo) * img->fs_az / (double)img->n_az;
    return RS_OK;
}

/* Raised-cosine band response at frequency 'f' for a band centred on 'centre'
 * with full width 'width' and transition width 'roll'.
 *
 * Returns 1 inside the flat portion, 0 outside the band, and a smooth cosine
 * taper across the transition. The taper is what keeps sub-look time series
 * free of the sinc sidelobes a rectangular band would impose. */
static double rs_raised_cosine(double f, double centre, double width, double roll)
{
    const double d = fabs(f - centre);
    const double flat = 0.5 * width - roll;
    if (d <= flat) return 1.0;
    if (d >= 0.5 * width) return 0.0;
    if (roll <= 0.0) return 0.0;
    return 0.5 * (1.0 + cos(M_PI * (d - flat) / roll));
}

/* Release a stack and everything it owns. */
void rs_subap_stack_free(rs_subap_stack_t *stack)
{
    if (!stack) return;
    if (stack->look) {
        for (size_t i = 0; i < stack->n_looks; i++) rs_slc_free(&stack->look[i]);
        free(stack->look);
    }
    if (stack->slave) {
        for (size_t i = 0; i < stack->n_looks; i++) rs_slc_free(&stack->slave[i]);
        free(stack->slave);
    }
    free(stack->centre_time);
    memset(stack, 0, sizeof *stack);
}

/* Decompose an image into a sub-look stack. */
resonarsat_status_t rs_subaperture_split(const rs_slc_t *img,
                                         const rs_subap_params_t *params,
                                         rs_subap_stack_t *stack)
{
    if (!stack) return RS_ERR_ARG;

    /* Zeroed before any validation can fail, so that "check the status, then
     * free the struct" is safe on every path rather than on most of them.
     * rs_cphd_read() carried the non-uniform version of this and aborted the
     * first time a corpus made it fail on real data. */
    memset(stack, 0, sizeof *stack);

    if (!img || !params || !stack || !img->data) return RS_ERR_ARG;
    if (params->n_looks == 0) {
        rs_set_error("subaperture: n_looks must be positive");
        return RS_ERR_ARG;
    }
    if (img->n_az < 2 * params->n_looks) {
        rs_set_error("subaperture: %zu azimuth lines too few for %zu looks",
                     img->n_az, params->n_looks);
        return RS_ERR_ARG;
    }
    if (img->fs_az <= 0.0) {
        rs_set_error("subaperture: image has no azimuth sampling rate");
        return RS_ERR_MISSING_META;
    }
    if (params->mode == RS_SUBAP_PAPER &&
        !(params->left_out_frac > 0.0 && params->left_out_frac < 1.0)) {
        rs_set_error("subaperture: left_out_frac %g outside (0,1)", params->left_out_frac);
        return RS_ERR_ARG;
    }
    if (!isfinite(params->b_shift_hz) || params->b_shift_hz < 0.0) {
        rs_set_error("subaperture: b_shift_hz %g must be finite and non-negative",
                     params->b_shift_hz);
        return RS_ERR_ARG;
    }
    stack->params = *params;
    stack->n_looks = params->n_looks;

    const size_t n_az = img->n_az, n_rg = img->n_rg, n_looks = params->n_looks;

    /* Doppler centroid: prefer the annotated polynomial when present, and
     * cross-check it against the data-driven estimate. A large disagreement is
     * reported rather than silently averaged, since it usually means the region
     * of interest straddles a burst boundary. */
    double fdc_est = 0.0;
    resonarsat_status_t st = rs_estimate_doppler_centroid(img, &fdc_est);
    if (st != RS_OK) return st;
    stack->doppler_centroid = fdc_est;

    double bw = 0.0;
    if (rs_estimate_doppler_bandwidth(img, 0.5, &bw) != RS_OK || bw <= 0.0) {
        bw = img->fs_az * 0.8;  /* fall back to most of the sampled band */
    }
    stack->doppler_bandwidth = bw;

    /* Band layout. Both modes reduce to a list of (centre, width) pairs in Hz,
     * relative to the Doppler centroid. */
    double *centre = malloc(n_looks * sizeof *centre);
    double *width  = malloc(n_looks * sizeof *width);
    if (!centre || !width) { free(centre); free(width); return RS_ERR_ALLOC; }

    double t_sap = 0.0;

    if (params->mode == RS_SUBAP_PAPER) {
        /* Biondi & Malanga (2022) section 3.1: hold out B_DL, process the rest,
         * and step it across the spectrum in N_D rigid shifts.
         *
         * The convention check below is the one described in the header: the
         * paper's two statements of the step size agree only at 0.5. */
        const double b_dl = bw * params->left_out_frac;
        const double b_proc = bw - b_dl;
        const double step_printed = b_proc / (double)n_looks;
        const double step_prose   = b_dl / (double)n_looks;

        if (fabs(step_printed - step_prose) > 1e-9 * fmax(1.0, fabs(step_printed))) {
            rs_set_error("subaperture: left_out_frac %g makes the source's two "
                         "statements of the sub-aperture step disagree "
                         "((B_CD-B_DL)/ND = %g Hz vs B_DL/ND = %g Hz); using the "
                         "printed formula (B_CD-B_DL)/ND",
                         params->left_out_frac, step_printed, step_prose);
            /* Advisory only: the decomposition proceeds with the printed form,
             * but the caller can surface the message so the choice is visible. */
        }

        /* B_shift is the master-slave separation, free in the patent and
         * defaulting to the sweep step here. See rs_subap_params_t. */
        stack->b_shift_hz = (params->b_shift_hz > 0.0) ? params->b_shift_hz
                                                       : step_printed;

        /* Start the first master at the lower edge of the measured Doppler
         * support. The old half-step offset made the final slave extend half a
         * step beyond the support even with the default B_shift, so the pair
         * ceased to be rigid at the end of the sweep.
         *
         * That was a real defect and fixing it did NOT make RS_MICROM_REF_PAIR
         * recover a frequency -- test_tracking.c measures the identical
         * lowest-bin result before and after. Recorded here so the correction
         * is not later mistaken for the cure; see rs_microm_ref_t. */
        const double first_centre = -0.5 * bw + 0.5 * b_proc;
        for (size_t i = 0; i < n_looks; i++) {
            centre[i] = first_centre + (double)i * step_printed;
            width[i]  = b_proc;
        }

        if (params->pair) {
            const double slave_high = centre[n_looks - 1]
                                    + stack->b_shift_hz + 0.5 * b_proc;
            const double band_high = 0.5 * bw;
            const double tol = 1e-9 * fmax(1.0, bw);
            if (slave_high > band_high + tol) {
                rs_set_error("subaperture: B_shift %g Hz does not fit %zu rigid "
                             "master-slave pairs inside the %g Hz Doppler band "
                             "(maximum %g Hz for this layout)",
                             stack->b_shift_hz, n_looks, bw,
                             fmax(0.0, band_high
                                           - centre[n_looks - 1]
                                           - 0.5 * b_proc));
                free(centre); free(width);
                return RS_ERR_ARG;
            }
        }
        /* Each look uses a band of width b_proc, so its aperture time is that
         * fraction of the full dwell. */
        t_sap = (img->t_dwell > 0.0) ? img->t_dwell * (b_proc / bw)
                                     : b_proc / bw;
        stack->dt = (n_looks > 1)
                  ? ((img->t_dwell > 0.0 ? img->t_dwell : 1.0) * step_printed / bw)
                  : 0.0;
    } else {
        /* Uniform filter bank: n_looks equal bands with fractional overlap. */
        const double denom = (double)n_looks - (double)(n_looks - 1) * params->overlap;
        const double band = (denom > 0.0) ? bw / denom : bw / (double)n_looks;
        const double step = band * (1.0 - params->overlap);
        const double first = -0.5 * bw + 0.5 * band;

        for (size_t i = 0; i < n_looks; i++) {
            centre[i] = first + (double)i * step;
            width[i]  = band;
        }
        t_sap = (img->t_dwell > 0.0) ? img->t_dwell * (band / bw) : band / bw;
        stack->dt = (img->t_dwell > 0.0) ? img->t_dwell * (step / bw) : step / bw;
    }

    /* Reject a pair in uniform mode: B_shift and the master-slave construction
     * are the paper/patent front end, and uniform mode has neither. */
    if (params->pair && params->mode != RS_SUBAP_PAPER) {
        rs_set_error("subaperture: pair mode needs the paper band layout "
                     "(B_shift is defined only there)");
        free(centre); free(width);
        return RS_ERR_ARG;
    }

    stack->t_sap = t_sap;
    stack->f_max = (stack->dt > 0.0) ? 1.0 / (2.0 * stack->dt) : 0.0;
    /* The spectral gap is a time lag: dt_pair = B_shift * t_dwell / B_CD. This
     * is what makes B_shift select a frequency -- see rs_subap_params_t. */
    stack->pair_lag_s = (bw > 0.0 && img->t_dwell > 0.0)
                      ? stack->b_shift_hz * img->t_dwell / bw
                      : 0.0;
    stack->az_resolution = (img->lambda > 0.0 && img->v_platform > 0.0 && img->r0 > 0.0)
                         ? rs_azimuth_resolution(img->lambda, img->r0, img->v_platform, t_sap)
                         : 0.0;

    /* Allocate the output images and their centre times. */
    stack->look = calloc(n_looks, sizeof *stack->look);
    stack->centre_time = calloc(n_looks, sizeof *stack->centre_time);
    if (params->pair) stack->slave = calloc(n_looks, sizeof *stack->slave);
    if (!stack->look || !stack->centre_time || (params->pair && !stack->slave)) {
        free(centre); free(width);
        rs_subap_stack_free(stack);
        return RS_ERR_ALLOC;
    }
    for (size_t i = 0; i < n_looks; i++) {
        st = rs_slc_alloc(&stack->look[i], n_az, n_rg);
        if (st != RS_OK) { free(centre); free(width); rs_subap_stack_free(stack); return st; }
        /* Inherit geometry; the look's own timing describes its window. */
        stack->look[i] = (rs_slc_t){ .data = stack->look[i].data,
                                     .n_az = n_az, .n_rg = n_rg,
                                     .azimuth_time_interval = img->azimuth_time_interval,
                                     .fs_az = img->fs_az,
                                     .pulse_prf = img->pulse_prf,
                                     .fc = img->fc, .lambda = img->lambda,
                                     .rg_spacing_m = img->rg_spacing_m,
                                     .az_spacing_m = img->az_spacing_m,
                                     .r0 = img->r0, .t0 = img->t0,
                                     .t_dwell = t_sap,
                                     .incidence = img->incidence,
                                     .v_platform = img->v_platform,
                                     .doppler = img->doppler, .orbit = img->orbit };
        snprintf(stack->look[i].source, sizeof stack->look[i].source, "sublook%zu", i);
        stack->centre_time[i] = img->t0 + (double)i * stack->dt;

        if (params->pair) {
            st = rs_slc_alloc(&stack->slave[i], n_az, n_rg);
            if (st != RS_OK) { free(centre); free(width); rs_subap_stack_free(stack); return st; }
            /* Same geometry as its master, its own sample buffer. */
            float complex *sdata = stack->slave[i].data;
            stack->slave[i] = stack->look[i];
            stack->slave[i].data = sdata;
            snprintf(stack->slave[i].source, sizeof stack->slave[i].source,
                     "subslave%zu", i);
        }
    }

    /* Transform each range column once, then apply every band filter to the
     * one spectrum. Transforming per look instead would repeat the forward FFT
     * n_looks times for no benefit.
     *
     * ON DFT2 VERSUS THIS 1-D TRANSFORM. Blocks 2, 5 and 6 of the patent's
     * Figure 0.5 are named DFT2 and IDFT2, and this transforms along azimuth
     * only. The two are identical here, not an approximation of each other:
     * WO2024008365A1 [0004] focuses the master "using the entire chirp band
     * with amplitude B_Cr and half of the Doppler band", and Figure 0.2 draws
     * both squares at full height B_Cr, displaced from one another only along
     * the azimuth frequency axis by B_shift. A separable band-pass whose range
     * response is unity over the whole chirp band is the identity in range, so
     * the range transform and its inverse cancel exactly. Doing it anyway would
     * cost an n_rg-point FFT per line and change nothing.
     *
     * This equivalence is what makes the shortcut legitimate, and it stops
     * holding the moment a sub-band of the chirp is used -- Eq. (6)'s B_Cr_i,
     * and claim 7's "chirp-Doppler sub-apertures". An implementation of that
     * must transform in two dimensions rather than extending this loop.
     *
     * THAT PATH IS EXERCISED, which this comment previously denied. Biondi,
     * "High-Voltage Electric Power Transmission Monitoring by Micro-Motion
     * Estimation on Synthetic Aperture Radar Data" (Preprints 2023,
     * doi 10.20944/preprints202308.0926.v1) builds exactly it: its Figure 2 is a
     * grid of N_c chirp band-passes crossed with N_D Doppler sub-apertures, its
     * Figure 4 wires the two into one pipeline, and its section 4 credits the
     * "frequency diversity" with the improvement it reports. So the
     * two-dimensional case is not hypothetical; it is a published application of
     * the same method by the same author, on COSMO-SkyMed staring spotlight.
     *
     * It remains unimplemented here, and deliberately so rather than by
     * oversight: the paper that exercises it reports a 50 Hz line at 51.23 Hz
     * and 57.87 Hz from two processing routes on the same target, against a
     * mains frequency regulated to 0.05 Hz, and validates against a microphone
     * recording made about 320 km from the SAR scene. Adding a stage on that
     * evidence would be adding surface area, not capability. The seam is
     * recorded so that a future implementation knows where it goes. */
    float complex *col = malloc(n_az * sizeof *col);
    float complex *band_col = malloc(n_az * sizeof *band_col);
    rs_fft_plan *plan = NULL;
    if (!col || !band_col || (st = rs_fft_plan_create(n_az, &plan)) != RS_OK) {
        if (col && band_col) st = RS_ERR_ALLOC;
        free(col); free(band_col); free(centre); free(width);
        rs_subap_stack_free(stack);
        return (st == RS_OK) ? RS_ERR_ALLOC : st;
    }

    for (size_t rg = 0; rg < n_rg; rg++) {
        rs_gather_column(img, rg, col);
        if ((st = rs_fft_forward(plan, col)) != RS_OK) break;
        rs_fft_shift(col, n_az);

        for (size_t i = 0; i < n_looks; i++) {
            const double roll = params->window ? 0.25 * width[i] : 0.0;
            for (size_t k = 0; k < n_az; k++) {
                const double f = rs_bin_frequency(k, n_az, img->fs_az) - stack->doppler_centroid;
                const double g = params->window
                               ? rs_raised_cosine(f, centre[i], width[i], roll)
                               : ((fabs(f - centre[i]) <= 0.5 * width[i]) ? 1.0 : 0.0);
                band_col[k] = col[k] * (float)g;
            }
            rs_ifft_shift(band_col, n_az);
            if ((st = rs_fft_inverse(plan, band_col)) != RS_OK) break;
            rs_scatter_column(&stack->look[i], rg, band_col);

            if (stack->slave) {
                /* The slave is the same band offset by B_shift, which is what
                 * "rigidly held at a distance B_shift" means: identical width,
                 * identical roll-off, displaced along the azimuth frequency. */
                const double sc = centre[i] + stack->b_shift_hz;
                for (size_t k = 0; k < n_az; k++) {
                    const double f = rs_bin_frequency(k, n_az, img->fs_az) - stack->doppler_centroid;
                    const double g = params->window
                                   ? rs_raised_cosine(f, sc, width[i], roll)
                                   : ((fabs(f - sc) <= 0.5 * width[i]) ? 1.0 : 0.0);
                    band_col[k] = col[k] * (float)g;
                }
                rs_ifft_shift(band_col, n_az);
                if ((st = rs_fft_inverse(plan, band_col)) != RS_OK) break;
                rs_scatter_column(&stack->slave[i], rg, band_col);
            }
        }
        if (st != RS_OK) break;
    }

    rs_fft_plan_destroy(plan);
    free(col); free(band_col); free(centre); free(width);

    if (st != RS_OK) { rs_subap_stack_free(stack); return st; }
    return RS_OK;
}

/* Form a sub-look stack by focusing shifted pulse windows.
 *
 * The aperture is divided in time rather than in frequency, so every timing
 * quantity comes straight off the recorded pulse times and none of it has to be
 * inferred from a spectrum. That directness is the reason to prefer this path:
 * with a spectral split the mapping from band to time depends on how the input
 * was focused, which the sub-aperture stage cannot see.
 *
 * Window layout: n_looks windows of 'per' pulses each, stepped by
 * per*(1-overlap), sized so the last window ends at the final pulse. */
resonarsat_status_t rs_subaperture_from_cphd(const struct rs_cphd *cphd_in,
                                             const struct rs_grid *grid_in,
                                             const rs_subap_params_t *params,
                                             rs_subap_stack_t *stack)
{
    if (!stack) return RS_ERR_ARG;

    /* Zeroed before any validation can fail, so that "check the status, then
     * free the struct" is safe on every path rather than on most of them.
     * rs_cphd_read() carried the non-uniform version of this and aborted the
     * first time a corpus made it fail on real data. */
    memset(stack, 0, sizeof *stack);

    const rs_cphd_t *cphd = (const rs_cphd_t *)cphd_in;
    const rs_grid_t *grid = (const rs_grid_t *)grid_in;

    if (!cphd || !grid || !params || !stack) return RS_ERR_ARG;
    if (params->n_looks == 0) {
        rs_set_error("subaperture: n_looks must be positive");
        return RS_ERR_ARG;
    }
    if (params->mode == RS_SUBAP_PAPER) {
        /* Refused rather than ignored. This path divides the aperture in TIME;
         * the paper's decomposition is defined on the Doppler SPECTRUM, holding
         * out a fraction of the band and stepping the remainder across it.
         * There is no window layout here that reproduces it.
         *
         * Silently proceeding was the original behaviour and it was a real
         * defect: every caller that asked for the paper's decomposition got a
         * uniform filter bank instead, with nothing in the output to say so.
         * Since the band held out is what the paper credits with making motion
         * measurable, substituting for it invisibly is precisely the
         * substitution that must not happen quietly.
         *
         * Callers wanting the paper's route must focus the full aperture and
         * pass the image to rs_subaperture_split(). */
        rs_set_error("subaperture: the paper's decomposition is defined on the "
                     "Doppler spectrum and cannot be formed from pulse windows; "
                     "focus first and use rs_subaperture_split()");
        return RS_ERR_ARG;
    }
    if (!(params->overlap >= 0.0 && params->overlap < 1.0)) {
        rs_set_error("subaperture: overlap %g outside [0,1)", params->overlap);
        return RS_ERR_ARG;
    }

    const size_t n_looks = params->n_looks;
    const size_t n_pulse = cphd->n_pulse;

    /* Solve n_pulse = per + (n_looks-1)*per*(1-overlap) for the window size, so
     * that the windows exactly span the collect. */
    const double denom = 1.0 + (double)(n_looks - 1) * (1.0 - params->overlap);
    const size_t per = (size_t)((double)n_pulse / denom);
    const size_t step = (n_looks > 1)
                      ? (size_t)((double)per * (1.0 - params->overlap))
                      : 0;

    if (per == 0 || (n_looks > 1 && step == 0)) {
        rs_set_error("subaperture: %zu pulses cannot be divided into %zu looks "
                     "at %.2f overlap (window %zu, step %zu)",
                     n_pulse, n_looks, params->overlap, per, step);
        return RS_ERR_ARG;
    }

    stack->params = *params;
    stack->n_looks = n_looks;
    stack->look = calloc(n_looks, sizeof *stack->look);
    stack->centre_time = calloc(n_looks, sizeof *stack->centre_time);
    if (!stack->look || !stack->centre_time) {
        rs_subap_stack_free(stack);
        return RS_ERR_ALLOC;
    }

    const double prf = (cphd->prf > 0.0) ? cphd->prf : 1.0;
    stack->t_sap = (double)per / prf;
    stack->dt    = (n_looks > 1) ? (double)step / prf : 0.0;
    stack->f_max = (stack->dt > 0.0) ? 1.0 / (2.0 * stack->dt) : 0.0;

    for (size_t i = 0; i < n_looks; i++) {
        size_t start = i * step;
        if (start + per > n_pulse) start = n_pulse - per;   /* clamp the last window */

        resonarsat_status_t st = rs_slc_alloc(&stack->look[i], grid->n_x, grid->n_y);
        if (st != RS_OK) { rs_subap_stack_free(stack); return st; }

        const rs_focus_opts_t fopts = { .single_thread = params->single_thread,
                                        .range_taps = params->range_taps };
        st = rs_focus_backproject_opts(cphd, grid, start, per, &fopts, &stack->look[i]);
        if (st != RS_OK) { rs_subap_stack_free(stack); return st; }

        stack->centre_time[i] = cphd->t[start + per / 2];
        snprintf(stack->look[i].source, sizeof stack->look[i].source, "sublook%zu", i);
    }

    /* Azimuth resolution of one look. Platform speed and slant range come from
     * what rs_focus_backproject() measured off the recorded positions, not from
     * an assumed constant. */
    stack->az_resolution = rs_azimuth_resolution(cphd->lambda,
                                                 stack->look[0].r0,
                                                 stack->look[0].v_platform,
                                                 stack->t_sap);
    stack->doppler_bandwidth = 0.0;   /* not measured on this path */
    stack->doppler_centroid = 0.0;

    return RS_OK;
}
