/* Per-window vibration spectra from tracked displacement series. */

#include "resonarsat/microm.h"
#include "resonarsat/fft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Release every array a spectrum result owns. */
void rs_spectrum_free(rs_spectrum_t *s)
{
    if (!s) return;
    free(s->psd);
    free(s->freq);
    free(s->dominant_freq);
    free(s->amplitude);
    free(s->excursion_px);
    free(s->quality);
    free(s->prominence);
    memset(s, 0, sizeof *s);
}

/* Hann-windowed periodogram of each window's line-of-sight displacement series.
 *
 * Three choices worth stating. A least-squares straight line is removed first,
 * not merely the mean. Removing the mean alone leaves any linear trend in the
 * record, and a trend is not a small effect here: temporal phase unwrapping
 * performs a random walk when the phase is noisy, which produces exactly such a
 * ramp and puts all the energy in the lowest bin. Every window would then report
 * the same spurious "dominant frequency" of one bin width, which looks like a
 * measurement and is not.
 * A Hann window is applied because the record is short and unwindowed leakage
 * from a strong low-frequency component would otherwise swamp weaker peaks. And
 * the zero bin is excluded when the dominant frequency is chosen, since a
 * residual trend is not a vibration mode.
 *
 * A plain periodogram is the right estimator at this record length. The source
 * papers invoke convex and atomic-norm methods that resolve closely spaced
 * peaks better, but with a few dozen samples those mostly offer new ways to
 * read structure into noise; a higher-resolution estimator belongs here only
 * once the record length justifies it. */
resonarsat_status_t rs_spectrum_compute(const rs_microm_t *m,
                                       rs_spectrum_source_t source,
                                       rs_spectrum_t *out)
{
    return rs_spectrum_compute_band(m, source, 0.0, out);
}

/* Spectra with a band floor. See the header. */
resonarsat_status_t rs_spectrum_compute_band(const rs_microm_t *m,
                                            rs_spectrum_source_t source,
                                            double f_min,
                                            rs_spectrum_t *out)
{
    return rs_spectrum_compute_opts(m, source, f_min, RS_DETREND_LINEAR, out);
}

/* Spectra with the detrend selectable. See the header on why it is. */
resonarsat_status_t rs_spectrum_compute_opts(const rs_microm_t *m,
                                             rs_spectrum_source_t source,
                                             double f_min,
                                             rs_detrend_t detrend,
                                             rs_spectrum_t *out)
{
    /* Zeroed before any validation can fail, so that "check the status, then
     * free the struct" is safe on every path rather than on most of them.
     * rs_cphd_read() carried the non-uniform version of this and aborted the
     * first time a corpus made it fail on real data. */
    memset(out, 0, sizeof *out);

    if (!m || !out || !m->disp_los || !m->vel_los) return RS_ERR_ARG;
    if (m->n_looks < 4) {
        rs_set_error("spectrum: %zu looks is too few to estimate a spectrum", m->n_looks);
        return RS_ERR_ARG;
    }
    if (!(m->dt > 0.0)) {
        rs_set_error("spectrum: sub-look sampling interval is unset");
        return RS_ERR_MISSING_META;
    }

    const size_t n = m->n_looks;
    const size_t n_freq = n / 2 + 1;      /* one-sided spectrum */
    const double fs = 1.0 / m->dt;

    out->psd           = calloc(m->n_win * n_freq, sizeof *out->psd);
    out->freq          = calloc(n_freq, sizeof *out->freq);
    out->dominant_freq = calloc(m->n_win, sizeof *out->dominant_freq);
    out->amplitude     = calloc(m->n_win, sizeof *out->amplitude);
    out->excursion_px  = calloc(m->n_win, sizeof *out->excursion_px);
    out->quality       = calloc(m->n_win, sizeof *out->quality);
    out->prominence    = calloc(m->n_win, sizeof *out->prominence);
    if (!out->psd || !out->freq || !out->dominant_freq || !out->amplitude ||
        !out->quality || !out->prominence || !out->excursion_px) {
        rs_spectrum_free(out);
        return RS_ERR_ALLOC;
    }

    out->n_win = m->n_win;
    out->n_win_az = m->n_win_az;
    out->n_win_rg = m->n_win_rg;
    out->n_freq = n_freq;
    out->df = fs / (double)n;
    /* Carried through so rs_spectrum_best_window() can evaluate the floor
     * without the caller having to pass the tracking parameters again. */
    out->quant_px = m->quant_px;

    for (size_t k = 0; k < n_freq; k++) out->freq[k] = (double)k * out->df;

    /* Precompute the window and its power, so the periodogram is normalised
     * consistently regardless of record length. */
    double *win = malloc(n * sizeof *win);
    float complex *buf = malloc(n * sizeof *buf);
    rs_fft_plan *plan = NULL;
    resonarsat_status_t st = RS_OK;

    if (!win || !buf || (st = rs_fft_plan_create(n, &plan)) != RS_OK) {
        free(win); free(buf); rs_fft_plan_destroy(plan);
        rs_spectrum_free(out);
        return (st == RS_OK) ? RS_ERR_ALLOC : st;
    }

    double win_power = 0.0;
    for (size_t i = 0; i < n; i++) {
        win[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)(n - 1)));
        win_power += win[i] * win[i];
    }
    if (win_power <= 0.0) win_power = 1.0;

    const double *all = (source == RS_SPEC_DISPLACEMENT) ? m->disp_los : m->vel_los;

    /* Sums that depend only on the sample index, hoisted out of the window
     * loop: the abscissa is the same for every window. */
    const double sx  = 0.5 * (double)n * (double)(n - 1);
    const double sxx = (double)(n - 1) * (double)n * (2.0 * (double)(n - 1) + 1.0) / 6.0;
    const double det = (double)n * sxx - sx * sx;

    for (size_t w = 0; w < m->n_win; w++) {
        const double *series = all + w * n;

        /* Peak-to-peak of the TRACKED SHIFT, in pixels, before any detrending.
         * Not of 'series': that may be a velocity or a line-of-sight distance
         * depending on the observable, and the quantisation floor it is
         * compared against lives in pixels. Taken raw because detrending can
         * only shrink an excursion, and the question here is what the
         * correlator resolved, not what survives conditioning. */
        {
            double lo = m->disp_az[w * n], hi = lo;
            for (size_t i = 1; i < n; i++) {
                const double v = m->disp_az[w * n + i];
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
            out->excursion_px[w] = hi - lo;
        }

        /* Fit y = a + b*i and subtract as much of it as the caller asked for.
         * RS_DETREND_NONE leaves a == b == 0, so the same expression below
         * covers every case without branching inside the sample loop. */
        double sy = 0.0, sxy = 0.0;
        for (size_t i = 0; i < n; i++) {
            sy  += series[i];
            sxy += (double)i * series[i];
        }
        double a = 0.0, b = 0.0;
        if (detrend == RS_DETREND_LINEAR && det != 0.0) {
            b = ((double)n * sxy - sx * sy) / det;
            a = (sy - b * sx) / (double)n;
        } else if (detrend == RS_DETREND_LINEAR || detrend == RS_DETREND_MEAN) {
            a = sy / (double)n;
        }

        for (size_t i = 0; i < n; i++) {
            buf[i] = (float)((series[i] - (a + b * (double)i)) * win[i]);
        }

        if (rs_fft_forward(plan, buf) != RS_OK) continue;

        double *psd = out->psd + w * n_freq;
        for (size_t k = 0; k < n_freq; k++) {
            const double mag = (double)cabsf(buf[k]);
            /* Scale to a one-sided density: double the interior bins, and
             * normalise by the window power and the sampling rate. */
            const double scale = (k == 0 || (n % 2 == 0 && k == n_freq - 1)) ? 1.0 : 2.0;
            psd[k] = scale * mag * mag / (win_power * fs);
        }

        /* Dominant peak, skipping the zero bin and anything below the band
         * floor. k_min is computed once per window rather than hoisted only for
         * clarity; it costs nothing beside the transform. */
        size_t k_min = 1;
        if (f_min > 0.0 && out->df > 0.0) {
            const double kf = ceil(f_min / out->df);
            if (kf > 1.0) k_min = (size_t)kf;
            if (k_min >= n_freq) k_min = n_freq - 1;
        }
        size_t best = k_min;
        for (size_t k = k_min + 1; k < n_freq; k++) if (psd[k] > psd[best]) best = k;

        out->dominant_freq[w] = out->freq[best];
        out->amplitude[w] = sqrt(psd[best]);   /* qualitative -- see the header */
        out->quality[w] = m->quality[w];

        /* Prominence: peak power against the mean of the rest. A window holding
         * a vibrating target concentrates its energy in one bin; a window of
         * tracking noise spreads it evenly, giving a prominence near 1. */
        double sum = 0.0;
        for (size_t k = k_min; k < n_freq; k++) sum += psd[k];
        const double mean_power = (n_freq > k_min)
                                ? sum / (double)(n_freq - k_min) : 0.0;
        out->prominence[w] = (mean_power > 0.0) ? psd[best] / mean_power : 0.0;
    }

    rs_fft_plan_destroy(plan);
    free(win);
    free(buf);
    return RS_OK;
}


/* Window with the most prominent spectral peak, among those that actually
 * tracked something.
 *
 * The quality gate is not optional, and the failure it prevents is instructive.
 * On a scene where most windows contain empty background, those windows track
 * nothing; their displacement series is a slow drift, which concentrates almost
 * all its energy in the lowest frequency bin and therefore scores a HIGHER
 * prominence than a genuine vibrating target whose energy is spread by noise
 * across neighbouring bins. Ranking on prominence alone hands the answer to the
 * emptiest part of the scene, reporting a confident-looking sub-bin frequency
 * from a window with nothing in it.
 *
 * The gate is relative rather than absolute -- windows must reach half the best
 * tracking coherence present in this scene -- because absolute coherence depends
 * on scene content, and an isolated point target on empty background scores far
 * lower than distributed clutter however well it tracks. */
resonarsat_status_t rs_spectrum_best_window(const rs_spectrum_t *spec,
                                            size_t *out_window,
                                            double *out_prominence,
                                            size_t *out_n_candidates)
{
    if (out_n_candidates) *out_n_candidates = 0;
    if (!spec || !spec->prominence || !spec->quality || spec->n_win == 0)
        return RS_ERR_ARG;

    double q_max = 0.0;
    for (size_t w = 0; w < spec->n_win; w++) {
        if (spec->quality[w] > q_max) q_max = spec->quality[w];
    }
    const double q_min = 0.5 * q_max;

    /* Three sigma of the quantisation noise, in pixels. See the header for the
     * derivation; zero disables the floor for estimators it does not describe. */
    const double floor_px = (spec->quant_px > 0.0) ? 2.449 * spec->quant_px : 0.0;

    size_t best = spec->n_win;
    size_t n_cand = 0;
    for (size_t w = 0; w < spec->n_win; w++) {
        if (spec->quality[w] < q_min) continue;
        /* The coherence gate above is necessary and not sufficient: a window
         * that never moved correlates with itself perfectly and clears it at
         * the top. This is the gate that removes it. */
        if (floor_px > 0.0 && spec->excursion_px &&
            spec->excursion_px[w] < floor_px) continue;
        n_cand++;
        if (best == spec->n_win || spec->prominence[w] > spec->prominence[best]) best = w;
    }
    if (out_n_candidates) *out_n_candidates = n_cand;

    if (best == spec->n_win) {
        /* Nothing resolved motion above the tracker's own resolution. Reported
         * rather than papered over: the previous behaviour fell back to window
         * zero, which handed the caller a frequency and a prominence computed
         * from rounding noise and no way to tell. An absent measurement and a
         * wrong one are not the same answer. */
        rs_set_error("spectrum: no window resolved motion above the %g px "
                     "quantisation floor (3 sigma of 1/%g px); the scene may "
                     "hold no detectable motion, or the sub-pixel refinement "
                     "may be too coarse for the excursion present",
                     floor_px, (spec->quant_px > 0.0) ? 1.0 / spec->quant_px : 0.0);
        return RS_ERR_RANGE;
    }

    if (out_window)     *out_window = best;
    if (out_prominence) *out_prominence = spec->prominence[best];
    return RS_OK;
}
/* The frequency the most windows agree on. See microm.h. */
resonarsat_status_t rs_spectrum_consensus(const rs_spectrum_t *spec,
                                          double *out_freq,
                                          size_t *out_n_agree,
                                          size_t *out_n_distinct,
                                          size_t *out_n_voting,
                                          size_t *out_n_contiguous)
{
    if (out_freq)         *out_freq = 0.0;
    if (out_n_agree)      *out_n_agree = 0;
    if (out_n_distinct)   *out_n_distinct = 0;
    if (out_n_voting)     *out_n_voting = 0;
    if (out_n_contiguous) *out_n_contiguous = 0;
    if (!spec || !spec->prominence || !spec->quality || !spec->dominant_freq ||
        spec->n_win == 0)
        return RS_ERR_ARG;

    /* The same two gates rs_spectrum_best_window() applies, so that the counts
     * the two functions return describe one population. Duplicated rather than
     * factored out because the gates are the contract here, not an
     * implementation detail: a change to one that silently changed the other
     * would make the counts incomparable without anything saying so. */
    double q_max = 0.0;
    for (size_t w = 0; w < spec->n_win; w++) {
        if (spec->quality[w] > q_max) q_max = spec->quality[w];
    }
    const double q_min = 0.5 * q_max;
    const double floor_px = (spec->quant_px > 0.0) ? 2.449 * spec->quant_px : 0.0;

    /* Half a bin, so two windows agree exactly when they picked the same bin.
     * Falls back to an absolute tolerance if the axis carries no spacing. */
    const double tol = (spec->df > 0.0) ? 0.5 * spec->df : 1e-9;

    size_t n_vote = 0, best_count = 0, n_distinct = 0;
    double best_freq = 0.0;

    for (size_t w = 0; w < spec->n_win; w++) {
        if (spec->quality[w] < q_min) continue;
        if (floor_px > 0.0 && spec->excursion_px &&
            spec->excursion_px[w] < floor_px) continue;
        n_vote++;

        /* Count how many gated windows share this one's bin. O(n_win^2), which
         * is nothing beside the tracking that produced the spectrum: n_win is
         * tens to low thousands and each comparison is a subtraction. */
        size_t count = 0;
        int seen_earlier = 0;
        for (size_t v = 0; v < spec->n_win; v++) {
            if (spec->quality[v] < q_min) continue;
            if (floor_px > 0.0 && spec->excursion_px &&
                spec->excursion_px[v] < floor_px) continue;
            if (fabs(spec->dominant_freq[v] - spec->dominant_freq[w]) <= tol) {
                count++;
                /* Distinct values are counted once, at their first occurrence. */
                if (v < w) seen_earlier = 1;
            }
        }
        if (!seen_earlier) n_distinct++;

        /* Ties go to the lower frequency, which is arbitrary but fixed: an
         * unstable choice would make the count meaningful and the frequency
         * beside it not. */
        if (count > best_count ||
            (count == best_count && best_count > 0 &&
             spec->dominant_freq[w] < best_freq)) {
            best_count = count;
            best_freq  = spec->dominant_freq[w];
        }
    }

    if (out_n_voting)   *out_n_voting = n_vote;
    if (out_n_distinct) *out_n_distinct = n_distinct;

    if (n_vote == 0) {
        rs_set_error("spectrum: no window passed the coherence gate and the "
                     "quantisation floor, so there is nothing to take a "
                     "consensus over");
        return RS_ERR_RANGE;
    }

    if (out_freq)    *out_freq = best_freq;
    if (out_n_agree) *out_n_agree = best_count;

    /* Largest 4-connected block of agreeing windows on the window grid.
     *
     * Iterative flood fill with an explicit stack rather than recursion: n_win
     * runs to thousands and a scattered agreement set would otherwise recurse
     * as deep as it is large. Allocation failure leaves the count at zero and
     * the rest of the answer intact, since contiguity refines the result rather
     * than constituting it. */
    if (out_n_contiguous && spec->n_win_az > 0 && spec->n_win_rg > 0 &&
        spec->n_win_az * spec->n_win_rg == spec->n_win) {
        unsigned char *mark = calloc(spec->n_win, 1);
        size_t *stack = malloc(spec->n_win * sizeof *stack);
        if (mark && stack) {
            for (size_t w = 0; w < spec->n_win; w++) {
                if (spec->quality[w] < q_min) continue;
                if (floor_px > 0.0 && spec->excursion_px &&
                    spec->excursion_px[w] < floor_px) continue;
                if (fabs(spec->dominant_freq[w] - best_freq) <= tol) mark[w] = 1;
            }
            size_t largest = 0;
            for (size_t seed = 0; seed < spec->n_win; seed++) {
                if (mark[seed] != 1) continue;
                size_t top = 0, size = 0;
                stack[top++] = seed; mark[seed] = 2;
                while (top > 0) {
                    const size_t c = stack[--top];
                    size++;
                    const size_t ia = c / spec->n_win_rg;
                    const size_t ir = c % spec->n_win_rg;
                    /* four neighbours on the window lattice */
                    const long da[4] = { -1, 1, 0, 0 };
                    const long dr[4] = { 0, 0, -1, 1 };
                    for (int k = 0; k < 4; k++) {
                        const long na = (long)ia + da[k];
                        const long nr = (long)ir + dr[k];
                        if (na < 0 || nr < 0 ||
                            (size_t)na >= spec->n_win_az ||
                            (size_t)nr >= spec->n_win_rg) continue;
                        const size_t nb = (size_t)na * spec->n_win_rg + (size_t)nr;
                        if (mark[nb] == 1) { mark[nb] = 2; stack[top++] = nb; }
                    }
                }
                if (size > largest) largest = size;
            }
            *out_n_contiguous = largest;
        }
        free(mark); free(stack);
    }
    return RS_OK;
}


/* Observation ratio for a measured frequency. See the header on why this is a
 * diagnostic quantity and not a validity threshold. */
double rs_spectrum_observation_ratio(double t_sap, double freq)
{
    if (!(freq > 0.0)) return 0.0;
    return t_sap * freq;
}

/* Sinc response of the sub-aperture averaging window. See microm.h. */
double rs_spectrum_subaperture_response(double t_sap, double freq)
{
    const double x = M_PI * freq * t_sap;
    if (!(t_sap > 0.0) || !(freq > 0.0)) return 1.0;
    /* The small-argument branch avoids 0/0; sinc(x) -> 1 - x^2/6 there, and the
     * cutoff is where the series and the quotient agree to double precision. */
    if (fabs(x) < 1.0e-8) return 1.0;
    return fabs(sin(x) / x);
}
