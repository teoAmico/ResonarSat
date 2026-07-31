/* Micro-motion extraction by sub-pixel offset tracking across a sub-look stack.
 * See include/resonarsat/microm.h for the contract and the accuracy envelope. */

#include "resonarsat/microm.h"
#include "resonarsat/coreg.h"
#include "resonarsat/phaselink.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Published working defaults. The patch size is a mid-range value from the
 * 51-131 pixel span reported in the validated literature, not a universal
 * optimum: too small a patch makes the tracker systematically underestimate
 * displacement, and the right size depends on the target's brightness and
 * apparent motion. Sweep it on real data. */
void rs_microm_params_default(rs_microm_params_t *params)
{
    if (!params) return;
    /* FIRST, not ADJACENT, despite ADJACENT correlating far better. Measured:
     * ADJACENT doubles per-window coherence (0.20 -> 0.50) and removes the
     * second-harmonic artefact, but accumulating its differentials integrates
     * tracking noise into a random walk whose low-frequency energy swamps the
     * signal -- on the single-target frequency sweep it takes recovery from
     * 5 of 6 to 0 of 6, with every window reporting the lowest spectral bin.
     * Until that accumulation is handled, FIRST is the better default. */
    params->estimator = RS_MICROM_EST_CORRELATION;
    params->reference = RS_MICROM_REF_FIRST;
    params->ref_lag = 1;
    params->win_az = 64;
    params->win_rg = 64;
    params->stride_az = 16;
    params->stride_rg = 16;
    params->upsample_az = 10;
    params->upsample_rg = 20;
    /* 0.4 is the literature value for distributed scenes. It is a real filter:
     * see rs_microm_track() for why sub-threshold windows are zeroed, and note
     * that synthetic point-target scenes score well below it. */
    params->coherence_min = 0.4;
    /* Off by default, on the same measured grounds as the reference mode above.
     * Common-mode removal is standard practice and sound where most of a scene
     * is genuinely static clutter, but the median it subtracts is only as good
     * as the windows it is taken over. On a scene dominated by EMPTY background
     * -- which every synthetic fixture here is -- most windows track noise, the
     * median is itself noisy, and subtracting it injects noise into the one
     * window that had signal: recovery falls from 5 of 6 to 3 of 6. Enable it
     * on real scenes with distributed clutter, where the premise holds. */
    params->remove_common_mode = 0;
    /* Optimised path by default.
     *
     * ASSIGNED EXPLICITLY, LIKE EVERY FIELD ABOVE, AND THAT IS LOAD-BEARING.
     * This function does not memset the struct, and its callers declare it
     * uninitialised on the stack -- so a field added to rs_microm_params_t and
     * not assigned here holds whatever was there. When this one was first added
     * and left out, the tracker read a garbage flag and took the exhaustive
     * correlation path on some runs and not others: test_nullmotion lost its
     * detections and test_tracking ran for over ten minutes, with nothing in
     * either output pointing at the cause. Any new field must be assigned here. */
    params->no_optimize = 0;
}

/* Median of a scratch array, which this reorders. */
static int rs_cmp_dbl(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Median of 'v', which this REORDERS in place. Callers must pass scratch they
 * do not need preserved. Returns 0 for an empty array. */
static double rs_median(double *v, size_t n)
{
    if (n == 0) return 0.0;
    qsort(v, n, sizeof *v, rs_cmp_dbl);
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* Subtract the scene-median shift from every window, look by look.
 *
 * See rs_microm_params_t.remove_common_mode for why this is not cosmetic. The
 * line-of-sight velocities are recomputed from the corrected shifts rather than
 * corrected separately, so the two cannot drift apart. */
static void rs_microm_remove_common_mode(rs_microm_t *m, double az_spacing,
                                         double v_plat, double r_slant,
                                         int can_velocity)
{
    double *scratch = malloc(m->n_win * sizeof *scratch);
    if (!scratch) return;

    for (size_t k = 0; k < m->n_looks; k++) {
        for (size_t w = 0; w < m->n_win; w++) scratch[w] = m->disp_az[w * m->n_looks + k];
        const double med_az = rs_median(scratch, m->n_win);

        for (size_t w = 0; w < m->n_win; w++) scratch[w] = m->disp_rg[w * m->n_looks + k];
        const double med_rg = rs_median(scratch, m->n_win);

        for (size_t w = 0; w < m->n_win; w++) {
            const size_t idx = w * m->n_looks + k;
            m->disp_az[idx] -= med_az;
            m->disp_rg[idx] -= med_rg;
            if (can_velocity) {
                m->vel_los[idx] = m->disp_az[idx] * az_spacing * v_plat / r_slant;
            }
        }
    }
    free(scratch);
}

/* Release every array a micro-motion result owns. */
void rs_microm_free(rs_microm_t *m)
{
    if (!m) return;
    free(m->disp_az);
    free(m->disp_rg);
    free(m->vel_los);
    free(m->disp_los);
    free(m->phase);
    free(m->quality);
    memset(m, 0, sizeof *m);
}

/* Track all windows across the stack. */
resonarsat_status_t rs_microm_track(const rs_subap_stack_t *stack,
                                    const rs_microm_params_t *params,
                                    rs_microm_t *out)
{
    if (!out) return RS_ERR_ARG;

    /* Zeroed before any validation can fail, so that "check the status, then
     * free the struct" is safe on every path rather than on most of them.
     * rs_cphd_read() carried the non-uniform version of this and aborted the
     * first time a corpus made it fail on real data. */
    memset(out, 0, sizeof *out);

    if (!stack || !params || !out || !stack->look) return RS_ERR_ARG;
    if (stack->n_looks < 2) {
        rs_set_error("microm: need at least 2 sub-looks, have %zu", stack->n_looks);
        return RS_ERR_ARG;
    }

    if (params->reference == RS_MICROM_REF_PAIR) {
        if (!stack->slave) {
            rs_set_error("microm: reference 'pair' needs a stack built with "
                         "rs_subap_params_t.pair (the slave bands are absent)");
            return RS_ERR_ARG;
        }
        if (params->estimator != RS_MICROM_EST_CORRELATION) {
            rs_set_error("microm: reference 'pair' measures a master-slave "
                         "offset, so it requires the correlation estimator");
            return RS_ERR_ARG;
        }
    }

    const rs_slc_t *ref_img = &stack->look[0];
    const size_t win_az = params->win_az, win_rg = params->win_rg;
    if (win_az == 0 || win_rg == 0 ||
        win_az > ref_img->n_az || win_rg > ref_img->n_rg) {
        rs_set_error("microm: window %zux%zu does not fit image %zux%zu",
                     win_az, win_rg, ref_img->n_az, ref_img->n_rg);
        return RS_ERR_ARG;
    }

    /* Reject an oversized exhaustive surface HERE, not in the tracking loop.
     *
     * The correlator's size ceiling depends only on the window and the upsample
     * factors, so if it is going to be exceeded it is exceeded for every window
     * identically. The loop body treats a failed correlation as "this window did
     * not track" and moves on, which is right for a degenerate patch and wrong
     * for this: it would return a full, well-formed result in which every single
     * window is zero, with no indication why. That is the project's signature
     * failure mode -- a plausible empty product -- so the condition is turned into
     * a refusal before anything is allocated. */
    if (params->no_optimize && params->estimator != RS_MICROM_EST_PHASE) {
        const resonarsat_status_t size_st =
            rs_coreg_surface_check(win_az, win_rg,
                                   params->upsample_az, params->upsample_rg);
        if (size_st != RS_OK) return size_st;
    }

    const size_t stride_az = params->stride_az ? params->stride_az : win_az;
    const size_t stride_rg = params->stride_rg ? params->stride_rg : win_rg;

    const size_t n_win_az = (ref_img->n_az - win_az) / stride_az + 1;
    const size_t n_win_rg = (ref_img->n_rg - win_rg) / stride_rg + 1;
    const size_t n_win = n_win_az * n_win_rg;
    const size_t n_looks = stack->n_looks;

    out->disp_az  = calloc(n_win * n_looks, sizeof *out->disp_az);
    out->disp_rg  = calloc(n_win * n_looks, sizeof *out->disp_rg);
    out->vel_los  = calloc(n_win * n_looks, sizeof *out->vel_los);
    out->disp_los = calloc(n_win * n_looks, sizeof *out->disp_los);
    out->phase    = calloc(n_win * n_looks, sizeof *out->phase);
    out->quality  = calloc(n_win, sizeof *out->quality);
    if (!out->disp_az || !out->disp_rg || !out->vel_los || !out->disp_los ||
        !out->phase || !out->quality) {
        rs_microm_free(out);
        return RS_ERR_ALLOC;
    }

    out->n_win_az = n_win_az;
    out->n_win_rg = n_win_rg;
    out->n_win = n_win;
    out->n_looks = n_looks;
    out->win_az = win_az;
    out->win_rg = win_rg;
    out->stride_az = stride_az;
    out->stride_rg = stride_rg;
    out->dt = stack->dt;
    out->f_max = stack->f_max;
    out->az_spacing_m = ref_img->az_spacing_m;
    out->rg_spacing_m = ref_img->rg_spacing_m;
    /* The correlator locates its peak on a grid of 1/upsample_az pixels, so
     * that is the finest excursion it can report and the floor below which a
     * series is rounding rather than motion. The phase estimator does not use
     * the correlation surface at all, so the concept does not apply and the
     * field stays zero -- see rs_microm_t.quant_px. */
    out->quant_px = (params->estimator == RS_MICROM_EST_PHASE || params->upsample_az == 0)
                  ? 0.0 : 1.0 / (double)params->upsample_az;

    const double lambda = (ref_img->lambda > 0.0) ? ref_img->lambda : 0.031;

    /* Geometry for converting a tracked azimuth shift into a line-of-sight
     * velocity: v_r = dx * V / R. All three come from what the focuser measured
     * off the recorded pulse positions. If any is missing the conversion is
     * skipped rather than guessed -- a plausible-looking velocity computed from
     * an assumed platform speed would be worse than none. */
    const double az_spacing = ref_img->az_spacing_m;
    const double v_plat     = ref_img->v_platform;
    const double r_slant    = ref_img->r0;
    const int    can_velocity = (az_spacing > 0.0 && v_plat > 0.0 && r_slant > 0.0);

    volatile resonarsat_status_t shared_st = RS_OK;

    /* Which extent the correlator searches. The exhaustive baseline is the half
     * of --no-optimize that can actually move a number; see
     * rs_microm_params_t.no_optimize. */
    const rs_coreg_refine_t refine = params->no_optimize
                                   ? RS_COREG_REFINE_EXHAUSTIVE
                                   : RS_COREG_REFINE_LOCAL;

    /* Windows are independent, so this is the pipeline's most naturally
     * parallel loop. Each thread allocates its own patch buffers; sharing them
     * would serialise the whole thing.
     *
     * Under --no-optimize the region is suppressed by the 'if' clause rather than
     * by a second copy of the loop, for the reason given in
     * rs_focus_backproject_opts(): one body means the arithmetic cannot drift
     * between the two modes. Note that 'dynamic' scheduling below already makes
     * the ORDER in which windows are visited unpredictable between runs -- and
     * that this is harmless, because each window writes only its own slots and
     * reads nothing another window wrote. The serial mode makes that visitation
     * order fixed, which is what allows a debugger or a profiler to be pointed at
     * a particular window. */
#ifdef _OPENMP
#pragma omp parallel if (!params->no_optimize)
#endif
    {
        float complex *pref = malloc(win_az * win_rg * sizeof *pref);
        float complex *pcur = malloc(win_az * win_rg * sizeof *pcur);
        /* Split-band estimation needs every look's patch simultaneously. */
        const int want_stack = (params->estimator == RS_MICROM_EST_SPLITBAND);
        float complex *pstack = want_stack
            ? malloc(n_looks * win_az * win_rg * sizeof *pstack) : NULL;
        double *pshift = want_stack ? malloc(n_looks * sizeof *pshift) : NULL;

        if (!pref || !pcur || (want_stack && (!pstack || !pshift))) {
            shared_st = RS_ERR_ALLOC;
        } else {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (long w = 0; w < (long)n_win; w++) {
                if (shared_st != RS_OK) continue;

                const size_t wa = (size_t)w / n_win_rg;
                const size_t wr = (size_t)w % n_win_rg;
                const size_t az0 = wa * stride_az;
                const size_t rg0 = wr * stride_rg;

                if (rs_coreg_extract(ref_img, az0, rg0, win_az, win_rg, pref) != RS_OK)
                    continue;

                if (params->estimator == RS_MICROM_EST_PHASE) {
                    /* The dominant scatterer's phase, read from each sub-look.
                     *
                     * A single pixel rather than the window's coherent sum: the
                     * sum averages scatterers that may be moving differently,
                     * which is the very thing being measured. The brightest
                     * pixel of the reference look is used throughout, so the
                     * series follows one scatterer rather than wandering to
                     * whichever happens to be brightest in each look. */
                    size_t best_i = 0;
                    double best_a = -1.0;
                    for (size_t i = 0; i < win_az * win_rg; i++) {
                        const double a = (double)cabsf(pref[i]);
                        if (a > best_a) { best_a = a; best_i = i; }
                    }
                    const size_t pa = az0 + best_i / win_rg;
                    const size_t pr = rg0 + best_i % win_rg;

                    /* PHASE AGAINST THE MEAN PHASOR, NOT ACCUMULATED FORWARD.
                     *
                     * This used to unwrap temporally: fold each step into
                     * (-pi, pi] and add it to a running total. That is the
                     * textbook way to extend the unambiguous range beyond
                     * lambda/4, and on a decorrelating aperture it is unusable,
                     * because the fold's error is added to the total and never
                     * heals.
                     *
                     * Measured on a real X-band spotlight collect: sub-look
                     * coherence is the fraction of pulses two looks share, so
                     * adjacent looks at 95 percent overlap sit at gamma = 0.85
                     * and their phase difference has a spread of 0.65 rad. Over
                     * 1548 looks the accumulated random walk is sigma*sqrt(N),
                     * about 25 radians -- four cycles of pure noise on a
                     * quantity whose whole unambiguous range is pi/2. The
                     * recovered series was a random walk in every window, on
                     * distributed clutter and on a scatterer 74x above its
                     * surroundings alike. See rs_microm_estimator_t, which
                     * carries the coherence-versus-lag measurements.
                     *
                     * Raising the overlap does not rescue it and makes it
                     * worse: per-step noise falls as sqrt(1-gamma) while the
                     * step count rises as 1/(1-gamma), and gamma levels off
                     * near 0.9 rather than reaching 1, so sigma stops falling
                     * while N keeps growing.
                     *
                     * So the accumulation is removed rather than tuned. Each
                     * look's phase is measured against the series' own mean
                     * phasor, which makes every sample independent: noise stays
                     * at sigma instead of growing as sigma*sqrt(N), and a
                     * periodic signal still averages down across the
                     * periodogram. The cost is the ambiguity the unwrap was
                     * there to buy -- motion beyond +/-lambda/4 in TOTAL now
                     * folds, where before it folded only beyond lambda/4
                     * BETWEEN looks. That is the right trade here: 1.7 mm of
                     * per-sample noise against an 8 mm range, versus a range
                     * that was unbounded in principle and swamped in practice.
                     * Use the correlation estimator for motion larger than
                     * that; it has no ambiguity at all. */
                    double amp_sum = 0.0, amp_sq = 0.0;
                    double ref_re = 0.0, ref_im = 0.0;
                    int ok = 1;

                    for (size_t k = 0; k < n_looks && ok; k++) {
                        const rs_slc_t *im = &stack->look[k];
                        if (pa >= im->n_az || pr >= im->n_rg) { ok = 0; break; }
                        const float complex z = im->data[pa * im->n_rg + pr];
                        const double a = (double)cabsf(z);
                        amp_sum += a;
                        amp_sq  += a * a;
                        /* Amplitude-weighted, so looks in which the scatterer
                         * faded contribute less to the reference they are all
                         * measured against. */
                        ref_re += (double)crealf(z);
                        ref_im += (double)cimagf(z);
                    }

                    if (!ok) goto next_window;

                    for (size_t k = 0; k < n_looks; k++) {
                        const rs_slc_t *im = &stack->look[k];
                        const float complex z = im->data[pa * im->n_rg + pr];
                        /* arg(z * conj(reference)), which is the phase relative
                         * to the mean and is already inside (-pi, pi] with no
                         * folding needed. */
                        const double pr_re = (double)crealf(z) * ref_re
                                           + (double)cimagf(z) * ref_im;
                        const double pr_im = (double)cimagf(z) * ref_re
                                           - (double)crealf(z) * ref_im;
                        const double psi = (pr_re != 0.0 || pr_im != 0.0)
                                         ? atan2(pr_im, pr_re) : 0.0;

                        const size_t idx = (size_t)w * n_looks + k;
                        out->phase[idx] = psi;
                        /* Line-of-sight displacement from phase: one wavelength
                         * of two-way path is 4*pi of phase. */
                        out->disp_los[idx] = -psi * lambda / (4.0 * M_PI);
                        out->disp_az[idx] = out->disp_rg[idx] = 0.0;
                    }

                    /* Velocity by differencing, so that the spectrum sees the
                     * same observable it does for the other estimators. */
                    for (size_t k = 0; k < n_looks; k++) {
                        const size_t idx = (size_t)w * n_looks + k;
                        if (k == 0 || stack->dt <= 0.0) { out->vel_los[idx] = 0.0; continue; }
                        out->vel_los[idx] = (out->disp_los[idx] - out->disp_los[idx - 1])
                                          / stack->dt;
                    }

                    /* Quality is amplitude STABILITY, not phase stability.
                     * Scoring a window by how constant its phase is would
                     * reward exactly the windows where nothing moves, and the
                     * gate would then select static ground. A persistent bright
                     * scatterer is the precondition for the phase meaning
                     * anything, and that is what this measures. */
                    const double mean_a = amp_sum / (double)n_looks;
                    double q = 0.0;
                    if (mean_a > 0.0) {
                        const double var = amp_sq / (double)n_looks - mean_a * mean_a;
                        const double cv = (var > 0.0) ? sqrt(var) / mean_a : 0.0;
                        q = 1.0 - cv;
                        if (q < 0.0) q = 0.0;
                        if (q > 1.0) q = 1.0;
                    }
                    out->quality[w] = q;

                    if (params->coherence_min > 0.0 && q < params->coherence_min) {
                        for (size_t k = 0; k < n_looks; k++) {
                            const size_t idx = (size_t)w * n_looks + k;
                            out->disp_los[idx] = out->vel_los[idx] = 0.0;
                            out->phase[idx] = 0.0;
                        }
                    }
                    goto next_window;
                }

                if (params->estimator == RS_MICROM_EST_SPLITBAND) {
                    /* Gather every look's patch, then estimate all shifts at
                     * once across the stack. The estimator is inherently
                     * multi-image: it cannot be decomposed into per-look calls
                     * without discarding the information that makes it better
                     * than correlation. */
                    for (size_t k = 0; k < n_looks; k++) {
                        if (rs_coreg_extract(&stack->look[k], az0, rg0, win_az, win_rg,
                                             pstack + k * win_az * win_rg) != RS_OK) {
                            goto next_window;
                        }
                    }

                    double coh = 0.0;
                    if (rs_splitband_shift(pstack, n_looks, win_az, win_rg,
                                           pshift, &coh) != RS_OK) {
                        goto next_window;
                    }

                    for (size_t k = 0; k < n_looks; k++) {
                        const size_t idx = (size_t)w * n_looks + k;
                        out->disp_az[idx] = pshift[k];
                        out->disp_rg[idx] = 0.0;
                        if (can_velocity) {
                            out->vel_los[idx] = pshift[k] * az_spacing * v_plat / r_slant;
                        }
                    }
                    out->quality[w] = coh;

                    if (params->coherence_min > 0.0 && coh < params->coherence_min) {
                        for (size_t k = 0; k < n_looks; k++) {
                            const size_t idx = (size_t)w * n_looks + k;
                            out->disp_az[idx] = out->disp_rg[idx] = 0.0;
                            out->vel_los[idx] = out->disp_los[idx] = 0.0;
                        }
                    }
                    goto next_window;
                }

                /* Summarise correlation across looks by the mean, not the
                 * minimum. With many looks the first and last share few or no
                 * pulses and correlate poorly by construction; taking the worst
                 * pair would then report near-zero quality for a window that
                 * tracks perfectly well over most of the aperture. */
                double qsum = 0.0;
                size_t qn = 0;

                /* Running absolute shift, used when accumulating differentials
                 * from adjacent looks. */
                double acc_az = 0.0, acc_rg = 0.0;

                for (size_t k = 0; k < n_looks; k++) {
                    const size_t idx = (size_t)w * n_looks + k;

                    /* PAIR measures slave(k) against master(k), so every k
                     * carries a sample -- there is no zero reference look. LAG
                     * has no sample until the lag is available. */
                    const size_t lag = (params->reference == RS_MICROM_REF_LAG)
                                     ? ((params->ref_lag > 0) ? params->ref_lag : 1)
                                     : 0;
                    if (params->reference == RS_MICROM_REF_LAG ? (k < lag)
                        : (k == 0 && params->reference != RS_MICROM_REF_PAIR)) {
                        out->disp_az[idx] = out->disp_rg[idx] = out->disp_los[idx] = 0.0;
                        continue;
                    }

                    if (params->reference == RS_MICROM_REF_LAG) {
                        /* The reference moves with the look, so it is re-extracted
                         * each time rather than held from look 0. That is the whole
                         * point of the mode: coherence is set by the lag alone. */
                        if (rs_coreg_extract(&stack->look[k - lag], az0, rg0,
                                             win_az, win_rg, pref) != RS_OK)
                            continue;
                        if (rs_coreg_extract(&stack->look[k], az0, rg0,
                                             win_az, win_rg, pcur) != RS_OK)
                            continue;
                    } else if (params->reference == RS_MICROM_REF_PAIR) {
                        if (rs_coreg_extract(&stack->look[k], az0, rg0,
                                             win_az, win_rg, pref) != RS_OK)
                            continue;
                        if (rs_coreg_extract(&stack->slave[k], az0, rg0,
                                             win_az, win_rg, pcur) != RS_OK)
                            continue;
                    } else if (rs_coreg_extract(&stack->look[k], az0, rg0,
                                                win_az, win_rg, pcur) != RS_OK) {
                        continue;
                    }

                    double sa = 0.0, sr = 0.0, pk = 0.0;
                    if (rs_coreg_shift_ex(pref, pcur, win_az, win_rg,
                                          params->upsample_az, params->upsample_rg,
                                          refine, &sa, &sr, &pk) != RS_OK) {
                        continue;
                    }

                    if (params->reference == RS_MICROM_REF_ADJACENT) {
                        /* 'sa' is the step since the previous look; the absolute
                         * shift is its running sum. The previous look becomes
                         * the reference for the next comparison. */
                        acc_az += sa;
                        acc_rg += sr;
                        sa = acc_az;
                        sr = acc_rg;
                        float complex *swap = pref; pref = pcur; pcur = swap;
                    }

                    out->disp_az[idx] = sa;
                    out->disp_rg[idx] = sr;
                    qsum += pk; qn++;

                    /* Azimuth shift to line-of-sight velocity. This observable
                     * does not wrap, which is why it is the primary one. */
                    if (can_velocity) {
                        out->vel_los[idx] = sa * az_spacing * v_plat / r_slant;
                    }

                    /* Phase refinement. The window-averaged interferometric
                     * phase between this look and the reference converts to a
                     * line-of-sight displacement, which is far finer than the
                     * tracked shift but ambiguous modulo lambda/2. Tracking
                     * above resolves the ambiguity; this supplies precision
                     * within it. */
                    double acc_re = 0.0, acc_im = 0.0;
                    const float complex *a =
                        (params->reference == RS_MICROM_REF_PAIR)
                            ? stack->look[k].data : stack->look[0].data;
                    const float complex *b =
                        (params->reference == RS_MICROM_REF_PAIR)
                            ? stack->slave[k].data : stack->look[k].data;
                    for (size_t i = 0; i < win_az; i++) {
                        const size_t row = (az0 + i) * ref_img->n_rg + rg0;
                        for (size_t j = 0; j < win_rg; j++) {
                            const float complex prod = b[row + j] * conjf(a[row + j]);
                            acc_re += (double)crealf(prod);
                            acc_im += (double)cimagf(prod);
                        }
                    }
                    out->phase[idx] = atan2(acc_im, acc_re);
                }

                /* Phase against the series' own mean phasor, for the same
                 * reason the phase estimator above no longer accumulates: a
                 * temporal unwrap adds each fold's error to a running total on
                 * an aperture whose sub-looks decorrelate, and the total is a
                 * random walk of sigma*sqrt(N) that swamps the signal. This is
                 * the interferometric phase between look and reference, which
                 * decorrelates the same way the single-pixel phase does.
                 *
                 * Averaging the phasor rather than the angle is what makes this
                 * well defined: angles near +/-pi average to nonsense, phasors
                 * do not. */
                {
                    double ref_re = 0.0, ref_im = 0.0;
                    for (size_t k = 0; k < n_looks; k++) {
                        const double p = out->phase[(size_t)w * n_looks + k];
                        ref_re += cos(p);
                        ref_im += sin(p);
                    }
                    for (size_t k = 0; k < n_looks; k++) {
                        const size_t idx = (size_t)w * n_looks + k;
                        const double c = cos(out->phase[idx]), s = sin(out->phase[idx]);
                        const double d_re = c * ref_re + s * ref_im;
                        const double d_im = s * ref_re - c * ref_im;
                        const double psi = (d_re != 0.0 || d_im != 0.0)
                                         ? atan2(d_im, d_re) : 0.0;
                        out->disp_los[idx] = -lambda / (4.0 * M_PI) * psi;
                    }
                }

                const double qmean = qn ? qsum / (double)qn : 0.0;
                out->quality[w] = qmean;

                /* Enforce the coherence mask here rather than leaving it to
                 * consumers. A window whose sub-looks do not correlate has no
                 * displacement measurement in it, only tracking noise -- and
                 * noise put through a periodogram yields a confident-looking
                 * spectral peak indistinguishable in appearance from a real
                 * structural mode. Zeroing the series makes such a window
                 * contribute a flat spectrum and no dominant frequency, while
                 * quality[] preserves why it was rejected.
                 *
                 * Note that isolated point targets on an otherwise empty scene
                 * score low here even when tracking perfectly, because most of
                 * the correlation window is background. Distributed clutter --
                 * which is what real structures sit in -- scores far higher.
                 * Set coherence_min to 0 to inspect an unmasked result. */
                if (params->coherence_min > 0.0 && qmean < params->coherence_min) {
                    for (size_t k = 0; k < n_looks; k++) {
                        const size_t idx = (size_t)w * n_looks + k;
                        out->disp_az[idx] = out->disp_rg[idx] = 0.0;
                        out->vel_los[idx] = out->disp_los[idx] = 0.0;
                    }
                }
            next_window: ;
            }
        }

        free(pref);
        free(pcur);
        free(pstack);
        free(pshift);
    }

    if (shared_st != RS_OK) {
        rs_microm_free(out);
        rs_set_error("microm: allocation failed while tracking");
        return shared_st;
    }

    if (params->remove_common_mode && out->n_win > 2) {
        rs_microm_remove_common_mode(out, az_spacing, v_plat, r_slant, can_velocity);
    }

    return RS_OK;
}


/* Largest measurable line-of-sight velocity before the shift wraps. */
double rs_microm_max_velocity(size_t win_px, double spacing_m,
                              double slant_range, double v_platform)
{
    if (win_px == 0 || spacing_m <= 0.0 || slant_range <= 0.0 || v_platform <= 0.0)
        return 0.0;
    const double max_shift_m = 0.5 * (double)win_px * spacing_m;
    return max_shift_m * v_platform / slant_range;
}


/* Sub-look count implied by the phase-ambiguity condition. */
size_t rs_microm_recommend_looks(double vib_freq, double amp_los,
                                 double t_dwell, double lambda, double overlap)
{
    if (!(vib_freq > 0.0) || !(amp_los > 0.0) || !(t_dwell > 0.0) ||
        !(lambda > 0.0) || !(overlap >= 0.0 && overlap < 1.0)) {
        return 0;
    }

    const double m_need = (4.0 * M_PI / 0.75) * vib_freq * amp_los * t_dwell / lambda;
    const double n = 1.0 + (m_need - 1.0) / (1.0 - overlap);
    if (!(n > 1.0)) return 2;
    return (size_t)ceil(n);
}
