/* Micro-motion extraction: sub-pixel offset tracking (SPOT) across a
 * sub-aperture stack, and the per-pixel vibration spectra it feeds.
 *
 * The technique is known in the literature as micro-Doppler SAR (MDSAR) and the
 * tracking step as SPOT; those names are used here so that results are
 * comparable with published work rather than privately defined.
 *
 * The realistic performance envelope, from an independent metrological
 * assessment against synchronous accelerometer ground truth on the same class
 * of X-band spotlight data this project targets (Vattulainen et al. 2026):
 * vibration frequencies 1-4 Hz, RMS radial displacement 10.43 mm down to
 * 0.10 mm, radial velocity error of order 1 mm/s, frequency resolution 0.06 Hz
 * from a 16 s acquisition.
 *
 * THAT BAND IS THE TARGET'S, NOT THE TECHNIQUE'S. It was measured on bridges,
 * which vibrate at a few hertz. The same group, on the same class of Umbra
 * X-band data, recovers 87 Hz from an idling van and 36 Hz from a ship's
 * engine, both confirmed against triaxial accelerometers (Clemente et al.,
 * EuRAD 2025). Quoting 1-4 Hz as a property of the method would understate its
 * reach by more than an order of magnitude.
 *
 * What sets the ceiling is the sub-aperture sampling rate, and reaching those
 * frequencies takes far more sub-looks than seems natural: that work uses
 * sub-apertures of 0.4167 s stepped by about 4 ms -- roughly 99 percent overlap
 * and some 3900 sub-apertures across a 16.4 s dwell, sampling near 240 Hz.
 * Configurations here have typically used tens of looks, which caps the
 * observable band in the single hertz and cannot represent an engine line at
 * all. The overlap is deliberate and is justified there exactly as in
 * time-frequency analysis: it restores time resolution without giving up
 * azimuth resolution.
 *
 * Note also that the observable differs. That work reads the PHASE of each
 * pixel in each sub-aperture directly; this interface defaults to
 * correlation-based offset tracking. See rs_microm_estimator_t.
 *
 * One limitation from that work shapes this interface: frequencies are
 * recovered reliably while relative amplitudes are not. The amplitude field
 * below is therefore a qualitative indicator, and callers must label it as
 * such wherever it is presented. */

#ifndef RESONARSAT_MICROM_H
#define RESONARSAT_MICROM_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/subaperture.h"

/* Tracking parameters.
 *
 * The defaults come from published working values rather than round numbers:
 * azimuth upsampling 10, range upsampling 20, and a patch tuned per target in
 * the range 51 to 131 pixels. Patch size is not a free parameter -- too small a
 * patch makes the tracker systematically underestimate displacement -- so it is
 * exposed prominently and rs_microm_params_default() picks a mid-range value
 * rather than pretending there is a universally correct one. */
/* Which sub-look each correlation is taken against.
 *
 *   RS_MICROM_REF_FIRST     Every look against look 0. Simple, and the shifts
 *                           come out absolute with no accumulation. But the
 *                           correlation surface of an N-pixel window is only
 *                           unambiguous over +/-N/2 pixels, and a target's
 *                           excursion relative to look 0 spans TWICE its
 *                           amplitude. Exceed that and the shift wraps, folding
 *                           the recovered series and putting its energy at twice
 *                           the true frequency -- a clean, confident-looking
 *                           second harmonic that is entirely an artefact.
 *
 *   RS_MICROM_REF_ADJACENT  Each look against its predecessor, accumulating the
 *                           differentials. Consecutive looks are separated by
 *                           one sampling interval, so the shift between them is
 *                           smaller than the full excursion by a factor of
 *                           2*pi*f*dt -- comfortably inside the unambiguous
 *                           range where the direct comparison is not. They also
 *                           share most of their pulses, so they correlate far
 *                           better than look 0 and look N do.
 *
 * MEASURED TRADE-OFF, and the reason FIRST remains the default. ADJACENT does
 * what it promises -- per-window coherence roughly doubles and the second
 * harmonic disappears from the target's own windows. But accumulating the
 * differentials integrates tracking noise into a random walk, whose 1/f^2
 * spectrum swamps the signal: on the single-target frequency sweep ADJACENT
 * takes recovery from 5 of 6 to 0 of 6, every window reporting the lowest bin.
 * The least-squares detrend in rs_spectrum_compute() removes only the linear
 * part of a random walk, which is not enough. Using ADJACENT profitably needs
 * the accumulation handled -- differencing the series back, or a high-pass --
 * and that is not yet implemented. PAIR below is the principled way to get the
 * same benefit, because it never accumulates in the first place.
 *
 *   RS_MICROM_REF_PAIR      Each look's slave against its own master, the pair
 *                           held B_shift apart. This is what the sources
 *                           describe: WO2024008365A1 [0004] sweeps two bands
 *                           "rigidly held at a distance B_shift" and tracks the
 *                           slave's position against its master, and the Giza
 *                           paper's Eqs. 4-5 give N_D masters and N_D slaves.
 *                           Needs a stack built with rs_subap_params_t.pair.
 *
 *                           Each sample is a displacement difference across one
 *                           FIXED lag dt = B_shift * t_dwell / B_CD, so the
 *                           series is a first difference of the displacement
 *                           rather than ADJACENT's running sum of differences.
 *                           Differencing is a high-pass with a known response,
 *                           |2 sin(pi f dt)|; summing is an integrator that
 *                           turns tracking noise into a random walk. So the
 *                           expectation was that PAIR would recover frequencies
 *                           where ADJACENT cannot.
 *
 *                           IT DOES NOT, AND THE REASON IS NOW UNDERSTOOD. On
 *                           the same single-target fixture that ADJACENT fails,
 *                           PAIR fails too -- and on the distributed-texture
 *                           fixture its series is exactly zero in every window
 *                           that holds the target: at the default B_shift the
 *                           slave-master offset never moves one quantisation
 *                           step, and the mode's occasional nonzero windows are
 *                           near-empty edge windows whose blips do not depend
 *                           on the injection: measured across five seeds and
 *                           six injected frequencies, on and off the fixture's
 *                           blind bins, the mode's answers pool to {0.1 Hz,
 *                           0.5 Hz} whatever is injected.
 *
 *                           Two structural reasons, not one defect. First, the
 *                           paper decomposition's sweep spans the master band's
 *                           own width, so the record length always equals the
 *                           sub-look duration (N*dt = t_sap) and every
 *                           resolvable frequency bin sits at an INTEGER
 *                           observation ratio -- each sub-look integrates a
 *                           whole number of cycles of any frequency the sweep
 *                           can resolve. A displacement-averaging observable
 *                           (which correlation tracking is) cannot both resolve
 *                           a frequency and retain it, at any left_out_frac.
 *                           Second, the pair's lag is one sweep step, so its
 *                           differential is attenuated by |2 sin(pi f dt)| on
 *                           top of that -- at 96 looks, under one quantisation
 *                           step for everything in band.
 *
 *                           The selection-artefact history is worth keeping:
 *                           PAIR used to return the lowest spectral bin for
 *                           every injection, which turned out to be
 *                           rs_spectrum_best_window() choosing windows whose
 *                           whole excursion was one sub-pixel step. The
 *                           quantisation floor removed that, and the remaining
 *                           "recoveries" (2 of 9 at 0.5 Hz) were then shown to
 *                           be coincidence: across five seeds the mode answers
 *                           0.5 Hz at a fixed per-seed rate whatever frequency
 *                           is injected, including with coherent whole-patch
 *                           motion injected at 0.35 Hz.
 *
 *                           AN EARLIER SUSPECT WAS CHECKED AND WAS NOT THE
 *                           CAUSE. rs_microm_ref_t once named a systematic
 *                           drift of the master-slave offset with sweep
 *                           position; the band layout did carry such a defect
 *                           (the sweep began half a step above the lower band
 *                           edge, clipping the last slave's filter), it was
 *                           fixed in rs_subaperture_split(), and the measured
 *                           result did not change. Real defect, wrong suspect.
 *
 *                           So this mode is FAITHFUL TO THE SOURCES AND NOT FIT
 *                           FOR MEASUREMENT. It is implemented and exposed
 *                           because the sources describe it and it should be
 *                           testable, not because it works. Do not read a
 *                           frequency out of it. test_tracking.c records the
 *                           behaviour so a fix would show up as that test
 *                           changing.
 *
 *                           If it is fixed, amplitudes will still be attenuated
 *                           by |2 sin(pi f dt)| and not comparable to FIRST's.
 *
 *   RS_MICROM_REF_LAG       Each look against the one 'ref_lag' places before
 *                           it, with NO accumulation. Untested; added to
 *                           separate the two defects the other three carry.
 *
 *                           FIRST decorrelates: looks i and j share pulses only
 *                           while |i-j| < 1/(1-overlap), so with a fixed
 *                           reference the coherent span is a handful of looks
 *                           out of hundreds. Measured on the single-target
 *                           fixture, the correlation peak of look 0 against the
 *                           rest averages 0.090 at zero overlap and only 0.310
 *                           at 0.90, while ADJACENT PAIRS at 0.90 correlate at
 *                           0.913. The coherence exists; a fixed reference
 *                           discards it.
 *
 *                           ADJACENT keeps that coherence and then integrates
 *                           tracking noise into a random walk. PAIR avoids the
 *                           integrator but inherits the published sweep's
 *                           geometry, where N*dt = t_sap puts every resolvable
 *                           bin at an integer observation ratio -- a null of the
 *                           averaging response.
 *
 *                           LAG is the remaining corner: a fixed short lag on
 *                           the PULSE route. Coherence is set by the lag alone
 *                           and stays high for small 'ref_lag'; nothing
 *                           accumulates, so there is no random walk; and the
 *                           record length is the whole dwell rather than one
 *                           sub-look, so df = 1/T and the bins do NOT land on
 *                           integer observation ratios the way the spectral
 *                           sweep forces them to.
 *
 *                           Each sample is a displacement difference across
 *                           ref_lag*dt, so the series is a first difference with
 *                           response |2 sin(pi f ref_lag dt)|. That nulls at
 *                           f = k/(ref_lag*dt); keep ref_lag small enough that
 *                           the first null sits above the band of interest.
 *                           Frequencies survive differencing, amplitudes are
 *                           attenuated by that response and are NOT comparable
 *                           to FIRST's.
 */
typedef enum {
    RS_MICROM_REF_FIRST = 0,
    RS_MICROM_REF_ADJACENT = 1,
    RS_MICROM_REF_PAIR = 2,
    RS_MICROM_REF_LAG = 3
} rs_microm_ref_t;

/* Which estimator computes the shifts.
 *
 *   RS_MICROM_EST_CORRELATION  Cross-correlation peak finding, refined to
 *                              sub-pixel by local upsampling. The classical
 *                              method, and characterised in the literature as
 *                              achieving "performances in the order of the
 *                              resolution element for a few independent
 *                              samples" -- which is what this project measures.
 *
 *   RS_MICROM_EST_SPLITBAND    Split-band Phase Linking across the whole stack
 *                              (see phaselink.h). Uses all N^2 interferograms
 *                              rather than the N-1 formed against one reference,
 *                              and comes within 0.5 dB of the Cramer-Rao bound.
 *                              Requires interferometric coherence between looks;
 *                              it has nothing to track if that is absent.
 *
 *   RS_MICROM_EST_PHASE        The phase of a single dominant pixel, read
 *                              directly from each sub-look and unwrapped across
 *                              them. This is the observable of Clemente et al.
 *                              (EuRAD 2025), the only one in this file with
 *                              published accelerometer validation on this class
 *                              of data.
 *
 *                              IT IS NOT A REFINEMENT OF THE OTHER TWO -- it
 *                              measures something different. Both of those
 *                              estimate WHERE a patch sits, averaged over the
 *                              sub-look, so they live in the displacement-
 *                              averaging regime and go blind at frequencies
 *                              near k/t_sap. Phase responds to sub-wavelength
 *                              line-of-sight motion without that averaging, and
 *                              the same work recovers 36 Hz at an observation
 *                              ratio of exactly 18, where the averaging model
 *                              predicts precisely zero. See
 *                              rs_spectrum_subaperture_response().
 *
 *                              The cost is ambiguity. Phase wraps every
 *                              lambda/2 of line-of-sight motion -- about 16 mm
 *                              at X band -- so a target moving further than
 *                              that between sub-looks unwraps wrongly and the
 *                              recovered series is nonsense. Correlation has no
 *                              such limit. Prefer phase for small fast motion
 *                              and correlation for large slow motion.
 *
 *                              SUB-LOOK COHERENCE IS SET BY PULSE SHARING, AND
 *                              THAT BOUNDS HOW LONG A SERIES CAN BE UNWRAPPED.
 *                              Measured on a real X-band spotlight collect, the
 *                              coherence between two sub-looks is very nearly
 *                              the fraction of pulses they have in common:
 *
 *                                shared  95%  90%  75%  50%   0%
 *                                gamma  0.85 0.78 0.61 0.39 0.07
 *
 *                              It reaches the noise floor exactly when the
 *                              windows stop overlapping. That is independent
 *                              speckle per pulse subset -- the coherence is the
 *                              shared fraction of the energy -- and it was the
 *                              same on distributed clutter and on a scatterer
 *                              74x above its surroundings, so it is not an SNR
 *                              effect and not a property of the target.
 *
 *                              THE CONSEQUENCE IS COUNTER-INTUITIVE: RAISING
 *                              THE OVERLAP MAKES THE UNWRAP WORSE. Per-step
 *                              phase noise falls as the overlap rises, but the
 *                              number of steps rises faster, and the unwrap
 *                              accumulates them as a random walk of
 *                              sigma*sqrt(N). Worse, gamma does not approach 1
 *                              as the overlap does -- it levels off near 0.9 --
 *                              so sigma stops falling while N keeps growing.
 *                              Over a 33 s dwell the accumulated error is tens
 *                              of radians at every overlap tried, where a
 *                              usable unwrap needs it well under pi.
 *
 *                              So do NOT reach for 99 percent overlap to
 *                              stabilise a phase series; it does the opposite.
 *                              Use the FEWEST looks that still sample the
 *                              frequency of interest, and treat any unwrapped
 *                              phase series spanning a full aperture as suspect
 *                              until its accumulated noise is checked. A
 *                              peak-to-peak line-of-sight velocity near
 *                              lambda/(2*dt) is the tell: that is the ceiling
 *                              the fold imposes, and hitting it means the
 *                              series wrapped rather than that the target moved.
 *
 *                              THE SHUFFLE NULL TEST IS INVALID FOR THIS
 *                              ESTIMATOR. rs_null_floor(), --shuffle-looks and
 *                              --null-trials destroy the sub-look time order and
 *                              hold everything else constant -- for a
 *                              correlation observable. For phase they do not:
 *                              reordering puts non-consecutive looks adjacent,
 *                              which is precisely where a phase series steps
 *                              furthest, so the shuffle inflates the per-step
 *                              noise it is supposed to preserve. Measured on the
 *                              Giza collect at 128 looks and 0.99 overlap, the
 *                              median largest step is 0.052 rad in order and
 *                              1.878 rad shuffled -- a factor of 36. A drifting
 *                              series therefore beats its own shuffles by
 *                              construction: that run cleared 32 of them at
 *                              p = 0.03 while every window on the grid reported
 *                              whatever bin the search was told to start at, and
 *                              8 simulated motionless collects reproduced the
 *                              same frequency at 99 percent of the same
 *                              prominence. Adjudicate a phase result with
 *                              rs_null_static() (--null-static), which a
 *                              motionless scene cannot walk over because it
 *                              carries the same overlap, unwrap and detrend.
 *                              See runs/giza/2026-07-30-uniform-phase-khufu/.
 */
typedef enum {
    RS_MICROM_EST_CORRELATION = 0,
    RS_MICROM_EST_SPLITBAND = 1,
    RS_MICROM_EST_PHASE = 2
} rs_microm_estimator_t;

typedef struct {
    rs_microm_estimator_t estimator;  /* how shifts are computed */
    rs_microm_ref_t reference;  /* which correlation reference (correlation only) */
    /* Lag in looks for RS_MICROM_REF_LAG, ignored by every other mode.
     *
     * Small keeps coherence high -- looks share pulses while the lag is under
     * 1/(1-overlap) -- and pushes the differencing null at 1/(ref_lag*dt) up
     * out of the band. Both want it small; the only thing wanting it large is
     * sensitivity, since |2 sin(pi f ref_lag dt)| grows with the lag until that
     * null. Default 1, which is ADJACENT's spacing without ADJACENT's
     * integrator. */
    size_t ref_lag;
    size_t win_az, win_rg;      /* correlation patch size in pixels */
    size_t stride_az, stride_rg;/* step between patch centres, pixels */
    size_t upsample_az;         /* sub-pixel refinement factor along azimuth */
    size_t upsample_rg;         /* sub-pixel refinement factor along range */
    double coherence_min;       /* discard windows below this correlation peak */

    /* Subtract the scene-median shift from every window, look by look.
     *
     * Sub-looks are focused from different pulse windows, so they can carry a
     * systematic offset relative to one another -- residual co-registration
     * rather than target motion. That offset is common to every window, so it
     * appears at the same frequency across the whole scene and can outrank a
     * genuine localised target: on a four-target test scene it produced a
     * prominent peak in static-ground windows that won window selection at every
     * injected frequency, giving the same wrong answer each time.
     *
     * The median is used rather than the mean so that a few genuinely moving
     * windows do not drag the estimate. Real localised motion survives; anything
     * shared by the whole scene cancels.
     *
     * MEASURED: off by default. The median is only as good as the windows it is
     * taken over, and on a scene dominated by empty background -- every
     * synthetic fixture here -- most windows track noise, so the median is noisy
     * and subtracting it injects noise into the one window that had signal:
     * recovery falls from 5 of 6 to 3 of 6. Enable on real scenes with
     * distributed clutter, where the premise that most of the scene is static
     * and well-tracked actually holds. */
    int remove_common_mode;

    /* Nonzero: run the tracker as an unoptimised reference. Set by --no-optimize.
     *
     * Two things change, and only one of them can change a number:
     *
     *   1. The correlator searches the WHOLE zero-padded surface for its peak
     *      instead of the neighbourhood of the strongest integer sample
     *      (RS_COREG_REFINE_EXHAUSTIVE). This CAN move a reported shift, and is
     *      the point of the flag. See rs_coreg_refine_t for the precise
     *      circumstance under which the two disagree -- it is narrower than it
     *      sounds, because the optimised path's integer peak is already a global
     *      maximum over the sampled surface.
     *
     *   2. The loop over windows runs on one thread. This cannot change a number.
     *      Windows are independent and each writes only its own output slots.
     *
     * MEASURED COST, WHICH IS FAR SMALLER THAN THE MECHANISM SUGGESTS. The
     * exhaustive correlator is 1.7x to 3.2x the optimised one per call across
     * every configuration tried: 2.5x at the tomo and mmotion defaults (32x32
     * window, 10x20 upsampling, a 320x640 padded surface at 1.56 MB), 1.7x at
     * 24x24 and 10x10, 3.2x at 64x64 and 10x20.
     *
     * That is not the ratio one guesses from "upsample the whole surface instead
     * of a neighbourhood", and the reason is worth knowing: the optimised path is
     * not cheap either. It evaluates (2*upsample_az+1)*(2*upsample_rg+1) points
     * and each one costs O(win_az*win_rg), so at the defaults it already does
     * about 880 thousand complex multiply-accumulates -- against roughly 3.6
     * million butterflies for the padded transform. Same order. The O(N) per
     * refinement point is what closes the gap.
     *
     * Serialising the window loop costs more than the search does: about 4x
     * wall-clock on eight cores. Both together put a full-scale run in single
     * digits, not orders of magnitude.
     *
     * A run made with this set is NOT a better measurement than one without. It
     * is a second measurement by a slower route, whose only use is comparison
     * with the first. Neither passes a null test on its own -- see README.md. */
    int no_optimize;
} rs_microm_params_t;

/* Per-window micro-motion result over the whole sub-look stack.
 *
 * 'n_win_az' by 'n_win_rg' windows, each carrying a displacement time series of
 * 'n_looks' samples. 'disp_az' and 'disp_rg' hold the tracked shift in pixels,
 * 'disp_los' the line-of-sight displacement in metres derived from phase, and
 * 'quality' the correlation peak value used as a mask.
 *
 * Series are stored window-major: window w's series begins at index
 * w * n_looks. */
typedef struct {
    double *disp_az, *disp_rg;  /* [n_win][n_looks], pixels */

    /* Line-of-sight VELOCITY, m/s, derived from the tracked azimuth shift.
     *
     * This is the primary observable, and the reason is worth stating. A target
     * moving radially is displaced in azimuth by dx = R * v_r / V -- the classic
     * moving-target azimuth shift -- so the tracked shift measures a velocity,
     * not a displacement. Being a geometric shift rather than a phase, it is
     * unambiguous: nothing here wraps, however large the motion, so long as the
     * shift stays inside the correlation window.
     *
     * This is also the quantity the independently validated literature reports
     * (velocity errors of order 1 mm/s against accelerometer ground truth), so
     * results computed from it are directly comparable with published work.
     *
     * A sinusoidal displacement at frequency f produces a sinusoidal velocity at
     * the same f, so frequency recovery is unaffected by working in velocity. */
    double *vel_los;            /* [n_win][n_looks], m/s */

    /* Line-of-sight DISPLACEMENT, m, from the phase taken against the series'
     * own mean phasor. Far finer than the velocity estimate, and ambiguous
     * beyond +/-lambda/4 of TOTAL motion -- about 8 mm at X band -- because it
     * is deliberately not unwrapped. See rs_microm_track() for why unwrapping
     * was removed rather than fixed. */
    double *disp_los;           /* [n_win][n_looks], metres */

    /* Phase relative to the series' mean phasor, rad, kept so a caller can
     * choose a different reference or attempt its own unwrap. Note this is no
     * longer an accumulated total: consecutive values are independent. */
    double *phase;              /* [n_win][n_looks], radians in (-pi, pi] */

    double *quality;            /* [n_win], mean correlation peak in [0,1] */

    size_t n_win_az, n_win_rg, n_win, n_looks;
    size_t win_az, win_rg, stride_az, stride_rg;

    double dt;      /* s between consecutive looks, from the stack */
    double f_max;   /* Hz, vibration Nyquist limit */

    /* Pixel spacings carried through from the sub-look geometry, so that
     * consumers can convert the tracked shifts from pixels to metres without
     * needing the image stack. */
    double az_spacing_m, rg_spacing_m;

    /* Sub-pixel quantisation of the tracked shift, in pixels: 1/upsample_az for
     * the correlation estimator, and zero when the concept does not apply.
     *
     * This is a FLOOR on what the tracker can report, not a precision figure,
     * and it has to travel with the result because a consumer cannot otherwise
     * tell a measurement from a rounding artefact. A window whose whole
     * excursion is one step returns a two-valued series, and a two-valued
     * series has its energy at low frequency whatever produced the transitions.
     * See rs_spectrum_best_window(), which refuses to select such a window.
     *
     * Zero for RS_MICROM_EST_PHASE, whose observable is pixel phase rather than
     * a correlation offset and whose limit is set by phase noise instead. */
    double quant_px;
} rs_microm_t;

/* Per-window vibration spectrum and the summary maps derived from it. */
typedef struct {
    double *psd;            /* [n_win][n_freq] power spectral density */
    double *freq;           /* [n_freq] frequency axis, Hz */
    double *dominant_freq;  /* [n_win] frequency of the largest peak, Hz */
    double *amplitude;      /* [n_win] QUALITATIVE peak strength, not calibrated */
    double *quality;        /* [n_win] copied from the tracking mask */

    /* [n_win] peak-to-peak excursion of the tracked azimuth shift, in PIXELS,
     * before detrending -- and 'quant_px' carried through from rs_microm_t.
     *
     * These two travel together because neither means anything alone. The
     * excursion says how far the tracker saw the patch move; the quantisation
     * says how far it could have seen it move by accident. Their ratio is what
     * decides whether a window holds a measurement, and it is the only quantity
     * here that is not scale-free -- prominence is a ratio of powers and
     * therefore identical for a strong peak and for rounding noise. */
    double *excursion_px;
    double  quant_px;

    /* [n_win] spectral prominence: the dominant peak's power divided by the mean
     * power of the rest of the spectrum.
     *
     * This is the metric that identifies which windows actually contain a
     * vibrating target, and it earns its place by being much better at it than
     * the obvious alternatives. Selecting the window with the largest
     * displacement excursion picks the NOISIEST window, because noise excursions
     * exceed real ones; selecting by tracking coherence picks the window that
     * correlates best, which is usually static ground. Prominence asks the
     * relevant question -- does this window's motion concentrate at one
     * frequency -- and on a synthetic sweep it recovers five injected
     * frequencies of six where excursion-based selection recovers two. */
    double *prominence;

    size_t n_win, n_win_az, n_win_rg, n_freq;
    double df;              /* Hz per spectral bin */
} rs_spectrum_t;

/* Fill 'params' with the published working defaults described above. */
void rs_microm_params_default(rs_microm_params_t *params);

/* Track every window across the sub-look stack and extract displacement series.
 *
 * Sub-look 0 is the reference. For each window and each subsequent look, a
 * two-dimensional normalised cross-correlation against the reference gives the
 * apparent shift, refined to sub-pixel precision by local upsampling of the
 * correlation surface around its peak (the Guizar-Sicairos approach: refine by
 * evaluating an upsampled inverse transform in a small neighbourhood, rather
 * than upsampling the whole surface).
 *
 * The tracked azimuth shift converts to a line-of-sight velocity through
 * v_r = dx * V / R, with dx the shift in metres, V the platform speed and R the
 * slant range -- all measured from the collection rather than assumed. This is
 * the unambiguous observable and the one to prefer.
 *
 * Line-of-sight displacement is separately obtained from the interferometric
 * phase between look and reference, averaged over the window, as
 * d = -lambda/(4*pi) * phi. Phase is far more precise than tracking but wraps
 * modulo lambda/2 -- 15.6 mm at X-band, which real structural motion exceeds.
 *
 * THE PHASE SERIES IS NOT UNWRAPPED. Each look's phase is expressed relative to
 * the series' own mean PHASOR -- psi = arg(z * conj(mean z)) -- which lands in
 * (-pi, pi] with nothing to fold and, crucially, makes every sample independent
 * of every other.
 *
 * It used to unwrap temporally, accumulating folded differences, which is the
 * textbook way to trade total range for per-step range. That is unusable on an
 * aperture whose sub-looks decorrelate, and the failure is not marginal. Sub-look
 * coherence is the fraction of pulses two looks share, so even 95 percent
 * overlap gives gamma = 0.85 and a per-step phase spread of 0.65 rad; over 1548
 * looks the accumulated random walk is sigma*sqrt(N) ~ 25 rad, against an
 * unambiguous range of pi/2. Measured on real data, the recovered series was a
 * random walk in every window, on distributed clutter and on a scatterer 74x
 * above its surroundings alike. Raising the overlap makes it worse, not better.
 * See rs_microm_estimator_t, which carries the measured coherence-versus-lag
 * numbers.
 *
 * WHAT THIS COSTS. The ambiguity the unwrap bought is gone: motion beyond
 * +/-lambda/4 in TOTAL now folds, where before it folded only beyond lambda/4
 * BETWEEN looks. At X band that is about 8 mm of line-of-sight motion. Against
 * it, per-sample noise is now bounded at sigma rather than growing without
 * limit, and a periodic signal still averages down across the periodogram.
 *
 * So prefer 'vel_los' -- the tracked shift -- for motion larger than lambda/4;
 * it has no ambiguity at any amplitude. Use 'disp_los' for small motion, where
 * it is far finer. 'phase' now holds psi rather than an accumulated total, so a
 * caller wanting a different reference or its own unwrap has the raw quantity.
 *
 * The loop over windows is parallelised with OpenMP when available. On success
 * '*out' owns its arrays and must be released with rs_microm_free().
 *
 * Returns RS_ERR_ARG if the stack has fewer than two looks, or if the window
 * size exceeds the image. */
resonarsat_status_t rs_microm_track(const rs_subap_stack_t *stack,
                                    const rs_microm_params_t *params,
                                    rs_microm_t *out);

/* Release everything a micro-motion result owns. */
void rs_microm_free(rs_microm_t *m);

/* Which observable a spectrum is computed from.
 *
 * Default to velocity: it does not wrap, and it is what the validated
 * literature reports. Displacement is finer but inherits the unwrapping caveat
 * on rs_microm_track(). Both give the same peak FREQUENCY for a sinusoid. */
typedef enum {
    RS_SPEC_VELOCITY = 0,
    RS_SPEC_DISPLACEMENT = 1
} rs_spectrum_source_t;

/* Compute per-window vibration spectra from tracked displacement series.
 *
 * A Hann-windowed periodogram of the selected observable's series. The
 * series are short -- one sample per sub-look -- so this is cheap, and a plain
 * periodogram is the honest choice at this record length; high-resolution
 * sparse estimators can resolve closely spaced peaks but invite reading
 * structure into noise when only a few dozen samples support them.
 *
 * The mean is removed before transforming, so a static offset does not appear
 * as a spurious zero-frequency peak, and the zero bin is excluded when the
 * dominant frequency is selected.
 *
 * On success '*out' owns its arrays and must be released with
 * rs_spectrum_free(). */
resonarsat_status_t rs_spectrum_compute(const rs_microm_t *m,
                                       rs_spectrum_source_t source,
                                       rs_spectrum_t *out);

/* As rs_spectrum_compute(), but ignoring every bin below 'f_min' hertz when
 * choosing the dominant peak and when computing prominence.
 *
 * The lowest bins are where a drift lives, and a drift is not a vibration. The
 * linear detrend applied to every series removes the straight-line part of one,
 * but a random walk or any curved trend survives it and piles its energy into
 * the first two or three bins -- see rs_spectrum_best_window() on how that can
 * outscore a genuine target. The quality gate normally keeps such windows out,
 * but the gate is relative, so it stops protecting when every window in the
 * scene tracks poorly.
 *
 * A band floor is the direct remedy and it is also a test. If a detection
 * survives excluding the lowest bins, whatever it is has structure at a real
 * frequency; if it evaporates, it was a trend. Passing 0 gives the unrestricted
 * behaviour, which is what rs_spectrum_compute() does. */
resonarsat_status_t rs_spectrum_compute_band(const rs_microm_t *m,
                                            rs_spectrum_source_t source,
                                            double f_min,
                                            rs_spectrum_t *out);

/* What is removed from each displacement series before its spectrum is taken.
 *
 *   RS_DETREND_LINEAR  Least-squares straight line. The default, and what the
 *                      two calls above use.
 *   RS_DETREND_MEAN    Mean only, leaving any ramp in the record.
 *   RS_DETREND_NONE    Nothing removed.
 *
 * WHY THIS IS SELECTABLE RATHER THAN FIXED. Detrending is not in the source
 * material; it was added here because temporal phase unwrapping random-walks
 * when the phase is noisy, producing a ramp that puts all the energy in the
 * lowest bin and makes every window report the same spurious "dominant
 * frequency" of one bin width.
 *
 * But it is not a neutral cleaning step, and under one of the models it is
 * actively destructive. A resonance interpretation maps depth as z = v/(2f), so
 * the LOWEST frequencies are the DEEPEST structure. Removing a straight line
 * preferentially attenuates exactly the part of the record that model reads as
 * deep, which means a linear detrend could suppress a real deep signal rather
 * than only a spurious one. Under Model A depth is linear in the transform's
 * bin index and no such asymmetry arises.
 *
 * A choice that helps one model and hurts another must not be a hidden default.
 * Turning it off is also a test in the same spirit as the band floor: if a
 * detection appears only without detrending and sits in the first bins, it is a
 * trend; if it survives both settings, it is not. */
typedef enum {
    RS_DETREND_LINEAR = 0,
    RS_DETREND_MEAN = 1,
    RS_DETREND_NONE = 2
} rs_detrend_t;

/* As rs_spectrum_compute_band(), with the detrend explicit. The two calls above
 * are wrappers passing RS_DETREND_LINEAR. */
resonarsat_status_t rs_spectrum_compute_opts(const rs_microm_t *m,
                                             rs_spectrum_source_t source,
                                             double f_min,
                                             rs_detrend_t detrend,
                                             rs_spectrum_t *out);

/* Release everything a spectrum result owns. */
void rs_spectrum_free(rs_spectrum_t *s);

/* Find the window whose spectrum shows the most prominent peak, among those
 * that resolved any motion at all, writing its index to '*out_window' and its
 * prominence to '*out_prominence' when those are non-NULL.
 *
 * This is how a caller should pick a window to report from a scene when it does
 * not already know where the target is. See rs_spectrum_t.prominence for why
 * the obvious alternatives are worse.
 *
 * THE QUANTISATION FLOOR, AND WHY PROMINENCE ALONE CANNOT SUPPLY IT.
 * Prominence is peak power over mean power -- a ratio, and therefore
 * scale-free. Multiply a window's whole series by any constant and its
 * prominence does not move. So prominence cannot distinguish a strong peak from
 * rounding noise, and the tracker's output is quantised: shifts are located to
 * 1/upsample of a pixel, and a window whose entire excursion is one such step
 * returns a two-valued series whose energy sits at low frequency whatever
 * produced the transitions. Such a window routinely out-scores real ones. It
 * also passes any coherence gate at the top, because a patch that never moves
 * correlates with itself perfectly.
 *
 * A window is therefore a candidate only if its excursion clears the
 * quantisation noise by three sigma. Quantisation error is uniform on +/-q/2
 * with RMS q/sqrt(12); a sinusoid of peak-to-peak A has RMS A/(2*sqrt(2)); so
 * the ratio is 1.2247*A/q and three sigma is A >= 2.449*q. That is a derived
 * limit rather than a tuned threshold, which matters because the alternative --
 * picking a fraction that works on the scene in front of you -- is how a
 * selection rule ends up fitted to one fixture.
 *
 * WHAT THIS DOES NOT DO. It does not make a marginal observable measurable.
 * On the distributed-texture fixture the floor moved RS_MICROM_REF_PAIR's
 * count from 0/9 to an apparent 2/9 -- and wider reproduction (five seeds,
 * frequencies on and off the fixture's blind bins) then showed those two, and
 * RS_MICROM_REF_FIRST's one, to be artefacts whose answers do not depend on
 * the injection at all: honest zero everywhere. The value is in what it
 * refuses, not in what it finds: without it the
 * answer is a confident wrong frequency, with it the answer is either honestly
 * absent or an artefact that repetition across seeds exposes.
 *
 * Returns RS_ERR_ARG if the spectrum is empty, and RS_ERR_RANGE if no window
 * cleared the floor -- which is a result, not a failure, and callers must
 * report it as "nothing resolved" rather than falling back to a window. When
 * 'quant_px' is zero the floor cannot be evaluated and every window is a
 * candidate, which is the behaviour for estimators whose limit is not a
 * correlation quantisation.
 *
 * HOW MANY WINDOWS WERE ELIGIBLE IS PART OF THE ANSWER, which is why
 * 'out_n_candidates' exists and why callers should present it beside the
 * selection rather than treating it as a diagnostic.
 *
 * It counts windows surviving BOTH gates -- the relative coherence gate and the
 * floor -- so a caller must not describe it as the floor's doing alone. On a
 * phase run, where 'quant_px' is zero and the floor never runs, every exclusion
 * is the coherence gate.
 *
 * The floor is three sigma for ONE window. It is then applied to every window
 * independently and the best survivor returned, with nothing accounting for how
 * many were tried -- so the chance that some window crosses it on quantisation
 * noise alone grows with the grid. Measured on real data: the same scene and
 * the same chain gave 225 of 225 windows at exactly zero excursion and an
 * honest RS_ERR_RANGE at one grid size, and 958 of 961 at zero with two
 * crossings and a confident-looking "0.183 Hz, prominence 29.9" at a larger
 * one. The observable was identical; only the number of opportunities changed.
 *
 * So a small count is the signature of that effect, and there is a
 * non-arbitrary place to draw the line. Windows are laid down at a stride,
 * typically half their width, so they OVERLAP: a target big enough to be
 * resolved at all falls inside a 2x2 block of them at minimum. Fewer than four
 * qualifying windows therefore cannot describe a spatially resolved mode, and
 * two scattered ones describe noise. That bound comes from the window geometry
 * rather than from a tuned constant, which matters here for the same reason it
 * mattered for the floor itself.
 *
 * STILL OPEN, and measured: on one real scene the same chain gave an honest
 * RS_ERR_RANGE at 225 windows and "0.183 Hz, prominence 29.9" at 961, off two
 * chance crossings, with the observable identically zero in 99.7 percent of
 * windows both times. A proper fix -- a multiplicity correction over the
 * EFFECTIVE number of independent windows, which is fewer than n_win because
 * overlapping windows are correlated, or ranking on the qualifying FRACTION
 * and its spatial contiguity instead of on one window's prominence -- is not
 * implemented.
 * Reporting the count is the honest interim: "2 of 961 cleared the floor" reads
 * as noise where "prominence 29.9" reads as a detection, and they are the same
 * result. Pass NULL if the count genuinely is not wanted. */
resonarsat_status_t rs_spectrum_best_window(const rs_spectrum_t *spec,
                                            size_t *out_window,
                                            double *out_prominence,
                                            size_t *out_n_candidates);

/* The frequency the most windows agree on, and how many agree.
 *
 * rs_spectrum_best_window() answers "which single window is most prominent";
 * this answers "what do the windows agree on", which is a different question and
 * the one a detection needs. A window can clear every gate and still report a
 * frequency no other window reports -- measured on a 3.000 Hz injection, 23 of
 * 49 windows made 2.995 Hz their top bin while 16 made 2.604 Hz theirs, and the
 * single most prominent window was one of the sixteen.
 *
 * Windows agree when their dominant frequencies fall in the same bin, tested at
 * half a bin. Only windows passing the SAME gates as rs_spectrum_best_window()
 * vote, so the two functions describe the same population and their counts are
 * comparable.
 *
 * 'out_n_agree' and 'out_n_distinct' ARE THE ANSWER, not diagnostics. The
 * measured behaviour on synthetic fixtures with known ground truth:
 *
 *     agreement      distinct winners     outcome
 *     47-61%         5-11                 3 of 3 recovered the injection
 *     14-24%         15-28                mostly wrong
 *     static scene   19-27                no motion present
 *
 * A fragmented vote and a motionless scene look alike, which is the property
 * this file has otherwise lacked: a single window's argmax is equally confident
 * whether or not anything moved, and that is why prominence turned out to be
 * anti-correlated with correctness. Nineteen distinct winners over 49 windows
 * should read as "no consensus" however prominent the leader is.
 *
 * THAT HOLDS FOR SCENE-DRIVEN NOISE ONLY, AND THE EXCEPTION MATTERS. Agreement
 * detects noise that is INDEPENDENT across windows. An artefact produced by the
 * processing rather than by the scene appears identically in every window, so
 * the windows agree about it unanimously and this statistic reports its highest
 * possible confidence. Measured: the phase estimator returns 0.407 Hz for every
 * injection from 0.2 to 0.7 Hz AND for a motionless scene, with 9 of 9 windows
 * agreeing in all seven cases. No threshold helps, because 100% is the ceiling.
 *
 * The only check that sees a common-mode artefact is a null control -- the same
 * processing over a scene known to be motionless, compared -- because such an
 * artefact is by definition identical whether or not anything moved. That is
 * rs_null_static() and mmotion's --null-static. This function does not replace
 * it and must not be presented as doing so: a run wants both, agreement for
 * scattered noise and a null control for coherent noise.
 *
 * NOT A THRESHOLD, and deliberately no threshold is applied here. The 40%
 * boundary above rests on three correct detections from one seed and one fixture
 * family; it is far too little to hard-code, and a caller that wants to gate on
 * it should do so where the choice is visible. Returning the counts and letting
 * the caller decide is the honest interim, exactly as reporting
 * 'out_n_candidates' is for the function above.
 *
 * 'out_n_contiguous' IS WHERE THE AGREEING WINDOWS ARE, not merely how many.
 * It is the size of the largest 4-connected block of agreeing windows on the
 * window grid, and it separates a mode from a coincidence in a way the count
 * alone cannot: twenty-three windows scattered across the scene and
 * twenty-three forming a patch are the same number and not the same evidence.
 * A vibrating structure occupies contiguous ground; chance crossings do not.
 *
 * There is a non-arbitrary floor for it, and it is the one
 * rs_spectrum_best_window()'s header already derives for the candidate count.
 * Windows are laid at a stride of typically half their width, so they overlap
 * and a target large enough to be resolved at all falls inside a 2x2 block at
 * minimum. A largest block below four therefore cannot describe a spatially
 * resolved mode whatever the agreement percentage says. That bound comes from
 * the window geometry rather than from a tuned constant.
 *
 * Returns RS_ERR_RANGE when no window passes the gates, matching
 * rs_spectrum_best_window(), and RS_ERR_ARG on a NULL or empty spectrum. */
resonarsat_status_t rs_spectrum_consensus(const rs_spectrum_t *spec,
                                          double *out_freq,
                                          size_t *out_n_agree,
                                          size_t *out_n_distinct,
                                          size_t *out_n_voting,
                                          size_t *out_n_contiguous);

/* Return the observation ratio implied by a sub-aperture duration and a measured
 * frequency: t_sap divided by that frequency's period, i.e. how many cycles of
 * the motion each sub-look integrates over.
 *
 * WHAT THIS DOES AND DOES NOT MEAN. This was previously documented as a
 * validity threshold -- above 0.5 the sub-look "smears the target away" and the
 * measurement "should not be trusted". That was wrong, and it was wrong in a
 * way that would have thrown away real results.
 *
 * Integrating a sinusoid over a window of length t_sap does not destroy it. It
 * multiplies its amplitude by |sinc(pi*f*t_sap)|, which falls off but stays
 * finite almost everywhere. The frequency survives; only the amplitude is
 * attenuated, and this project already reports amplitude as qualitative for
 * independent reasons.
 *
 * The correction came from checking the rule against ground truth rather than
 * against intuition. Clemente et al. (EuRAD 2025) recover 87 Hz from an idling
 * van using 0.4167 s sub-apertures, and 36 Hz from a ship using 0.5 s ones,
 * both confirmed against triaxial accelerometers. Those are observation ratios
 * of 36.3 and 18.0 -- seventy times what the old threshold permitted. A rule
 * that rejects two accelerometer-validated measurements is not a rule about the
 * physics. Note that 18.0 is an exact integer, which also rules out the sinc
 * response below as a universal law; see there.
 *
 * USE rs_spectrum_subaperture_response() INSTEAD to judge a measurement. The
 * ratio remains useful as the quantity that response is computed from, and
 * because integer values of it are exactly where the response nulls. */
double rs_spectrum_observation_ratio(double t_sap, double freq);

/* Amplitude the sub-aperture window passes at 'freq', as a factor in [0,1].
 *
 * A sub-look averages the motion over its own duration, which for a sinusoid is
 * a sinc weighting: |sin(pi*f*t_sap) / (pi*f*t_sap)|. Multiply a recovered
 * amplitude by the reciprocal to undo it, if amplitudes are wanted at all.
 *
 * THIS APPLIES TO THE OFFSET-TRACKING OBSERVABLE ONLY, and the boundary is
 * worth stating carefully, because getting it wrong twice is how this function
 * came to be documented at such length.
 *
 * The sinc is the response of a window that AVERAGES DISPLACEMENT. That is what
 * a correlation tracker measures: where the patch sits, averaged over the
 * sub-look. In that picture the response nulls at every integer observation
 * ratio, f = k/t_sap, where a whole number of cycles averages to exactly zero.
 * For 0.4167 s sub-apertures that is a comb at 2.4, 4.8, 7.2 Hz and upward.
 *
 * That picture does not describe phase-based micro-Doppler, and the literature
 * says so plainly. Clemente et al. (EuRAD 2025) recover 36 Hz from a ship using
 * 0.5 s sub-apertures. The observation ratio is exactly 18.0 -- an integer,
 * where the sinc above predicts zero response -- and the measurement is
 * confirmed against triaxial accelerometers. A vibration modulates the PHASE,
 * which puts its energy into micro-Doppler sidebands in the azimuth spectrum;
 * far above 1/t_sap that sideband structure, not the averaged displacement, is
 * what carries the frequency.
 *
 * So the honest reading is a regime boundary rather than a limit:
 *
 *   f well below 1/t_sap   the sub-look tracks quasi-static displacement, the
 *                          sinc describes the amplitude, and this pipeline's
 *                          correlation tracker is in its element.
 *
 *   f near k/t_sap         the AVERAGING observable is blind there. A peak from
 *                          RS_MICROM_EST_CORRELATION beside such a null is
 *                          suspect; a phase measurement need not be.
 *
 *   f well above 1/t_sap   sidebands dominate and the averaging model stops
 *                          applying. Reaching this regime needs the phase
 *                          observable, which this interface does not yet
 *                          provide -- see the note at the top of this file.
 *
 * Returns 1.0 at zero frequency and for a non-positive t_sap, since neither
 * attenuates anything. */
double rs_spectrum_subaperture_response(double t_sap, double freq);

/* Return the number of sub-looks needed for a target of the given vibration
 * frequency and line-of-sight amplitude to sit inside the phase-ambiguity
 * interval, for a collect of the given dwell and wavelength.
 *
 * The condition comes from requiring the peak azimuth shift to fall within three
 * quarters of a sub-look resolution cell:
 *
 *     peak shift  = 2*pi*f*A * R/V
 *     resolution  = lambda*R*M / (2*V*T_dwell),  M = 1 + (N-1)(1-overlap)
 *     require       peak shift < 0.75 * resolution
 *     giving        M > (4*pi/0.75) * f * A * T_dwell / lambda
 *
 * and this returns the N implied for the given overlap. Note the direction,
 * which is the opposite of intuition: SHORTER sub-looks are needed, because a
 * coarser resolution cell makes the ambiguity interval WIDER in metres while the
 * target's shift is fixed by the geometry. Overlap works against the condition,
 * so zero overlap needs the fewest looks.
 *
 * This is a lower bound from one constraint, not an optimum. Pushing the look
 * count far past it eventually starves each sub-look of bandwidth until the
 * resolution cell exceeds the correlation window, and measured performance falls
 * away again. Use rs_vibration_fmax() and rs_observation_ratio() to check the
 * other two constraints, and sweep around the result.
 *
 * Returns 0 for degenerate inputs. */
size_t rs_microm_recommend_looks(double vib_freq, double amp_los,
                                 double t_dwell, double lambda, double overlap);

/* Return the largest line-of-sight velocity in m/s whose azimuth shift a
 * correlation window of 'win_px' pixels can measure without wrapping, given the
 * pixel spacing, slant range and platform speed.
 *
 * The azimuth shift of a radially moving target is dx = R*v_r/V, and a
 * correlation surface over an N-pixel window is unambiguous only over +/-N/2
 * pixels. Beyond that the measurement folds and reports a second harmonic. Which
 * limit applies depends on the reference mode: against look 0 the full excursion
 * must fit, whereas against the adjacent look only the change between samples
 * must, which is smaller by 2*pi*f*dt.
 *
 * Returns 0.0 for degenerate geometry. */
double rs_microm_max_velocity(size_t win_px, double spacing_m,
                              double slant_range, double v_platform);

#endif /* RESONARSAT_MICROM_H */
