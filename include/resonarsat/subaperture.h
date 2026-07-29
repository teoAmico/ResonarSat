/* Doppler sub-aperture decomposition.
 *
 * Each azimuth Doppler sub-band is an image of the scene seen at a slightly
 * different time within the synthetic aperture, so a stack of them is a movie:
 * a vibrating structure modulates its apparent position from look to look. This
 * is the stage that turns one SAR image into a time series, and everything
 * downstream measures that series. */

#ifndef RESONARSAT_SUBAPERTURE_H
#define RESONARSAT_SUBAPERTURE_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/slc.h"

/* Forward declarations: the phase-history path needs these types but this
 * header must not pull in focus.h, which includes this one. */
struct rs_cphd;
struct rs_grid;

/* How the aperture is divided.
 *
 * MASTER-SLAVE PAIRING -- THIS DEPARTS FROM THE SOURCES. The paper's Eqs. 4 and
 * 5 define two matrices, master and slave, each of N_D elements; Figure 7 draws
 * them as two rectangles offset by B_shift stepping across the azimuth band
 * together; and the patent's Figure 0.5 makes blocks 3-6 two band-pass filters
 * and two IDFT2s feeding one pixel-tracking block. All say the displacement at
 * step m is master(m) against slave(m).
 *
 * This instead forms one stepped sequence and compares every look to look 0.
 * That reading rested on one phrase in reference [52] -- "the first SAR image
 * (master)" -- which is most plausibly the first OF THE PAIR: the same paper's
 * method sentence one paragraph earlier correlates "adjacent Doppler
 * sub-apertures", and ITS source names the strategy outright as "the
 * 'small-temporal baseline' with the sliding master" (Biondi et al. 2019, RS
 * 11(14):1637, sec. 2.3). The patent read directly (WO2024008365A1 [0004]) is
 * unambiguous: "N_D rigid shifts of the master-slave system", the pair "rigidly
 * held at a distance B_shift", tracking the slave's position against its own
 * master. No source describes a fixed common reference. The sources are against
 * what this does; it is kept because changing it changes what the displacement
 * series IS, and that is a decision to take deliberately rather than silently.
 *
 * The two choices are coupled, and the coupling is total rather than
 * approximate. With B_shift pinned to the step (below), slave(m) is centred
 * exactly where master(m+1) is, and both are cut from the same spectrum with
 * the same width and taper, so the two are the SAME IMAGE sample for sample --
 * measured at zero difference in test_subaperture.c, not argued. The sources'
 * scheme therefore does not merely resemble tracking consecutive looks under
 * this default; it is that, minus the accumulation.
 *
 * AND THE LAYOUT CANNOT ESCAPE UPWARD. N_D masters of width B_CD-B_DL stepping
 * across the band leave exactly one step of headroom at the top, so a B_shift
 * larger than the step puts the last slave outside the Doppler support and is
 * refused. Every separation that makes the slave an independent band is
 * therefore SMALLER than the step -- a shorter lag, and a differential
 * attenuated further by |2 sin(pi f dt)|. That runs against the one thing
 * WO2024008365A1 [0004] says about choosing the parameter: "the higher this
 * parameter, the lower the mechanical frequency observed". The patent's sweep
 * geometry and its stated tuning knob are in conflict as soon as B_shift
 * exceeds the step, and it never reconciles them.
 *
 * B_SHIFT IS FREE IN THE SOURCES AND DERIVED HERE. WO2024008365A1 [0004] keeps
 * the two bands "rigidly held at a distance B_shift", chosen to select "the
 * precise vibrational frequency one wishes to observe", adding only that the
 * higher it is, the lower the frequency observed. The patent treats it as
 * independent of the sweep, which advances the rigid pair separately in N_D
 * steps. But it gives no value and no selection rule for B_shift anywhere,
 * including in its claims -- it states a best value only for the held-out band,
 * B_CL = B_CD/2 -- so there is no defensible value to set it to.
 *
 * This therefore pins B_shift to the printed Doppler step (B_CD - B_DL)/N_D,
 * the one spacing the sources do define, and records it in rs_subap_stack_t.
 * The cost is that the pair cannot be tuned to a chosen frequency, which may be
 * the mechanism the method depends on.
 *
 *   RS_SUBAP_PAPER    Biondi & Malanga (2022) section 3.1 and figure 7. A
 *                     fraction B_DL of the total Doppler band is deliberately
 *                     left out of the matched filter "to obtain a sufficient
 *                     sensitivity to estimate target motions", and the
 *                     remaining band is stepped across the spectrum in N_D
 *                     rigid shifts. Sub-looks overlap heavily by construction.
 *
 *   RS_SUBAP_UNIFORM  A plain filter bank: N equal sub-bands with a specified
 *                     fractional overlap. Simpler to reason about, matches the
 *                     independently validated processing literature, and is the
 *                     right default for development and cross-checking.
 */
typedef enum {
    RS_SUBAP_PAPER = 0,
    RS_SUBAP_UNIFORM = 1
} rs_subap_mode_t;

/* Sub-aperture decomposition parameters.
 *
 * 'left_out_frac' applies to RS_SUBAP_PAPER only and defaults to 0.5, the
 * value the paper uses. See the warning in rs_subaperture_split() about why
 * changing it is not as innocuous as it looks. */
typedef struct {
    rs_subap_mode_t mode;
    size_t n_looks;         /* N_D: number of sub-looks to form */
    double overlap;         /* uniform mode only: fractional band overlap [0,1) */
    double left_out_frac;   /* paper mode only: B_DL / B_CD, default 0.5 */
    int    window;          /* nonzero: raised-cosine band edges (strongly advised) */

    /* B_shift, the master-slave separation in Hz. Paper mode only.
     *
     * WO2024008365A1 [0004] makes this free and independent of the sweep: the
     * slave is the master "slightly shifted along the azimuth frequency by
     * B_shift, which represents the precise vibrational frequency one wishes to
     * observe. The higher this parameter, the lower the mechanical frequency
     * observed." The sweep is a separate motion, advancing the rigid pair in
     * N_D steps of (B_CD-B_DL)/N_D.
     *
     * Zero means derive it from the step, which is what this code did
     * exclusively before the patent was read directly. That remains the default
     * because no source gives a value or a selection rule for B_shift -- the
     * patent states a best value only for the held-out band, B_CL = B_CD/2.
     *
     * The frequency this selects follows from the spectral gap being a time lag,
     * dt = B_shift * t_dwell / B_CD: a pair separated by dt is blind at DC and
     * most sensitive near 1/(2*dt), so the observed frequency goes as 1/B_shift,
     * which is the patent's sentence. */
    double b_shift_hz;

    /* Nonzero: also build the slave band of each pair, offset by B_shift, so a
     * tracker can measure slave(i) against master(i) as the sources describe.
     * Paper mode only; doubles the stack's memory. */
    int    pair;
} rs_subap_params_t;

/* A stack of sub-look images sharing one geometry.
 *
 * 'look' holds 'n_looks' images of identical dimensions. 'centre_time' gives
 * each look's centre time within the aperture in seconds, which is what makes
 * the stack a time series rather than an unordered set. 'dt' is the sampling
 * interval between consecutive looks, 'f_max' the Nyquist limit it implies, and
 * 'az_resolution' the azimuth resolution each individual look achieves.
 *
 * Those last two travel together deliberately. Reporting an observable
 * vibration band without the resolution it cost is how a 30-metre result gets
 * mistaken for a 1-metre one. */
typedef struct {
    rs_slc_t *look;
    /* The slave of each pair, offset from look[i] by b_shift_hz, or NULL when
     * params.pair was not requested. Same length as 'look'. */
    rs_slc_t *slave;
    double   *centre_time;
    size_t    n_looks;

    double dt;             /* s between consecutive look centres */
    double f_max;          /* Hz, vibration Nyquist limit = 1/(2*dt) */
    double t_sap;          /* s, duration of one sub-aperture */
    double az_resolution;  /* m, azimuth resolution of one sub-look */
    double doppler_centroid;  /* Hz, centroid removed before splitting */
    double doppler_bandwidth; /* Hz, total band the split covered */
    double b_shift_hz;        /* Hz, paper-mode B_shift actually used, or 0 */
    double pair_lag_s;        /* s, time lag the B_shift gap corresponds to, or 0 */

    rs_subap_params_t params;
} rs_subap_stack_t;

/* Fill 'params' with the defaults this project recommends: uniform mode, 16
 * looks, 40 percent overlap, raised-cosine windowing on. The overlap figure
 * comes from the published working range of 38 to 43 percent. */
void rs_subap_params_default(rs_subap_params_t *params);

/* Decompose a focused image into a sub-look stack.
 *
 * The image is transformed along azimuth column by column, the Doppler centroid
 * is estimated and removed, band filters are applied per the selected mode, and
 * each band is transformed back. On success '*stack' owns every sub-look image
 * and must be released with rs_subap_stack_free().
 *
 * Band edges are raised-cosine when params->window is set, and callers should
 * leave it set. A rectangular band produces sinc sidelobes in the sub-look
 * time series, which alias directly into the vibration band the next stage
 * measures and manufacture spectral peaks that look exactly like structural
 * modes. This is one of the few places in the pipeline where the cheap choice
 * is actively misleading rather than merely imprecise.
 *
 * On the paper mode's step size: the source states both that B_DL "is divided
 * into N_D equally distributed bandwidth steps" and that the shift per step is
 * (B_CD - B_DL)/N_D. These agree only when B_DL = B_CD/2, which is the value it
 * uses. This implementation takes the printed formula, (B_CD-B_DL)/N_D, and
 * warns via rs_set_error() if a caller supplies a left_out_frac other than 0.5
 * where the two readings would diverge -- silently picking one would make the
 * decomposition depend on an ambiguity the caller cannot see.
 *
 * Returns RS_ERR_ARG for a stack of zero looks or an image too short in azimuth
 * to divide, and RS_ERR_MISSING_META if the image lacks azimuth timing. */
resonarsat_status_t rs_subaperture_split(const rs_slc_t *img,
                                         const rs_subap_params_t *params,
                                         rs_subap_stack_t *stack);

/* Release every sub-look image and the stack's own arrays. Safe on a
 * zero-initialised or partially constructed stack. */
void rs_subap_stack_free(rs_subap_stack_t *stack);

/* Form a sub-look stack directly from phase history, by focusing shifted pulse
 * windows.
 *
 * This is the faithful implementation of the method and the one to prefer
 * whenever phase history is available. Biondi & Malanga (2022) section 3.1
 * describes sub-bands as "generated by focusing the SAR image, where the
 * matched filter is set to exploit a range-azimuth bandwidth equal to
 * (B_cr, B_cD - B_DL)": the looks come from re-focusing over restricted
 * apertures, and that is exactly what this does -- N_D calls to
 * rs_focus_backproject() over windows of pulses stepped along the collect.
 *
 * Splitting an already-focused image spectrally, as rs_subaperture_split()
 * does, only approximates this. It is the right function when the input is a
 * focused product and there is no phase history to go back to, but it inherits
 * whatever the original focusing did, and the mapping from spectral band to
 * time is indirect. Here the timing is exact by construction: each look's
 * duration and centre time are read straight off the pulse times.
 *
 * 'params' supplies n_looks and, in uniform mode, the fractional overlap
 * between consecutive pulse windows; the mode's band-layout fields are not used
 * because the aperture is divided in time rather than in frequency.
 *
 * Cost is n_looks backprojections, so it is the most expensive operation in the
 * pipeline by a wide margin. Work on cropped patches.
 *
 * Returns RS_ERR_ARG if the requested overlap or look count cannot divide the
 * available pulses into windows of at least one pulse. */
resonarsat_status_t rs_subaperture_from_cphd(const struct rs_cphd *cphd,
                                             const struct rs_grid *grid,
                                             const rs_subap_params_t *params,
                                             rs_subap_stack_t *stack);

/* Estimate the Doppler centroid of an image in Hz, by locating the centre of
 * mass of the azimuth power spectrum averaged over range.
 *
 * This is the fallback for products whose annotated Doppler polynomial is
 * absent or untrustworthy, and the cross-check for those where it is present.
 * A large disagreement between the two is worth surfacing rather than
 * averaging away: it usually means the region of interest straddles a burst
 * boundary, which invalidates the sub-aperture assumption entirely.
 *
 * Returns the estimate through 'out'. */
resonarsat_status_t rs_estimate_doppler_centroid(const rs_slc_t *img, double *out);

/* Measure the Doppler bandwidth of an image in Hz at the given power fraction
 * below peak ('frac' of 0.5 gives the conventional -3 dB width).
 *
 * Averages the azimuth power spectrum over range, then walks outward from the
 * peak to the first crossing on each side. Reported to let the caller check a
 * product against its expected bandwidth before trusting any decomposition of
 * it. */
resonarsat_status_t rs_estimate_doppler_bandwidth(const rs_slc_t *img, double frac, double *out);

#endif /* RESONARSAT_SUBAPERTURE_H */
