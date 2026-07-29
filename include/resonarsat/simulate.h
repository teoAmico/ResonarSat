/* A static scene with a real collect's geometry, for null testing.
 *
 * WHY THIS EXISTS, AND WHY THE EXISTING NULL WAS NOT ENOUGH.
 *
 * The shuffled null permutes the time order of ALREADY FORMED sub-looks. That
 * controls for a great deal -- scene, brightness, coherence and look count are
 * preserved exactly -- and it is the right test for a decomposition whose
 * sub-looks are independent.
 *
 * It is the wrong test for one whose sub-looks overlap. The decomposition of
 * Biondi & Malanga (2022) section 3.1 holds out a fraction of the Doppler band
 * and steps the remainder across the spectrum in rigid shifts, so adjacent
 * sub-looks share most of their bandwidth and therefore most of their speckle
 * realisation. Their measured shifts are correlated BEFORE any target moves.
 * That correlation makes the shift series smooth, a smooth series has its energy
 * at low frequency, and shuffling destroys the smoothness. So the unshuffled
 * series outscores its shuffles for a reason that has nothing to do with the
 * ground, and the test reports significance on a scene where nothing moves.
 *
 * Measured on the Panama collect: the paper decomposition returned p = 0.08 at
 * a frequency whose period is three times shorter than the sub-aperture that
 * measured it. A sub-look integrating three full cycles averages that motion
 * away, so the peak cannot be the motion it appears to be. The overlap is the
 * candidate explanation, and the shuffled null cannot rule it in or out.
 *
 * WHAT THIS DOES INSTEAD. It synthesises phase history for a scene of purely
 * STATIC scatterers, using the reference collect's own pulse times, platform
 * positions, carrier and range sampling. Running the identical chain over it --
 * same sub-aperture route, same tracker, same spectral estimator -- produces the
 * prominence a perfectly motionless world yields through this exact processing.
 * Overlap, tracker bias and estimator behaviour are all inherited, because the
 * same code computes them. Only the motion is absent.
 *
 * A measurement that does not beat that distribution has not demonstrated
 * motion, whatever it does against a shuffle. */

#ifndef RESONARSAT_SIMULATE_H
#define RESONARSAT_SIMULATE_H

#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/readers.h"

/* Synthesise a static-scene collect with 'ref's geometry.
 *
 * Pulse times, platform positions, per-pulse reference ranges, carrier and
 * range bin spacing are copied from 'ref', so the simulated aperture, dwell and
 * Doppler history are the real ones rather than a plausible substitute. The
 * scene is 'n_target' point scatterers placed pseudo-randomly within
 * 'extent_m' metres of the grid origin, with Rayleigh-distributed reflectivity
 * so the focused image has realistic speckle statistics. None of them moves.
 *
 * 'n_rbin' caps the simulated range extent, which the caller should set just
 * large enough to cover its processing grid: copying the reference's full swath
 * would allocate gigabytes to synthesise range bins no window ever reads. Pass
 * 0 to use the reference's own count. The near range is recentred so the grid
 * origin sits mid-swath whatever the count.
 *
 * 'seed' selects the realisation. Distinct seeds give independent speckle over
 * the same geometry, which is what makes a distribution of prominences rather
 * than a single number.
 *
 * On success '*out' owns its arrays and must be released with rs_cphd_free().
 * Returns RS_ERR_ARG on a NULL argument or a reference carrying no pulses, and
 * RS_ERR_ALLOC if the phase history cannot be sized. */
resonarsat_status_t rs_simulate_static_like(const rs_cphd_t *ref, unsigned seed,
                                            size_t n_target, const double centre[2],
                                            double extent_m,
                                            size_t n_rbin, rs_cphd_t *out);

#endif /* RESONARSAT_SIMULATE_H */
