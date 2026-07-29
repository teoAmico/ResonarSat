/* Split-band Phase Linking: shift estimation across a stack of sub-looks.
 *
 * This implements the estimator of De Zan, "Coherent Shift Estimation for Stacks
 * of SAR Images" (IEEE GRSL 8, 1095, 2011), which is specified for exactly the
 * configuration this project has -- N partially coherent images of one scene,
 * shifts wanted between them -- and which comes within 0.5 dB (a loss factor of
 * 8/9) of the Cramer-Rao bound.
 *
 * WHY THIS EXISTS ALONGSIDE THE CORRELATION TRACKER. The same paper
 * characterises cross-correlation peak finding, which src/core/coreg.c does, as
 * achieving "performances in the order of the resolution element for a few
 * independent samples". That matches this project's measurements: the tracking
 * noise is about 1.2 sub-look resolution cells. The noise is therefore not a
 * defect to be hunted -- it is what the method delivers -- and reaching further
 * requires a different estimator rather than a better-tuned one.
 *
 * THE METHOD, in the paper's terms:
 *   1. Filter the lower and upper third of each signal's bandwidth, giving two
 *      sets of N vectors that each retain a third of the independent samples.
 *   2. From each set estimate N phases (N-1 really, one being the reference) by
 *      Phase Linking: the maximum-likelihood phase estimator over all N^2
 *      interferograms formable from the stack, given the coherence matrix.
 *   3. Difference the two phase vectors and scale by 3/(4*pi), which inverts the
 *      phase-delay relation phi = 2*pi*f*d for a sub-band separation of 2/3 of
 *      the bandwidth.
 *
 * WHAT THIS CANNOT DO, and it is a property of the problem rather than of the
 * implementation: only RELATIVE shifts are estimable. The Fisher information
 * matrix is rank deficient, with the all-ones vector always an eigenvector of
 * eigenvalue zero, so no information exists about a delay affecting every image
 * equally. Every shift returned here is relative to the reference look. */

#ifndef RESONARSAT_PHASELINK_H
#define RESONARSAT_PHASELINK_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"

/* Estimate the phase of each of 'n_sig' signals relative to signal 0, by Phase
 * Linking over the stack's coherence matrix.
 *
 * 'sig' holds 'n_sig' vectors of 'n_samp' complex samples each, laid out
 * signal-major: signal i occupies sig[i*n_samp .. i*n_samp + n_samp).
 *
 * The sample coherence matrix is formed from all pairs, then the phase vector is
 * found by the fixed-point iteration x_n <- exp(j*arg(sum_m Gamma_nm x_m)),
 * started from the all-ones vector. That iteration is the standard practical
 * route to the maximum-likelihood phase estimate over a stack; it converges in a
 * handful of steps for the coherence levels seen here, and the iteration count
 * is capped so a pathological input cannot spin.
 *
 * The estimate uses every pair, not just pairs against the reference. That is
 * the whole point: with N looks there are N^2 interferograms carrying phase
 * information, and using only the N-1 formed against one reference discards most
 * of it -- which is precisely what the correlation tracker does.
 *
 * Phases are written to 'phase_out' in radians, with phase_out[0] == 0 by
 * construction. Returns RS_ERR_ARG for degenerate sizes and RS_ERR_ALLOC on
 * memory failure. A signal with no energy contributes nothing and receives a
 * phase of zero rather than a NaN. */
resonarsat_status_t rs_phase_link(const float complex *sig,
                                  size_t n_sig, size_t n_samp,
                                  double *phase_out);

/* Estimate the shift of each sub-look relative to sub-look 0, in samples, by the
 * split-band Phase Linking method.
 *
 * 'patch' holds 'n_look' patches of 'n_az' by 'n_rg' complex samples, patch-major
 * and row-major within a patch. Shifts are measured along the azimuth (first)
 * axis, which is the along-track direction and the one carrying micro-motion.
 *
 * A subtlety the source does not have to address: it assumes signals scaled to
 * unit bandwidth, whereas these sub-looks are heavily OVERSAMPLED -- a 2.5 m
 * resolution sampled at 0.5 m occupies only about a fifth of the sampled band.
 * Splitting the sampled band into thirds would put the outer thirds in empty
 * spectrum and return noise. The occupied fraction is therefore measured from
 * the data (the band holding 'RS_PL_ENERGY_FRAC' of the azimuth power, averaged
 * over looks and range) and the thirds are taken within that, with the delay
 * scaling adjusted by the same fraction.
 *
 * Returns the shifts in 'shift_out' (n_look doubles, shift_out[0] == 0) and, when
 * non-NULL, the mean off-diagonal coherence in 'coherence_out' as a quality
 * measure comparable to the correlation tracker's peak value.
 *
 * Phase ambiguity bites when the shift exceeds three quarters of a resolution
 * cell of the sub-band separation; for the geometries here that is metres, well
 * above the motion of interest, but it is not unlimited. */
resonarsat_status_t rs_splitband_shift(const float complex *patch,
                                       size_t n_look, size_t n_az, size_t n_rg,
                                       double *shift_out,
                                       double *coherence_out);

/* Fraction of azimuth spectral energy used to delimit the occupied band. */
#define RS_PL_ENERGY_FRAC 0.90

#endif /* RESONARSAT_PHASELINK_H */
