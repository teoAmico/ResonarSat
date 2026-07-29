/* FFT wrapper.
 *
 * This is the only interface through which the rest of ResonarSat performs a
 * Fourier transform. No file outside src/core/fft.c may include, link against
 * or name the underlying library. That rule exists for two reasons: it keeps a
 * backend swap (to a GPU library, or back to FFTW for benchmarking) a one-file
 * change, and it keeps the project's licence honest -- ResonarSat is MIT, and
 * linking a GPL FFT into it would encumber the distributed binary while the
 * LICENSE file promised otherwise. The current backend is PocketFFT
 * (BSD-3-Clause), vendored under third_party/. */

#ifndef RESONARSAT_FFT_H
#define RESONARSAT_FFT_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"

/* Opaque handle for a transform of one fixed length. Creating a plan is the
 * expensive part; a plan is reusable across any number of transforms of that
 * length and callers are expected to hoist creation out of their loops. */
typedef struct rs_fft_plan rs_fft_plan;

/* Create a plan for complex transforms of length 'n' and store it in '*out'.
 *
 * 'n' may be any positive length, not only a power of two; the backend handles
 * arbitrary factorisations, though highly composite lengths are fastest. On
 * success '*out' owns memory that rs_fft_plan_destroy() must release. On
 * failure '*out' is set to NULL.
 *
 * A plan is NOT thread-safe for concurrent transforms. Under OpenMP, give each
 * thread its own plan rather than sharing one; plan creation is cheap relative
 * to the per-pixel work in this pipeline. */
resonarsat_status_t rs_fft_plan_create(size_t n, rs_fft_plan **out);

/* Release a plan. Accepts NULL, so it is safe on an error path where creation
 * may not have happened. */
void rs_fft_plan_destroy(rs_fft_plan *plan);

/* Return the transform length a plan was created for. */
size_t rs_fft_plan_length(const rs_fft_plan *plan);

/* Transform 'data' in place, forward (time to frequency), unnormalised.
 *
 * 'data' must hold exactly rs_fft_plan_length(plan) elements. The transform is
 * unnormalised in the forward direction and normalised by 1/n in the inverse,
 * so that a forward followed by an inverse returns the original samples; this
 * is the convention the rest of the pipeline assumes.
 *
 * Note on precision: the pipeline carries single-precision complex data, while
 * the current backend transforms in double. This function converts at its
 * boundary. The conversion costs memory traffic but keeps arithmetic accuracy
 * comfortably above what the data warrants, and it is invisible to callers --
 * which is the point of the wrapper. */
resonarsat_status_t rs_fft_forward(rs_fft_plan *plan, float complex *data);

/* Transform 'data' in place, inverse (frequency to time), normalised by 1/n.
 * See rs_fft_forward() for the length requirement and the normalisation
 * convention. */
resonarsat_status_t rs_fft_inverse(rs_fft_plan *plan, float complex *data);

/* Rotate a spectrum in place so that the zero-frequency component moves from
 * index 0 to the centre, matching the layout most people expect when they plot
 * a spectrum and the layout this project's band filters assume.
 *
 * For even 'n' this is a simple halves swap and is its own inverse. For odd 'n'
 * it is not: use rs_ifft_shift() to undo it. Getting this wrong on one axis of
 * a two-dimensional transform, and not the other, is a classic source of
 * results that look almost right. */
void rs_fft_shift(float complex *data, size_t n);

/* Exact inverse of rs_fft_shift(), moving the centre component back to index 0.
 * Identical to rs_fft_shift() for even 'n' and different for odd 'n'. */
void rs_ifft_shift(float complex *data, size_t n);

/* Transform a 2-D array of 'n_row' by 'n_col' complex samples in place, row
 * major, forward if 'inverse' is zero and inverse otherwise.
 *
 * Rows are transformed first, then columns. The inverse direction is normalised
 * by 1/(n_row*n_col) in total, so a forward-then-inverse round trip is the
 * identity. Used by the cross-correlation path; the azimuth-only transforms of
 * the sub-aperture stage use rs_fft_forward() on strided columns instead.
 *
 * Allocates a scratch column buffer internally and returns RS_ERR_ALLOC if that
 * fails, leaving 'data' partially transformed -- callers on an error path
 * should discard the array rather than trying to recover it. */
resonarsat_status_t rs_fft2(float complex *data, size_t n_row, size_t n_col, int inverse);

#endif /* RESONARSAT_FFT_H */
