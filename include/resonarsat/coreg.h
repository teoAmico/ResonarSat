/* Cross-correlation primitives used by the micro-motion stage. */

#ifndef RESONARSAT_COREG_H
#define RESONARSAT_COREG_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/slc.h"

/* How the sub-pixel peak is located on the correlation surface.
 *
 * The two modes search the SAME lattice -- 1/upsample_az by 1/upsample_rg
 * pixels -- and differ only in extent. That is what makes them comparable: a
 * disagreement is a disagreement about where the peak is, not about how finely
 * the surface was sampled.
 *
 * RS_COREG_REFINE_LOCAL is the production path. It takes the integer peak of the
 * full surface and then evaluates the fine lattice within one pixel of it.
 *
 * RS_COREG_REFINE_EXHAUSTIVE is the audit baseline reached by --no-optimize. It
 * evaluates every point of the fine lattice over the WHOLE surface, by
 * zero-padding the cross-spectrum to (n_az*upsample_az) by (n_rg*upsample_rg)
 * and inverse transforming, then taking the global maximum. This is the
 * textbook, uninteresting way to do it and it costs the memory to match:
 * upsample factors of 10 by 20 on a 24 by 24 patch build a 240 by 480 surface
 * per call per thread.
 *
 * WHAT THE COMPARISON CAN AND CANNOT SHOW. The local path's integer peak is
 * already the GLOBAL maximum of the sampled surface -- it is found by a full
 * inverse transform and a scan over every bin (see rs_coreg_shift()), not by a
 * local search. So the exhaustive mode cannot detect a missed distant lobe in
 * the sampled surface; it can only detect the case where the true continuous
 * peak lies more than one pixel from the strongest SAMPLE. That happens when the
 * surface is not adequately sampled: two comparable lobes a few pixels apart,
 * where interpolation between the samples crests higher near the weaker one.
 * On decorrelated sub-looks with several competing scatterers this is a real
 * possibility, and before this mode existed there was no way to measure how
 * often it occurred. Expect the two modes to agree on almost every window; the
 * windows where they do not are the interesting output. */
typedef enum {
    RS_COREG_REFINE_LOCAL = 0,  /* one pixel about the integer peak (default) */
    RS_COREG_REFINE_EXHAUSTIVE  /* whole zero-padded surface, global maximum */
} rs_coreg_refine_t;

/* Estimate the shift between two equally sized complex patches.
 *
 * 'ref' and 'img' are 'n_az' by 'n_rg' row-major patches with their means
 * already removed (rs_coreg_extract() does that). The shift that best aligns
 * 'img' onto 'ref' is written to '*shift_az' and '*shift_rg' in pixels, and the
 * normalised correlation coefficient at that shift to '*peak', in [0, 1].
 *
 * The estimate is refined below one pixel by evaluating the correlation surface
 * on a grid of 1/upsample_az by 1/upsample_rg spacing spanning one pixel either
 * side of the integer peak, following Guizar-Sicairos et al. (2008): only the
 * neighbourhood of the peak is upsampled, never the whole surface. Published
 * working values for this pipeline's data are 10 in azimuth and 20 in range.
 *
 * Patches with no variance yield a zero shift and a zero peak rather than a
 * division by zero, so a blank region of a scene masks itself out.
 *
 * Equivalent to rs_coreg_shift_ex() with RS_COREG_REFINE_LOCAL. */
resonarsat_status_t rs_coreg_shift(const float complex *ref, const float complex *img,
                                   size_t n_az, size_t n_rg,
                                   size_t upsample_az, size_t upsample_rg,
                                   double *shift_az, double *shift_rg, double *peak);

/* As rs_coreg_shift(), with the peak-search extent selectable.
 *
 * 'refine' picks the search strategy; see rs_coreg_refine_t for what the choice
 * does and does not buy. Both modes report the shift on the same 1/upsample
 * lattice and normalise '*peak' identically, so a caller may switch between them
 * and compare the numbers directly.
 *
 * RS_COREG_REFINE_EXHAUSTIVE returns RS_ERR_RANGE if the padded surface would
 * exceed RS_COREG_MAX_SURFACE elements, rather than attempting the allocation --
 * the product of four caller-supplied sizes is easy to make enormous by
 * accident. */
resonarsat_status_t rs_coreg_shift_ex(const float complex *ref, const float complex *img,
                                      size_t n_az, size_t n_rg,
                                      size_t upsample_az, size_t upsample_rg,
                                      rs_coreg_refine_t refine,
                                      double *shift_az, double *shift_rg, double *peak);

/* Report whether an exhaustive search of this size is permitted.
 *
 * Returns RS_OK if the zero-padded surface for 'n_az' by 'n_rg' patches at
 * 'upsample_az' by 'upsample_rg' fits within RS_COREG_MAX_SURFACE, and
 * RS_ERR_RANGE with a description naming the offending sizes otherwise.
 *
 * Exposed so a caller can refuse an impossible configuration up front instead of
 * discovering it once per window inside a tracking loop. That distinction
 * matters: a tracker treats a failed correlation as "this window did not track",
 * so a size error raised per window returns a complete result in which
 * everything is zero and nothing says why. rs_microm_track() calls this before
 * it allocates. */
resonarsat_status_t rs_coreg_surface_check(size_t n_az, size_t n_rg,
                                          size_t upsample_az, size_t upsample_rg);

/* Ceiling on the zero-padded correlation surface, in complex samples.
 *
 * 2^24 elements is 128 MB at single precision, per call and per thread. Well
 * above anything the pipeline asks for (10 by 20 upsampling on a 24 by 24 patch
 * needs 115200) and low enough that a mistyped upsample factor fails with a
 * message instead of exhausting memory. */
#define RS_COREG_MAX_SURFACE ((size_t)1 << 24)

/* Copy the 'n_az' by 'n_rg' patch whose top-left corner is at (az0, rg0) out of
 * an image into 'patch', subtracting the patch mean as it goes.
 *
 * Mean removal turns the subsequent correlation into a covariance, so that a
 * brightness offset between two sub-looks does not bias the peak. Returns
 * RS_ERR_ARG, with a description naming the offending extent, if any part of
 * the patch would fall outside the image. */
resonarsat_status_t rs_coreg_extract(const rs_slc_t *img, size_t az0, size_t rg0,
                                     size_t n_az, size_t n_rg, float complex *patch);

#endif /* RESONARSAT_COREG_H */
