/* Cross-correlation primitives used by the micro-motion stage. */

#ifndef RESONARSAT_COREG_H
#define RESONARSAT_COREG_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/slc.h"

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
 * division by zero, so a blank region of a scene masks itself out. */
resonarsat_status_t rs_coreg_shift(const float complex *ref, const float complex *img,
                                   size_t n_az, size_t n_rg,
                                   size_t upsample_az, size_t upsample_rg,
                                   double *shift_az, double *shift_rg, double *peak);

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
