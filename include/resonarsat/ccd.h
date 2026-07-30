/* Scale-invariant coherent change detection over a sub-aperture stack.
 *
 * This stage LOCATES micro-motion. It does not measure it. Everything else in
 * this project produces a displacement or velocity series and a spectrum from
 * it; this produces one map answering "where in this scene is something
 * moving", with no statement about frequency or amplitude.
 *
 * That is worth having because every measurement made here so far has swept
 * hundreds of windows with no prior on where a peak should be, and the number of
 * tries is itself a source of false peaks -- see the warning beside the
 * eligible-window count in mmotion, where the identical chain refused honestly
 * at 225 windows and reported a confident frequency at 961 off two threshold
 * crossings. A locator supplies the prior, so the expensive estimator can be
 * pointed rather than swept.
 *
 * THE METHOD. Rollo, Vattulainen, Ilioudis, Milillo and Clemente, "Scale
 * Invariant Coherent Change Detection to Locate Micro-Motion in Single Pass SAR
 * Images", apply the maximal-invariant change detector of Carotenuto et al.
 * (IEEE TGRS 2015, 2016) to consecutive sub-apertures in place of polarimetric
 * channels. For a W x W window around each pixel and each sub-aperture index m:
 *
 *     X = channels {look[m-1], look[m]}     Y = channels {look[m], look[m+1]}
 *     R_X = 2 x W^2, columns are the per-pixel channel vectors
 *     S_X = R_X R_X^H  (2x2 Hermitian)      S_Y likewise
 *     lambda_1 >= lambda_2 = eigenvalues of S_X S_Y^-1
 *     G = lambda_1 / lambda_2
 *
 * accumulated over m.
 *
 * WHY SCALE INVARIANCE IS THE POINT. The null hypothesis is Sigma_X = gamma *
 * Sigma_Y for ANY gamma > 0, so a target whose brightness merely changes between
 * sub-apertures does not trigger. Every bright stationary target does change
 * that way, because each sub-look sees it from a different part of the aperture,
 * and that is the false alarm an amplitude-based detector cannot suppress. The
 * algebra is immediate: if S_X = gamma S_Y then S_X S_Y^-1 = gamma I, both
 * eigenvalues are gamma, and G is exactly 1 whatever gamma is.
 *
 * BE PRECISE ABOUT WHAT THAT DOES NOT SAY, because the obvious test of it is
 * wrong and was written here before being caught. The invariance is between the
 * two COVARIANCE MATRICES, not between individual sub-looks. Scaling one look
 * multiplies a single CHANNEL of each matrix it appears in, which is not a
 * scalar multiple of the matrix and does not leave G unchanged -- and since the
 * X and Y pairs share look[m] by construction, no single-look scaling can probe
 * the identity at all. What is exactly testable at map level is a scaling of the
 * WHOLE stack: every covariance scales together, S_X S_Y^-1 is untouched, and
 * the map is unchanged bit for bit. That is the calibration invariance a caller
 * actually relies on, and it is what tests/test_ccd.c asserts.
 *
 * WHAT A MAP FROM THIS IS NOT. The source paper implements no detection
 * threshold -- it is explicitly a proof of concept, with no ROC and no
 * false-alarm rate -- and its two trials use corner reflectors in open fields.
 * So a bright pixel here has no calibrated meaning on its own. The instrument
 * for supplying one already exists: rs_null_static() simulates a motionless
 * scene through the identical chain, and the map it produces is the floor this
 * one has to clear. A CCD map without that comparison is not evidence. See
 * docs/CCD-MICROMOTION-LOCATOR.md.
 */

#ifndef RESONARSAT_CCD_H
#define RESONARSAT_CCD_H

#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/subaperture.h"

/* Which maximal-invariant statistic to accumulate.
 *
 * Only the two-channel case is implemented. The three-channel statistics of the
 * source paper (its Eqs. 5 and 7-9, including a cubic solve for gamma_opt) are
 * not: G_1,2 is the simplest complete statistic and the only one whose
 * invariance can be checked in closed form, which is what makes it testable
 * rather than merely plausible. */
typedef enum {
    RS_CCD_G12 = 0   /* lambda_1 / lambda_2, two channels */
} rs_ccd_stat_t;

typedef struct {
    size_t win;          /* sliding window side in pixels; the paper uses 5 */
    rs_ccd_stat_t stat;

    /* Diagonal loading, as a fraction of the STACK'S OWN MEAN PIXEL POWER.
     *
     * NOT IN THE SOURCE PAPER, and needed for a reason the paper does not
     * discuss -- its scenes are real clutter, where no window is empty. Where a
     * window holds no scatterers, both covariances are noise, the ratio is free
     * to take any value, and the detector reports its brightest response over
     * bare ground. Measured on the synthetic two-target scene: 75612 over
     * background against 1.5 at the targets, before this was added.
     *
     * The floor must be ABSOLUTE, not a fraction of each window's own trace: a
     * trace-relative term vanishes exactly where it is needed. Derived from the
     * scene's mean power so it still scales with the data and preserves the
     * whole-stack scaling identity.
     *
     * Exposed rather than hard-coded so its effect can be measured; zero gives
     * the unregularised behaviour, which is worth looking at once. */
    double loading;
} rs_ccd_params_t;

/* An accumulated statistic map over the stack's image grid.
 *
 * 'map' is n_row x n_col in row-major order, on the same grid as the sub-look
 * images, so it can be overlaid on a focused image directly. Pixels within half
 * a window of the border are not computed and hold zero: a partial window would
 * give a covariance from fewer samples and a systematically different
 * statistic, which reads as a bright frame around the scene.
 *
 * 'n_triples' records how many sub-aperture triples contributed, which is
 * n_looks - 2. It is carried because the map is a MEAN over them, and comparing
 * two maps means knowing they averaged comparable numbers of terms. */
typedef struct {
    double *map;
    size_t  n_row, n_col;
    size_t  n_triples;
    rs_ccd_params_t params;
} rs_ccd_t;

/* Fill 'params' with the source paper's window and a conservative loading. */
void rs_ccd_params_default(rs_ccd_params_t *params);

/* Accumulate the change statistic over every sub-aperture triple in 'stack'.
 *
 * Needs at least three looks. The result is the MEAN of G over the triples
 * rather than the source paper's Eq. 10 sum divided by the channel count; the
 * two differ by a constant factor and the mean makes maps from different look
 * counts comparable, which the sum does not.
 *
 * Returns RS_ERR_ARG on a NULL argument, a stack of fewer than three looks, or a
 * window larger than the image; RS_ERR_ALLOC if the map cannot be sized. On
 * success '*out' owns its map and must be released with rs_ccd_free(). */
resonarsat_status_t rs_ccd_locate(const rs_subap_stack_t *stack,
                                  const rs_ccd_params_t *params,
                                  rs_ccd_t *out);

/* Release a map and everything it owns. */
void rs_ccd_free(rs_ccd_t *ccd);

#endif /* RESONARSAT_CCD_H */
