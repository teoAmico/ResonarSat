/* The complex image type every stage of the pipeline flows through, plus the
 * orbit and Doppler metadata that travels with it. */

#ifndef RESONARSAT_SLC_H
#define RESONARSAT_SLC_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#include "resonarsat/geocode.h"
#include "resonarsat/resonarsat.h"

#define RS_DOPP_POLY_MAX 6   /* highest polynomial order we retain */
#define RS_ORBIT_MAX     64  /* state vectors kept per product */

/* A Doppler centroid or FM-rate polynomial in slant-range time.
 *
 * Evaluated as sum(coeff[i] * (tau - tau_ref)^i) for i in [0, n_coeff), with
 * 'tau' the two-way slant-range time in seconds. Products that carry only a
 * constant Doppler centroid set n_coeff = 1. */
typedef struct {
    double coeff[RS_DOPP_POLY_MAX];
    size_t n_coeff;
    double tau_ref;   /* reference slant-range time, s */
} rs_dopp_poly_t;

/* One orbit state vector: time in seconds since the product epoch, position
 * and velocity in metres and metres per second in an earth-fixed frame. */
typedef struct {
    double t;
    double pos[3];
    double vel[3];
} rs_state_vector_t;

/* Platform trajectory as a short table of state vectors, interpolatable. */
typedef struct {
    rs_state_vector_t sv[RS_ORBIT_MAX];
    size_t            n;
} rs_orbit_t;

/* A focused single-look complex image and everything needed to process it.
 *
 * Azimuth sampling deserves the comment it gets below. For Sentinel-1 IW the
 * transmit PRF (about 1451.6 Hz) is shared across the three IW sub-swaths and
 * is roughly three times the per-swath azimuth line rate (about 486 Hz for
 * IW2). Storing "the PRF" as the azimuth sampling rate scales every Doppler
 * axis in the pipeline by that factor, and the error hides itself: a Doppler
 * bandwidth measured against an axis built from the same wrong constant still
 * lands on its expected value. So the line period is the stored quantity, the
 * sampling rate is derived from it, and the transmit PRF is carried separately
 * where no signal-processing code will reach for it by accident. */
typedef struct {
    float complex *data;      /* row-major [azimuth][range], owned */
    size_t n_az, n_rg;

    double azimuth_time_interval;  /* s per azimuth line */
    double fs_az;                  /* Hz, == 1/azimuth_time_interval (derived) */
    double pulse_prf;              /* Hz, transmit PRF -- diagnostics only */

    double fc;                     /* radar centre frequency, Hz */
    double lambda;                 /* radar wavelength, m (derived from fc) */
    double rg_spacing_m;           /* slant-range sample spacing, m */
    double az_spacing_m;           /* azimuth sample spacing, m */
    double r0;                     /* slant range of first range sample, m */
    double t0;                     /* azimuth time of first line, s */
    double t_dwell;                /* target illumination time, s */
    double incidence;              /* incidence angle at scene centre, rad */
    double v_platform;             /* platform speed, m/s */

    rs_dopp_poly_t doppler;
    rs_orbit_t     orbit;

    /* Where this image sits on the earth, when the product says. Zeroed and
     * marked invalid for products that carry no plane definition. */
    rs_geo_plane_t plane;

    char source[64];               /* e.g. "UAVSAR", "SICD", for provenance */
} rs_slc_t;

/* Allocate an image of 'n_az' by 'n_rg' complex samples and zero it.
 *
 * Only the sample buffer and the dimensions are set; every metadata field is
 * left zeroed for the caller or reader to populate. Returns RS_ERR_ALLOC if the
 * buffer cannot be obtained, and RS_ERR_ARG if either dimension is zero or the
 * product would overflow. On failure 'img' is left untouched. */
resonarsat_status_t rs_slc_alloc(rs_slc_t *img, size_t n_az, size_t n_rg);

/* Release an image's sample buffer and zero its dimensions. Accepts an image
 * whose data pointer is already NULL, so it is safe on any error path, and
 * leaves the struct reusable for a subsequent rs_slc_alloc(). */
void rs_slc_free(rs_slc_t *img);

/* Set the derived fields that must never be parsed independently of the values
 * they come from: 'fs_az' from 'azimuth_time_interval', and 'lambda' from 'fc'.
 *
 * Readers call this once after populating the primary metadata. Keeping the
 * derivation in one place is what stops a reader from filling 'fs_az' straight
 * out of a product's "PRF" field, which is the specific mistake described in
 * the struct comment above.
 *
 * Returns RS_ERR_MISSING_META if either source field is still zero, since an
 * image with no azimuth timing or no carrier frequency cannot be processed and
 * failing here is far cheaper than discovering it three stages downstream. */
resonarsat_status_t rs_slc_finalise_metadata(rs_slc_t *img);

/* Check that an image's metadata is self-consistent and physically plausible,
 * writing a description of the first problem found via rs_set_error().
 *
 * Applies bounds a real spaceborne or airborne SAR product must satisfy: the
 * azimuth sampling rate lies in [1, 100000] Hz, the wavelength in [0.001, 1] m,
 * sample spacings in [0.01, 1000] m, and the incidence angle in (0, pi/2). It
 * additionally cross-checks the azimuth sampling rate against platform speed
 * and azimuth spacing when both are known, since those three cannot vary
 * independently, and a threefold disagreement there is the signature of the
 * PRF-for-line-rate substitution.
 *
 * This is advisory: it is called by the CLI's info path and by readers after
 * parsing, and it catches transcription errors rather than proving a product
 * good. Returns RS_OK when every check passes. */
resonarsat_status_t rs_slc_validate(const rs_slc_t *img);

/* Evaluate a Doppler polynomial at slant-range time 'tau' in seconds. An unset
 * polynomial evaluates to zero, which is correct for products shipping no
 * Doppler annotation -- the pipeline then estimates the centroid from data. */
double rs_dopp_poly_eval(const rs_dopp_poly_t *poly, double tau);

/* Interpolate the platform state at time 't' seconds into 'out', using Lagrange
 * interpolation over the four state vectors nearest 't'. Extrapolates rather
 * than failing just outside the tabulated span, since sub-aperture centre times
 * can fall marginally outside it. Returns RS_ERR_MISSING_META on an empty
 * orbit table. */
resonarsat_status_t rs_orbit_interp(const rs_orbit_t *orbit, double t,
                                    rs_state_vector_t *out);

/* Return a pointer to the first sample of azimuth line 'az'. No bounds check;
 * callers iterate over known dimensions. */
static inline float complex *rs_slc_row(rs_slc_t *img, size_t az)
{
    return img->data + az * img->n_rg;
}

#endif /* RESONARSAT_SLC_H */
