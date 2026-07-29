/* Putting a pixel on the earth.
 *
 * Every stage before this works in image coordinates -- azimuth and range
 * samples, or metres on a scene-local plane. That is enough to process a product
 * and not enough to say anything about a place, which is the difference between
 * "there is a feature at pixel 1200" and "there is a feature at 46.87 N,
 * 97.26 W". It is also what makes two products of the same ground comparable:
 * a focused CPHD grid lies in the ground plane while a SICD grid lies in the
 * slant plane, and until both are expressed in a common frame, cropping one
 * against the other compares the projection rather than the scene.
 *
 * The frame used throughout is earth-centred earth-fixed (ECF) on WGS 84, which
 * is what both formats already carry, so nothing here has to guess a datum. */

#ifndef RESONARSAT_GEOCODE_H
#define RESONARSAT_GEOCODE_H

#include "resonarsat/resonarsat.h"

/* The plane a product's image coordinates live on.
 *
 * Both formats describe their image grid the same way: a reference point in ECF
 * plus two orthogonal unit vectors spanning the plane the samples are laid out
 * on. SICD calls them the scene centre point and the Row and Col direction
 * vectors; CPHD calls them the image area reference point and uIAX/uIAY. The
 * distinction that matters is which plane they span -- SICD's is the slant
 * plane, CPHD's the ground plane -- and carrying the vectors rather than an
 * assumed orientation is what lets that difference be handled instead of
 * ignored.
 *
 * 'u_x' and 'u_y' are unit vectors; 'origin' is the ECF position of image
 * coordinate (0,0). A zeroed struct means the product carried no plane
 * definition, which every consumer must treat as "cannot geocode" rather than
 * as a plane at the centre of the earth. */
typedef struct {
    double origin[3];   /* ECF metres */
    double u_x[3];      /* unit vector along the first image axis */
    double u_y[3];      /* unit vector along the second */
    /* Sensor position at aperture centre, and the height of the plane's own
     * reference point. Together these are what rs_geo_project_to_height() needs
     * to turn a slant-plane coordinate into a ground one, and keeping them with
     * the plane means a caller cannot geolocate without the means to do it
     * correctly close at hand. */
    double sensor[3];     /* ECF metres; zeroed when unknown */
    double sensor_vel[3]; /* ECF m/s, needed for the Doppler cone */
    double ref_hae;     /* height of the reference point above the ellipsoid, m */

    int    valid;       /* zero when the product carried no plane */
    int    is_slant;    /* non-zero when the plane is the slant plane, so that
                         * geolocating it without projection gives points on a
                         * tilted surface rather than on the ground */
} rs_geo_plane_t;

/* Convert an ECF position to geodetic latitude, longitude and height above the
 * WGS 84 ellipsoid.
 *
 * Uses Bowring's closed-form solution rather than iterating: for terrestrial
 * heights it is accurate to well under a millimetre, and a closed form cannot
 * fail to converge on a degenerate input the way a loop can. Latitude and
 * longitude are returned in degrees, height in metres.
 *
 * Returns RS_ERR_ARG on a NULL argument or a position at the earth's centre,
 * where longitude is undefined. */
resonarsat_status_t rs_geo_ecf_to_llh(const double ecf[3],
                                      double *lat_deg, double *lon_deg, double *hae_m);

/* Convert geodetic coordinates to ECF. The exact inverse of the above, provided
 * so that a round trip can be tested rather than assumed. */
resonarsat_status_t rs_geo_llh_to_ecf(double lat_deg, double lon_deg, double hae_m,
                                      double ecf[3]);

/* ECF position of the point 'x' metres along the plane's first axis and 'y'
 * along its second, measured from the plane origin.
 *
 * This is a planar approximation: it treats the image grid as flat, which is
 * what both formats' own grid definitions assume over a scene of a few
 * kilometres. Over a 5 km scene the earth's curvature displaces the corners by
 * about half a metre, comparable to a resolution cell, so the approximation is
 * consistent with the products rather than a shortcut past them. */
resonarsat_status_t rs_geo_plane_point(const rs_geo_plane_t *plane,
                                       double x_m, double y_m, double ecf[3]);

/* Geodetic position of a point on the plane, combining the two calls above.
 * Returns RS_ERR_MISSING_META when the plane is not valid, since a product
 * without a plane definition cannot be located and a fabricated answer would be
 * worse than none. */
resonarsat_status_t rs_geo_plane_llh(const rs_geo_plane_t *plane,
                                     double x_m, double y_m,
                                     double *lat_deg, double *lon_deg, double *hae_m);

/* The inverse of rs_geo_plane_llh(): where a known place falls on the plane, in
 * the metres that --offset speaks.
 *
 * WHY THIS EXISTS, which is worth stating because the arithmetic is three lines
 * and doing it by hand has already gone wrong on this project. A processing
 * grid is positioned by an offset in the product's own planar frame, and that
 * frame is the file's, not a compass: CPHD's uIAX/uIAY and SICD's row/column
 * vectors point where the collector chose. An offset derived by assuming the
 * axes are north and east, or azimuth and ground range, lands somewhere
 * plausible and wrong, and the resulting image looks like a perfectly ordinary
 * image OF THE WRONG GROUND -- there is nothing in it to say so. That is how a
 * run over the Giza plateau came to be centred 900 m off the Great Pyramid,
 * with its offset recorded in a manifest as though measured.
 *
 * Taking the axes from the product removes the assumption entirely. 'hae_m' is
 * the target's height above the ellipsoid; pass the plane's own 'ref_hae' when
 * the target sits on the reference surface, since a height error displaces the
 * result along the line of sight.
 *
 * Any component of the position out of the plane is dropped, which is the same
 * planar approximation rs_geo_plane_point() documents. Returns
 * RS_ERR_MISSING_META when the plane is not valid. */
resonarsat_status_t rs_geo_plane_offset(const rs_geo_plane_t *plane,
                                        double lat_deg, double lon_deg,
                                        double hae_m,
                                        double *x_m, double *y_m);

/* Project a point onto a surface of constant height above the ellipsoid, along
 * the line of sight from the sensor.
 *
 * This is the step that turns a plane coordinate into a ground coordinate, and
 * omitting it is a silent error rather than a loud one. A SICD image grid lies
 * in the SLANT plane, which is tilted from the ground by the grazing angle: on
 * a 4.1 km scene at 29.9 degrees grazing, the far edge of the grid sits about
 * 2 km higher than the near edge. Geolocating grid points directly returns
 * positions on that tilted plane -- correct as geometry, wrong as an answer to
 * "where on the ground is this" -- and the resulting latitudes and longitudes
 * look entirely plausible while being displaced by kilometres.
 *
 * The ray runs from 'sensor_ecf' through 'point_ecf' and is followed until it
 * reaches 'target_hae' metres above the ellipsoid.
 *
 * WHAT THIS IS NOT SUFFICIENT FOR. It cannot geolocate a slant-plane grid, and
 * the failure is silent. In such a grid the range axis runs essentially along
 * the line of sight -- on the reference product both sit about 30 degrees from
 * horizontal -- so following the ray back to a fixed height maps every range
 * position onto nearly the same ground point. Tried on a 4.1 km scene it
 * collapsed the whole range axis into about 10 metres while returning heights
 * uniformly equal to the reference height, which reads as a correct and
 * consistent answer and is nothing of the sort.
 *
 * A slant-plane grid needs range-Doppler geolocation: the intersection of the
 * range sphere about the sensor, the Doppler cone about its velocity vector, and
 * the earth surface. That is not implemented here, and until it is, SICD grids
 * are not geolocated rather than being geolocated wrongly.
 *
 * What this IS for: a ground-plane grid, such as a CPHD image area reference
 * surface, where the projection is well conditioned because the grid axes are
 * not parallel to the look direction.
 *
 * Returns RS_ERR_RANGE if the ray does not reach the target height. */
resonarsat_status_t rs_geo_project_to_height(const double sensor_ecf[3],
                                             const double point_ecf[3],
                                             double target_hae,
                                             double out_ecf[3]);

/* Locate a slant-plane grid point on the ground, by range and Doppler.
 *
 * A point on a slant-plane grid already carries the two quantities that fix a
 * target's position: its distance from the sensor, and its angle to the velocity
 * vector. What it does not carry is height -- it sits on a tilted plane rather
 * than on the terrain. So the range and the Doppler cone angle are read off the
 * plane point, and the ground point is the place with the SAME range and Doppler
 * at the wanted height.
 *
 * Geometrically: the constant-range surface is a sphere about the sensor, the
 * constant-Doppler surface is a cone about the velocity vector, and the two meet
 * in a circle. The answer is where that circle crosses the target height.
 *
 * The circle generally crosses twice, once either side of the ground track, and
 * choosing wrongly mirrors the scene across the flight path. The root nearest
 * the input point is taken, which keeps the side of track the product was
 * imaged on without needing to interpret a left/right flag.
 *
 * This replaces projecting along the line of sight, which is degenerate here: a
 * slant-plane range axis runs nearly parallel to the look direction, so ray
 * projection collapses the whole range extent onto a point while returning
 * heights that look perfectly correct.
 *
 * Returns RS_ERR_RANGE if the circle never reaches 'target_hae', which happens
 * for a geometry that does not see the ground at that height. */
resonarsat_status_t rs_geo_slant_to_ground(const double sensor_ecf[3],
                                           const double sensor_vel[3],
                                           const double slant_ecf[3],
                                           double target_hae,
                                           double out_ecf[3]);

#endif /* RESONARSAT_GEOCODE_H */
