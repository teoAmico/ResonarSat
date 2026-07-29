/* Geodetic conversion and image-plane geolocation. See geocode.h. */

#include "resonarsat/geocode.h"

#include <math.h>
#include <string.h>

/* WGS 84 defining parameters. The semi-major axis and the flattening are the
 * two the standard defines; everything else here is derived from them, so a
 * transcription error in a derived constant cannot silently disagree with the
 * datum. */
#define RS_WGS84_A  6378137.0
#define RS_WGS84_F  (1.0 / 298.257223563)

/* Latitude and height from an ECF position. See geocode.h on the choice of a
 * closed form. */
resonarsat_status_t rs_geo_ecf_to_llh(const double ecf[3],
                                      double *lat_deg, double *lon_deg, double *hae_m)
{
    if (!ecf) return RS_ERR_ARG;

    const double a  = RS_WGS84_A;
    const double f  = RS_WGS84_F;
    const double b  = a * (1.0 - f);
    const double e2 = f * (2.0 - f);                 /* first eccentricity squared */
    const double ep2 = (a * a - b * b) / (b * b);    /* second, for Bowring */

    const double x = ecf[0], y = ecf[1], z = ecf[2];
    const double p = sqrt(x * x + y * y);

    if (p < 1.0 && fabs(z) < 1.0) {
        /* At the centre of the earth longitude has no value and latitude no
         * meaning. Returning zeros here would look like a position off West
         * Africa, which is the classic way a missing coordinate is mistaken for
         * a real one. */
        return RS_ERR_ARG;
    }

    /* Bowring's auxiliary angle, then latitude in one step. */
    const double theta = atan2(z * a, p * b);
    const double st = sin(theta), ct = cos(theta);
    const double lat = atan2(z + ep2 * b * st * st * st,
                             p - e2  * a * ct * ct * ct);
    const double lon = atan2(y, x);

    const double sl = sin(lat);
    const double N = a / sqrt(1.0 - e2 * sl * sl);

    double h;
    if (fabs(cos(lat)) > 1.0e-9) {
        h = p / cos(lat) - N;
    } else {
        /* Near the poles p tends to zero and p/cos(lat) loses all precision;
         * the polar form uses z instead. */
        h = fabs(z) - b;
    }

    if (lat_deg) *lat_deg = lat * 180.0 / M_PI;
    if (lon_deg) *lon_deg = lon * 180.0 / M_PI;
    if (hae_m)   *hae_m   = h;
    return RS_OK;
}

/* Geodetic to ECF. See geocode.h. */
resonarsat_status_t rs_geo_llh_to_ecf(double lat_deg, double lon_deg, double hae_m,
                                      double ecf[3])
{
    if (!ecf) return RS_ERR_ARG;
    if (lat_deg < -90.0 || lat_deg > 90.0) return RS_ERR_RANGE;

    const double a  = RS_WGS84_A;
    const double f  = RS_WGS84_F;
    const double e2 = f * (2.0 - f);

    const double lat = lat_deg * M_PI / 180.0;
    const double lon = lon_deg * M_PI / 180.0;
    const double sl = sin(lat), cl = cos(lat);
    const double N = a / sqrt(1.0 - e2 * sl * sl);

    ecf[0] = (N + hae_m) * cl * cos(lon);
    ecf[1] = (N + hae_m) * cl * sin(lon);
    ecf[2] = (N * (1.0 - e2) + hae_m) * sl;
    return RS_OK;
}

/* Plane coordinate to ECF. See geocode.h. */
resonarsat_status_t rs_geo_plane_point(const rs_geo_plane_t *plane,
                                       double x_m, double y_m, double ecf[3])
{
    if (!plane || !ecf) return RS_ERR_ARG;
    if (!plane->valid) return RS_ERR_MISSING_META;

    for (int i = 0; i < 3; i++) {
        ecf[i] = plane->origin[i] + x_m * plane->u_x[i] + y_m * plane->u_y[i];
    }
    return RS_OK;
}

/* Plane coordinate to geodetic. See geocode.h. */
resonarsat_status_t rs_geo_plane_llh(const rs_geo_plane_t *plane,
                                     double x_m, double y_m,
                                     double *lat_deg, double *lon_deg, double *hae_m)
{
    double ecf[3];
    const resonarsat_status_t st = rs_geo_plane_point(plane, x_m, y_m, ecf);
    if (st != RS_OK) return st;
    return rs_geo_ecf_to_llh(ecf, lat_deg, lon_deg, hae_m);
}

/* Where a geodetic position falls on the plane. See the header for why this is
 * a function rather than three lines at each call site. */
resonarsat_status_t rs_geo_plane_offset(const rs_geo_plane_t *plane,
                                        double lat_deg, double lon_deg,
                                        double hae_m,
                                        double *x_m, double *y_m)
{
    if (!plane) return RS_ERR_ARG;
    if (!plane->valid) return RS_ERR_MISSING_META;

    double ecf[3];
    const resonarsat_status_t st = rs_geo_llh_to_ecf(lat_deg, lon_deg, hae_m, ecf);
    if (st != RS_OK) return st;

    /* The axes are orthonormal, so projecting the displacement onto each is the
     * whole inverse; no solve is needed. */
    const double d[3] = { ecf[0] - plane->origin[0],
                          ecf[1] - plane->origin[1],
                          ecf[2] - plane->origin[2] };
    if (x_m) *x_m = d[0] * plane->u_x[0] + d[1] * plane->u_x[1] + d[2] * plane->u_x[2];
    if (y_m) *y_m = d[0] * plane->u_y[0] + d[1] * plane->u_y[1] + d[2] * plane->u_y[2];
    return RS_OK;
}

/* Follow the look ray to a surface of constant height. See geocode.h. */
resonarsat_status_t rs_geo_project_to_height(const double sensor_ecf[3],
                                             const double point_ecf[3],
                                             double target_hae,
                                             double out_ecf[3])
{
    if (!sensor_ecf || !point_ecf || !out_ecf) return RS_ERR_ARG;

    double dir[3];
    double norm = 0.0;
    for (int i = 0; i < 3; i++) {
        dir[i] = point_ecf[i] - sensor_ecf[i];
        norm += dir[i] * dir[i];
    }
    norm = sqrt(norm);
    if (norm < 1.0) return RS_ERR_ARG;
    for (int i = 0; i < 3; i++) dir[i] /= norm;

    /* Height as a function of distance along the ray, sampled at the point
     * itself and well beyond it, so the bracket spans the target. */
    double t_lo = 0.0, t_hi = 2.0 * norm;
    double p[3], h_lo = 0.0, h_hi = 0.0;

    for (int i = 0; i < 3; i++) p[i] = sensor_ecf[i] + t_lo * dir[i];
    if (rs_geo_ecf_to_llh(p, NULL, NULL, &h_lo) != RS_OK) return RS_ERR_ARG;
    for (int i = 0; i < 3; i++) p[i] = sensor_ecf[i] + t_hi * dir[i];
    if (rs_geo_ecf_to_llh(p, NULL, NULL, &h_hi) != RS_OK) return RS_ERR_ARG;

    if ((h_lo - target_hae) * (h_hi - target_hae) > 0.0) {
        rs_set_error("geocode: the look ray spans %.1f to %.1f m and never reaches %.1f m",
                     h_lo, h_hi, target_hae);
        return RS_ERR_RANGE;
    }

    for (int it = 0; it < 60; it++) {
        const double t_mid = 0.5 * (t_lo + t_hi);
        double h_mid;
        for (int i = 0; i < 3; i++) p[i] = sensor_ecf[i] + t_mid * dir[i];
        if (rs_geo_ecf_to_llh(p, NULL, NULL, &h_mid) != RS_OK) return RS_ERR_ARG;
        if ((h_lo - target_hae) * (h_mid - target_hae) <= 0.0) {
            t_hi = t_mid;
        } else {
            t_lo = t_mid; h_lo = h_mid;
        }
    }
    const double t = 0.5 * (t_lo + t_hi);
    for (int i = 0; i < 3; i++) out_ecf[i] = sensor_ecf[i] + t * dir[i];
    return RS_OK;
}

/* Height above the ellipsoid of a point on the range-Doppler circle at angle
 * 'phi', used as the function whose root is the ground position. */
static int rs_geo_circle_height(const double centre[3], double radius,
                                const double e1[3], const double e2[3],
                                double phi, double *h, double p_out[3])
{
    const double c = cos(phi), s = sin(phi);
    double p[3];
    for (int i = 0; i < 3; i++) p[i] = centre[i] + radius * (c * e1[i] + s * e2[i]);
    if (p_out) memcpy(p_out, p, 3 * sizeof *p);
    return rs_geo_ecf_to_llh(p, NULL, NULL, h) == RS_OK;
}

/* Range-Doppler geolocation. See geocode.h. */
resonarsat_status_t rs_geo_slant_to_ground(const double sensor_ecf[3],
                                           const double sensor_vel[3],
                                           const double slant_ecf[3],
                                           double target_hae,
                                           double out_ecf[3])
{
    if (!sensor_ecf || !sensor_vel || !slant_ecf || !out_ecf) return RS_ERR_ARG;

    double uv[3], vn = 0.0;
    for (int i = 0; i < 3; i++) vn += sensor_vel[i] * sensor_vel[i];
    vn = sqrt(vn);
    if (vn < 1.0) {
        rs_set_error("geocode: sensor velocity is zero; Doppler is undefined");
        return RS_ERR_ARG;
    }
    for (int i = 0; i < 3; i++) uv[i] = sensor_vel[i] / vn;

    /* Range and Doppler, read off the slant-plane point. */
    double d[3], R = 0.0, along = 0.0;
    for (int i = 0; i < 3; i++) { d[i] = slant_ecf[i] - sensor_ecf[i]; R += d[i] * d[i]; }
    R = sqrt(R);
    if (R < 1.0) return RS_ERR_ARG;
    for (int i = 0; i < 3; i++) along += d[i] * uv[i];

    /* The circle where the range sphere meets the Doppler cone. */
    double centre[3], e1[3], e1n = 0.0;
    for (int i = 0; i < 3; i++) {
        centre[i] = sensor_ecf[i] + along * uv[i];
        e1[i] = d[i] - along * uv[i];
        e1n += e1[i] * e1[i];
    }
    e1n = sqrt(e1n);
    if (e1n < 1.0) {
        rs_set_error("geocode: the look direction is along the velocity vector; "
                     "the range-Doppler circle degenerates to a point");
        return RS_ERR_RANGE;
    }
    const double radius = e1n;
    for (int i = 0; i < 3; i++) e1[i] /= e1n;

    double e2[3];
    e2[0] = uv[1] * e1[2] - uv[2] * e1[1];
    e2[1] = uv[2] * e1[0] - uv[0] * e1[2];
    e2[2] = uv[0] * e1[1] - uv[1] * e1[0];

    /* Scan for sign changes of (height - target), then bisect the one nearest
     * phi = 0, which is the input point and therefore the correct side of the
     * ground track. */
    const int N = 720;
    double best_lo = 0.0, best_hi = 0.0;
    int found = 0;
    double prev_h;
    if (!rs_geo_circle_height(centre, radius, e1, e2, -M_PI, &prev_h, NULL)) return RS_ERR_ARG;

    for (int k = 1; k <= N; k++) {
        const double phi = -M_PI + 2.0 * M_PI * (double)k / (double)N;
        double h;
        if (!rs_geo_circle_height(centre, radius, e1, e2, phi, &h, NULL)) return RS_ERR_ARG;
        if ((prev_h - target_hae) * (h - target_hae) <= 0.0) {
            const double lo = -M_PI + 2.0 * M_PI * (double)(k - 1) / (double)N;
            if (!found || fabs(0.5 * (lo + phi)) < fabs(0.5 * (best_lo + best_hi))) {
                best_lo = lo; best_hi = phi; found = 1;
            }
        }
        prev_h = h;
    }
    if (!found) {
        rs_set_error("geocode: the range-Doppler circle never reaches %.1f m", target_hae);
        return RS_ERR_RANGE;
    }

    double h_lo;
    (void)rs_geo_circle_height(centre, radius, e1, e2, best_lo, &h_lo, NULL);
    for (int it = 0; it < 60; it++) {
        const double mid = 0.5 * (best_lo + best_hi);
        double h_mid;
        if (!rs_geo_circle_height(centre, radius, e1, e2, mid, &h_mid, NULL)) return RS_ERR_ARG;
        if ((h_lo - target_hae) * (h_mid - target_hae) <= 0.0) {
            best_hi = mid;
        } else {
            best_lo = mid; h_lo = h_mid;
        }
    }
    double h_final;
    (void)rs_geo_circle_height(centre, radius, e1, e2,
                               0.5 * (best_lo + best_hi), &h_final, out_ecf);
    return RS_OK;
}
