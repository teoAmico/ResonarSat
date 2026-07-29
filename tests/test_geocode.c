/* Geodetic conversion and range-Doppler geolocation.
 *
 * The conversions are checked against a point whose position is stated twice,
 * independently, by a real product: the reference point of an Umbra collect
 * appears in its metadata both as an earth-centred vector and as a latitude and
 * longitude. Agreement between them is a real external check rather than a round
 * trip through our own arithmetic, which would pass even if the datum were
 * wrong.
 *
 * The geolocation is checked by its invariants instead, because there is no
 * published answer to compare against. A correct range-Doppler solution moves a
 * point onto the terrain WITHOUT changing either quantity that located it: the
 * distance to the sensor and the angle to the velocity vector must both survive.
 * That is a stronger test than it sounds, and it is exactly what the previous
 * attempt failed -- projecting along the line of sight preserved neither, and
 * collapsed a 4.8 km range extent to ten metres while returning heights that
 * looked perfectly correct. */

#include "resonarsat/geocode.h"
#include "rs_test.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

int main(void)
{
    /* The image area reference point of a real collect, as its metadata gives
     * it in both forms. */
    const double ref_ecf[3] = { -552223.3276367188, -4332875.0, 4632556.640625 };
    const double ref_lat = 46.876176372949395;
    const double ref_lon = -97.26316876232174;
    const double ref_hae = 261.3063053144142;

    RS_CASE("ECF converts to the geodetic position the product states");
    {
        double lat, lon, hae;
        RS_CHECK_OK(rs_geo_ecf_to_llh(ref_ecf, &lat, &lon, &hae));
        printf("    %.9f, %.9f, %.4f m\n", lat, lon, hae);
        RS_CHECK_NEAR(lat, ref_lat, 1e-9);
        RS_CHECK_NEAR(lon, ref_lon, 1e-9);
        RS_CHECK_NEAR(hae, ref_hae, 1e-3);
    }

    RS_CASE("the conversion round trips");
    {
        double ecf[3];
        RS_CHECK_OK(rs_geo_llh_to_ecf(ref_lat, ref_lon, ref_hae, ecf));
        const double e = sqrt((ecf[0]-ref_ecf[0])*(ecf[0]-ref_ecf[0]) +
                              (ecf[1]-ref_ecf[1])*(ecf[1]-ref_ecf[1]) +
                              (ecf[2]-ref_ecf[2])*(ecf[2]-ref_ecf[2]));
        printf("    round trip error %.6g m\n", e);
        RS_CHECK(e < 1e-6);
    }

    RS_CASE("a position at the earth's centre is refused, not answered");
    {
        /* Zeros here would convert to a plausible-looking point in the Gulf of
         * Guinea, which is the classic way an unset coordinate passes for a
         * real one. */
        const double centre[3] = { 0.0, 0.0, 0.0 };
        double lat, lon, hae;
        RS_CHECK_ERR(rs_geo_ecf_to_llh(centre, &lat, &lon, &hae), RS_ERR_ARG);
        RS_CHECK_ERR(rs_geo_ecf_to_llh(NULL, &lat, &lon, &hae), RS_ERR_ARG);
    }

    RS_CASE("an invalid plane cannot be geolocated");
    {
        rs_geo_plane_t p;
        memset(&p, 0, sizeof p);
        double lat, lon, hae;
        RS_CHECK_ERR(rs_geo_plane_llh(&p, 0.0, 0.0, &lat, &lon, &hae),
                     RS_ERR_MISSING_META);
        double x, y;
        RS_CHECK_ERR(rs_geo_plane_offset(&p, 0.0, 0.0, 0.0, &x, &y),
                     RS_ERR_MISSING_META);
    }

    /* The inverse is what positions a processing grid on a named place, so it
     * is checked against the forward direction rather than against a table of
     * expected numbers: whatever rs_geo_plane_llh() says is at (x,y) must map
     * back to (x,y). A one-sided test would pass just as happily with the two
     * axes swapped, which is the error that actually occurs and the one that
     * puts a grid on the wrong ground. */
    RS_CASE("plane offset inverts plane geolocation, including across the axes");
    {
        /* A plane built the way a reader builds one: origin at the reference
         * point, two orthonormal axes spanning the local horizontal. They are
         * deliberately NOT north/east -- the real uIAX/uIAY point wherever the
         * collector chose, and an inverse that only works on axis-aligned
         * frames would pass a north/east test and fail in the field. */
        double up[3] = { ref_ecf[0], ref_ecf[1], ref_ecf[2] };
        const double un = sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
        for (int i = 0; i < 3; i++) up[i] /= un;

        /* An arbitrary vector, orthogonalised against up, then rotated. */
        double ax[3] = { 0.37, -0.81, 0.45 };
        double dot = ax[0]*up[0] + ax[1]*up[1] + ax[2]*up[2];
        for (int i = 0; i < 3; i++) ax[i] -= dot * up[i];
        double an = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
        for (int i = 0; i < 3; i++) ax[i] /= an;

        rs_geo_plane_t plane;
        memset(&plane, 0, sizeof plane);
        for (int i = 0; i < 3; i++) {
            plane.origin[i] = ref_ecf[i];
            plane.u_x[i] = ax[i];
        }
        /* u_y = up x u_x, so the two axes are orthonormal and in the plane. */
        plane.u_y[0] = up[1]*ax[2] - up[2]*ax[1];
        plane.u_y[1] = up[2]*ax[0] - up[0]*ax[2];
        plane.u_y[2] = up[0]*ax[1] - up[1]*ax[0];
        plane.ref_hae = ref_hae;
        plane.valid = 1;

        const double probes[][2] = {
            {    0.0,    0.0 },
            {  500.0,    0.0 },      /* pure u_x -- catches a swap */
            {    0.0,  500.0 },      /* pure u_y */
            { -152.0, -552.0 },      /* the Giza case that motivated this */
            { 2000.0, -3000.0 },
        };
        for (size_t i = 0; i < sizeof probes / sizeof probes[0]; i++) {
            const double x_in = probes[i][0], y_in = probes[i][1];
            double lat, lon, hae;
            RS_CHECK_OK(rs_geo_plane_llh(&plane, x_in, y_in, &lat, &lon, &hae));

            double x_out, y_out;
            RS_CHECK_OK(rs_geo_plane_offset(&plane, lat, lon, hae, &x_out, &y_out));
            printf("    (%+8.1f, %+8.1f) -> %.6f, %.6f -> (%+8.1f, %+8.1f)\n",
                   x_in, y_in, lat, lon, x_out, y_out);
            RS_CHECK_NEAR(x_out, x_in, 1e-3);
            RS_CHECK_NEAR(y_out, y_in, 1e-3);
        }
    }

    /* ------------------------------------------------------------------
     * Range-Doppler geolocation, checked by what it must preserve.
     * ------------------------------------------------------------------ */
    RS_CASE("range-Doppler preserves range and Doppler while reaching the ground");
    {
        /* A sensor above the reference point, looking down at an angle, with a
         * slant point deliberately placed off the terrain. */
        const double sensor[3] = { 210216.5459579905, -4613027.927174058,
                                   5102894.483334189 };
        const double vel[3]    = { -2060.9627780307974, 5456.22495774894,
                                   5002.713301654648 };

        for (int k = -2; k <= 2; k++) {
            /* Points spread along the look direction, all off the ground. */
            double slant[3];
            for (int i = 0; i < 3; i++) {
                slant[i] = ref_ecf[i] + (double)k * 800.0 *
                           ((i == 0) ? -0.9868703498505056 :
                            (i == 1) ? -0.04207154083997011 : -0.15593876410275698);
            }

            double ground[3];
            RS_CHECK_OK(rs_geo_slant_to_ground(sensor, vel, slant, ref_hae, ground));

            /* Range to the sensor must be unchanged. */
            double r_in = 0.0, r_out = 0.0, dot_in = 0.0, dot_out = 0.0, vn = 0.0;
            for (int i = 0; i < 3; i++) {
                const double a = slant[i]  - sensor[i];
                const double b = ground[i] - sensor[i];
                r_in  += a * a;
                r_out += b * b;
                dot_in  += a * vel[i];
                dot_out += b * vel[i];
                vn += vel[i] * vel[i];
            }
            r_in = sqrt(r_in); r_out = sqrt(r_out); vn = sqrt(vn);
            const double cone_in  = acos(dot_in  / (r_in  * vn)) * 180.0 / M_PI;
            const double cone_out = acos(dot_out / (r_out * vn)) * 180.0 / M_PI;

            double h;
            RS_CHECK_OK(rs_geo_ecf_to_llh(ground, NULL, NULL, &h));

            printf("    offset %+5d m: range %.4f -> %.4f m, cone %.6f -> %.6f deg, "
                   "height %.3f m\n",
                   k * 800, r_in, r_out, cone_in, cone_out, h);

            RS_CHECK_NEAR(r_out, r_in, 1e-3);          /* range preserved */
            RS_CHECK_NEAR(cone_out, cone_in, 1e-6);    /* Doppler preserved */
            RS_CHECK_NEAR(h, ref_hae, 1e-3);           /* and it reached the ground */
        }
    }

    RS_CASE("a degenerate geometry is refused");
    {
        const double sensor[3] = { 7000000.0, 0.0, 0.0 };
        const double vel[3]    = { 0.0, 0.0, 0.0 };     /* no velocity, no Doppler */
        double out[3];
        RS_CHECK_ERR(rs_geo_slant_to_ground(sensor, vel, ref_ecf, 0.0, out), RS_ERR_ARG);
    }

    RS_TEST_END();
}
