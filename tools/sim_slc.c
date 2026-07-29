/* Synthetic focused-image generator.
 *
 * Where sim_cphd produces phase history for the focusing stage to consume, this
 * produces an already-focused image directly, by placing analytic point spread
 * functions whose azimuth phase carries a known vibration. It exists so that
 * the sub-aperture and tracking stages can be tested without running
 * backprojection first, which keeps their unit tests fast and isolates a
 * failure to the stage under test. */

#include "resonarsat/geom.h"
#include "resonarsat/raster.h"
#include "resonarsat/slc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Fill an image with point targets, one of which oscillates in azimuth.
 *
 * Each target is a separable Gaussian in range and azimuth. The vibrating one
 * has its azimuth centre displaced sinusoidally as a function of the azimuth
 * coordinate, which stands in for the target moving during the aperture: after
 * sub-aperture decomposition each look sees it at a different point in its
 * cycle, which is exactly the signal the tracker measures. */
static void rs_sim_fill(rs_slc_t *img, double vib_freq, double vib_amp_px)
{
    const double sigma = 1.5;
    struct { double az, rg, amp, f, a; } tg[] = {
        { 0.25, 0.30, 1.0, 0.0,      0.0        },
        { 0.50, 0.50, 1.0, vib_freq, vib_amp_px },
        { 0.75, 0.70, 0.8, 0.0,      0.0        },
    };

    for (size_t t = 0; t < sizeof tg / sizeof tg[0]; t++) {
        const double cz = tg[t].az * (double)img->n_az;
        const double cr = tg[t].rg * (double)img->n_rg;

        for (size_t a = 0; a < img->n_az; a++) {
            /* Azimuth position is time within the aperture. */
            const double t_s = (double)a * img->azimuth_time_interval;
            const double shift = (tg[t].f > 0.0)
                               ? tg[t].a * sin(2.0 * M_PI * tg[t].f * t_s)
                               : 0.0;
            const double da = ((double)a - (cz + shift)) / sigma;
            if (fabs(da) > 4.0) continue;

            for (size_t r = 0; r < img->n_rg; r++) {
                const double dr = ((double)r - cr) / sigma;
                if (fabs(dr) > 4.0) continue;
                const double env = tg[t].amp * exp(-0.5 * (da * da + dr * dr));
                img->data[a * img->n_rg + r] += (float)env;
            }
        }
    }
}

/* Command-line entry point. */
int main(int argc, char **argv)
{
    const char *out = (argc > 1) ? argv[1] : "sim_slc.pgm";
    const size_t n = (argc > 2) ? (size_t)atol(argv[2]) : 256;
    const double f = (argc > 3) ? atof(argv[3]) : 2.0;
    const double a = (argc > 4) ? atof(argv[4]) : 2.0;

    rs_slc_t img;
    resonarsat_status_t st = rs_slc_alloc(&img, n, n);
    if (st != RS_OK) { rs_report_error("sim_slc", st); return 1; }

    img.fc = 9.6e9;
    img.azimuth_time_interval = 10.0 / (double)n;
    img.t_dwell = 10.0;
    img.v_platform = 7500.0;
    img.r0 = 500000.0;
    img.az_spacing_m = img.rg_spacing_m = 1.0;
    rs_slc_finalise_metadata(&img);

    rs_sim_fill(&img, f, a);

    st = rs_raster_write_quicklook(&img, out, 40.0);
    if (st != RS_OK) rs_report_error("sim_slc", st);
    else printf("wrote %s: %zux%zu, vibrating target %g Hz at %g px\n", out, n, n, f, a);

    rs_slc_free(&img);
    return st == RS_OK ? 0 : 1;
}
