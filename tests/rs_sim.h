/* Shared synthetic-scene helper for the test suite.
 *
 * Builds phase history for a straight-line spotlight collect over point
 * targets, some vibrating, so that several test binaries can share one
 * definition of ground truth rather than each rolling its own and drifting. */

#ifndef RS_SIM_H
#define RS_SIM_H

#include "resonarsat/focus.h"
#include "resonarsat/geom.h"

#include <math.h>
#include <stdlib.h>

/* A simulated point target. */
typedef struct {
    double x, y, z, rcs, vib_freq, vib_amp, vib_phase;
} rs_sim_tgt_t;

/* Generate range-compressed phase history for a target list.
 *
 * Geometry: platform along +x at 'v_platform', offset cross-track by
 * 'range_offset' and at height 'height', staring at the origin. The compressed
 * pulse is a Gaussian of width 'range_res' carrying the exact propagation
 * phase, which is all any stage downstream reads.
 *
 * Returns RS_OK, or an allocation failure from rs_cphd_alloc(). */
static resonarsat_status_t rs_sim_scene(rs_cphd_t *cphd,
                                        const rs_sim_tgt_t *tg, size_t n_tgt,
                                        double t_dwell, double prf,
                                        size_t n_rbin, double dr)
{
    const double fc = 9.6e9;
    const double v_platform = 7500.0;
    const double height = 500000.0;
    const double range_offset = 350000.0;
    const double range_res = 1.0;

    const size_t n_pulse = (size_t)(prf * t_dwell);
    resonarsat_status_t st = rs_cphd_alloc(cphd, n_pulse, n_rbin);
    if (st != RS_OK) return st;

    cphd->fc = fc;
    cphd->lambda = RS_C_LIGHT / fc;
    cphd->prf = prf;
    cphd->dr = dr;

    const double r_centre = sqrt(height * height + range_offset * range_offset);
    cphd->r_near = r_centre - 0.5 * (double)n_rbin * dr;

    const double k_phase = 4.0 * M_PI / cphd->lambda;
    const double sigma = range_res / 2.355;

    for (size_t i = 0; i < n_pulse; i++) {
        const double t = (double)i / prf;
        cphd->t[i] = t;
        cphd->pos[3 * i + 0] = v_platform * (t - 0.5 * t_dwell);
        cphd->pos[3 * i + 1] = range_offset;
        cphd->pos[3 * i + 2] = height;

        /* Motion compensation: the receive window follows the scene reference
         * point (the origin), so bin n_rbin/2 always sits on it. Without this
         * a long dwell would walk the target clean out of a narrow swath. */
        cphd->r_ref[i] = sqrt(cphd->pos[3 * i + 0] * cphd->pos[3 * i + 0]
                            + cphd->pos[3 * i + 1] * cphd->pos[3 * i + 1]
                            + cphd->pos[3 * i + 2] * cphd->pos[3 * i + 2]);

        float complex *row = cphd->signal + i * n_rbin;

        for (size_t g = 0; g < n_tgt; g++) {
            double dz = 0.0;
            if (tg[g].vib_freq > 0.0 && tg[g].vib_amp != 0.0) {
                dz = tg[g].vib_amp * sin(2.0 * M_PI * tg[g].vib_freq * t + tg[g].vib_phase);
            }
            const double dx = cphd->pos[3 * i + 0] - tg[g].x;
            const double dy = cphd->pos[3 * i + 1] - tg[g].y;
            const double dzz = cphd->pos[3 * i + 2] - (tg[g].z + dz);
            const double R = sqrt(dx * dx + dy * dy + dzz * dzz);

            const double fbin = (R - cphd->r_ref[i]) / dr + 0.5 * (double)n_rbin;
            if (fbin < 0.0 || fbin >= (double)n_rbin) continue;

            const long lo = (long)floor(fbin - 4.0 * sigma / dr);
            const long hi = (long)ceil(fbin + 4.0 * sigma / dr);
            const double ph = -k_phase * R;
            const double cr = cos(ph), ci = sin(ph);

            for (long b = lo; b <= hi; b++) {
                if (b < 0 || b >= (long)n_rbin) continue;
                const double d = ((double)b - fbin) * dr;
                const double env = tg[g].rcs * exp(-0.5 * (d * d) / (sigma * sigma));
                row[b] += (float)(env * cr) + (float)(env * ci) * I;
            }
        }
    }
    return RS_OK;
}

#endif /* RS_SIM_H */
