/* Cut a small patch out of a large scene.
 *
 * All development happens on patches of a couple of thousand pixels centred on
 * the structure of interest, never on multi-gigabyte scenes: backprojection and
 * per-window tracking both scale with area, and a full scene turns a two-minute
 * experiment into an overnight one. */

#include "resonarsat/readers.h"
#include "resonarsat/raster.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Copy a rectangular region out of one image into another, which the caller
 * must free. Bounds are clamped to the source rather than rejected, so a patch
 * requested near an edge yields the available part instead of an error. */
static resonarsat_status_t rs_crop(const rs_slc_t *src, size_t az0, size_t rg0,
                                   size_t n_az, size_t n_rg, rs_slc_t *dst)
{
    if (!src || !dst || !src->data) return RS_ERR_ARG;
    if (az0 >= src->n_az || rg0 >= src->n_rg) {
        rs_set_error("crop: origin (%zu,%zu) outside image %zux%zu",
                     az0, rg0, src->n_az, src->n_rg);
        return RS_ERR_ARG;
    }
    if (az0 + n_az > src->n_az) n_az = src->n_az - az0;
    if (rg0 + n_rg > src->n_rg) n_rg = src->n_rg - rg0;

    resonarsat_status_t st = rs_slc_alloc(dst, n_az, n_rg);
    if (st != RS_OK) return st;

    for (size_t a = 0; a < n_az; a++) {
        memcpy(dst->data + a * n_rg,
               src->data + (az0 + a) * src->n_rg + rg0,
               n_rg * sizeof *dst->data);
    }

    /* Carry the geometry across, adjusting only what the offset changes. */
    const size_t keep_az = dst->n_az, keep_rg = dst->n_rg;
    float complex *keep = dst->data;
    *dst = *src;
    dst->data = keep;
    dst->n_az = keep_az;
    dst->n_rg = keep_rg;
    dst->r0 = src->r0 + (double)rg0 * src->rg_spacing_m;
    dst->t0 = src->t0 + (double)az0 * src->azimuth_time_interval;
    dst->t_dwell = (double)n_az * src->azimuth_time_interval;
    return RS_OK;
}

/* Command-line entry point. */
int main(int argc, char **argv)
{
    if (argc < 7) {
        printf("usage: crop_slc SLC ANN AZ0 RG0 SIZE OUT.pgm\n");
        return 1;
    }

    rs_slc_t src;
    resonarsat_status_t st = rs_read_uavsar(argv[1], argv[2], &src);
    if (st != RS_OK) { rs_report_error("crop_slc", st); return 1; }

    const size_t az0 = (size_t)atol(argv[3]);
    const size_t rg0 = (size_t)atol(argv[4]);
    const size_t size = (size_t)atol(argv[5]);

    rs_slc_t dst;
    st = rs_crop(&src, az0, rg0, size, size, &dst);
    if (st != RS_OK) { rs_report_error("crop_slc", st); rs_slc_free(&src); return 1; }

    st = rs_raster_write_quicklook(&dst, argv[6], 40.0);
    if (st != RS_OK) rs_report_error("crop_slc", st);
    else printf("wrote %s: %zux%zu\n", argv[6], dst.n_az, dst.n_rg);

    rs_slc_free(&dst);
    rs_slc_free(&src);
    return st == RS_OK ? 0 : 1;
}
