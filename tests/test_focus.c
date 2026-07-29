/* Backprojection: does a point target focus where it was put?
 *
 * This is the foundational ground-truth test. It depends on none of the
 * contested physics downstream -- it checks only that image formation is
 * geometrically correct -- so a failure here invalidates everything after it,
 * and a pass makes every later test interpretable. */

#include "resonarsat/focus.h"
#include "resonarsat/geom.h"
#include "rs_sim.h"
#include "rs_test.h"

/* Locate the brightest cell of an image. */
static void rs_peak_cell(const rs_slc_t *img, size_t *py, size_t *px, double *pmag)
{
    size_t best = 0;
    double best_mag = -1.0;
    for (size_t i = 0; i < img->n_az * img->n_rg; i++) {
        const double m = (double)cabsf(img->data[i]);
        if (m > best_mag) { best_mag = m; best = i; }
    }
    /* Row is along-track (azimuth), column is cross-track (range). */
    *py = best / img->n_rg;
    *px = best % img->n_rg;
    *pmag = best_mag;
}

int main(void)
{
    /* A short dwell keeps the test fast; the geometry is unchanged. */
    const double t_dwell = 1.0, prf = 800.0;
    const size_t n_rbin = 256;
    const double dr = 0.5;

    /* One target displaced 10 m along track from the scene centre. */
    const rs_sim_tgt_t tg[] = { { .x = 10.0, .y = 0.0, .z = 0.0, .rcs = 1.0 } };

    rs_cphd_t cphd;
    RS_CASE("simulate phase history");
    RS_CHECK_OK(rs_sim_scene(&cphd, tg, 1, t_dwell, prf, n_rbin, dr));

    rs_grid_t grid = { .origin = {0,0,0}, .n_x = 64, .n_y = 64,
                       .dx = 1.0, .dy = 1.0, .height = 0.0 };

    rs_slc_t img;
    RS_CASE("focus the full aperture");
    RS_CHECK_OK(rs_focus_full(&cphd, &grid, &img));

    /* The target must appear at its injected position. Grid cell (n/2 + 10) in
     * x corresponds to x = +10 m given the centred grid convention. */
    RS_CASE("point target focuses at its injected position");
    size_t py = 0, px = 0;
    double pmag = 0.0;
    rs_peak_cell(&img, &py, &px, &pmag);

    /* Along-track position comes from the ROW index. */
    const double x_found = ((double)py - 0.5 * (double)(grid.n_x - 1)) * grid.dx;
    RS_CHECK_NEAR(x_found, 10.0, 1.5);
    RS_CHECK(pmag > 0.0);

    /* Energy must concentrate: a correctly focused point is far brighter than
     * the mean. A defocused image -- the signature of a wrong phase sign --
     * spreads energy and fails this even though the peak may sit in the right
     * place. */
    RS_CASE("focused energy concentrates at the peak");
    double sum = 0.0;
    for (size_t i = 0; i < img.n_az * img.n_rg; i++) sum += (double)cabsf(img.data[i]);
    const double mean = sum / (double)(img.n_az * img.n_rg);
    RS_CHECK(pmag > 20.0 * mean);

    /* A sub-aperture must still focus the target at the same place, with less
     * resolution. This is what the sub-aperture stage relies on. */
    RS_CASE("a pulse sub-window focuses at the same position");
    rs_slc_t sub;
    RS_CHECK_OK(rs_slc_alloc(&sub, grid.n_x, grid.n_y));
    RS_CHECK_OK(rs_focus_backproject(&cphd, &grid, 0, cphd.n_pulse / 4, &sub));

    size_t sy = 0, sx = 0;
    double smag = 0.0;
    rs_peak_cell(&sub, &sy, &sx, &smag);
    const double x_sub = ((double)sy - 0.5 * (double)(grid.n_x - 1)) * grid.dx;
    RS_CHECK_NEAR(x_sub, 10.0, 3.0);

    RS_CASE("invalid pulse windows are refused");
    RS_CHECK_ERR(rs_focus_backproject(&cphd, &grid, 0, cphd.n_pulse + 1, &sub), RS_ERR_ARG);
    RS_CHECK_ERR(rs_focus_backproject(&cphd, &grid, cphd.n_pulse, 1, &sub), RS_ERR_ARG);
    RS_CHECK_ERR(rs_focus_backproject(&cphd, &grid, 0, 0, &sub), RS_ERR_ARG);

    rs_slc_free(&sub);
    rs_slc_free(&img);
    rs_cphd_free(&cphd);
    RS_TEST_END();
}
