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

    /* THE CLAIM --no-optimize MAKES ABOUT BACKPROJECTION, TESTED RATHER THAN
     * ASSERTED.
     *
     * rs_focus_opts_t states that the single-threaded mode produces bitwise
     * identical samples, because the parallel loop is over cells and each cell
     * accumulates privately over pulses in chronological order. An untested
     * statement of that kind is worth nothing -- it is exactly the sort of claim
     * that stays true until someone adds a reduction clause -- so it is checked
     * to the bit, not to a tolerance. A tolerance would pass on a build that HAD
     * acquired threading drift, which is the one thing this must detect.
     *
     * If this ever fails, the flag has become load-bearing and rs_focus_opts_t's
     * documentation is wrong and must be corrected before the flag is cited.
     *
     * Note what it does NOT prove: with OpenMP absent, or with one core, both
     * calls run the same code and the comparison is vacuous. It is still worth
     * having -- it is not vacuous on the machines that matter, and CMakeLists.txt
     * prints whether OpenMP was found. */
    RS_CASE("single-threaded focusing is bitwise identical to threaded");
    {
        rs_slc_t threaded, serial;
        RS_CHECK_OK(rs_slc_alloc(&threaded, grid.n_x, grid.n_y));
        RS_CHECK_OK(rs_slc_alloc(&serial, grid.n_x, grid.n_y));

        const rs_focus_opts_t par = { .single_thread = 0 };
        const rs_focus_opts_t seq = { .single_thread = 1 };
        RS_CHECK_OK(rs_focus_backproject_opts(&cphd, &grid, 0, cphd.n_pulse, &par,
                                              &threaded));
        RS_CHECK_OK(rs_focus_backproject_opts(&cphd, &grid, 0, cphd.n_pulse, &seq,
                                              &serial));

        size_t differing = 0;
        for (size_t i = 0; i < threaded.n_az * threaded.n_rg; i++) {
            if (crealf(threaded.data[i]) != crealf(serial.data[i]) ||
                cimagf(threaded.data[i]) != cimagf(serial.data[i])) {
                differing++;
            }
        }
        RS_CHECK(differing == 0);

        /* A NULL 'opts' must mean the default, so the wrapper cannot quietly
         * change behaviour for its existing callers. */
        rs_slc_t defaulted;
        RS_CHECK_OK(rs_slc_alloc(&defaulted, grid.n_x, grid.n_y));
        RS_CHECK_OK(rs_focus_backproject_opts(&cphd, &grid, 0, cphd.n_pulse, NULL,
                                              &defaulted));
        size_t default_diff = 0;
        for (size_t i = 0; i < threaded.n_az * threaded.n_rg; i++) {
            if (crealf(defaulted.data[i]) != crealf(threaded.data[i]) ||
                cimagf(defaulted.data[i]) != cimagf(threaded.data[i])) {
                default_diff++;
            }
        }
        RS_CHECK(default_diff == 0);

        rs_slc_free(&defaulted);
        rs_slc_free(&serial);
        rs_slc_free(&threaded);
    }

    /* The range interpolator, both kernels, on the same phase history.
     *
     * MEASURED ON REAL DATA FIRST, and the numbers are recorded here because a
     * synthetic fixture cannot reproduce them. Focusing the same 256x256 patch
     * of the Giza collect at 0.5 m cells, changing nothing but the kernel:
     *
     *   mean amplitude   +1.43 dB      peak amplitude  +1.06 dB
     *   median           +1.43 dB      top 20 targets  +1.32 dB
     *   total power      +1.42 dB      contrast (sd/mean)  0.5382 -> 0.5378
     *
     * So linear interpolation was costing about 1.4 dB of amplitude uniformly,
     * which is real and modest. What it was NOT doing is generating the
     * cross-range artefacts NGA's reference warns about: image contrast is
     * unchanged to four significant figures, so the correction is a gain, not a
     * cleanup. That distinction is why the wider kernel is an option rather
     * than the default -- see rs_focus_opts_t.
     *
     * WHAT THIS FIXTURE CANNOT CHECK, and a trap worth recording. The
     * simulator deposits a GAUSSIAN range response of sigma = range_res/2.355,
     * sampled at 'dr'. At the defaults that is 0.425 m sampled every 0.5 m --
     * 0.85 samples per sigma -- so the fixture's range signal is ALIASED. Sinc
     * interpolation assumes a band-limited signal and faithfully reconstructs
     * whatever is there, aliasing included, where linear interpolation
     * incidentally low-passes it. On this fixture the wider kernel therefore
     * raises total energy while LOWERING the peak, at 256 cells: peak 265.9,
     * 228.3, 161.2 for 0, 4 and 8 taps against total power 1.32e6, 1.44e6,
     * 1.43e6. On real range-compressed data, band-limited by the chirp, every
     * statistic including the peak rose instead.
     *
     * So an amplitude assertion here would pin the simulator's aliasing rather
     * than the interpolator's quality, and would have to be inverted against
     * what real data does. What the fixture does support is agreement about
     * WHERE the energy is, which any interpolator must give and a bug would
     * break. */
    RS_CASE("the wider range kernel agrees with linear on target position");
    {
        const rs_sim_tgt_t tgt_one = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
        rs_cphd_t c2;
        RS_CHECK_OK(rs_sim_scene(&c2, &tgt_one, 1, 4.0, 400.0, 256, 0.5));

        rs_grid_t g = { .n_x = 64, .n_y = 64, .dx = 0.5, .dy = 0.5, .height = 0.0 };
        g.origin[0] = g.origin[1] = g.origin[2] = 0.0;

        rs_slc_t lin, sinc;
        RS_CHECK_OK(rs_slc_alloc(&lin, g.n_x, g.n_y));
        RS_CHECK_OK(rs_slc_alloc(&sinc, g.n_x, g.n_y));

        const rs_focus_opts_t o_lin = { .single_thread = 0, .range_taps = 0 };
        const rs_focus_opts_t o_sinc = { .single_thread = 0, .range_taps = 8 };
        RS_CHECK_OK(rs_focus_backproject_opts(&c2, &g, 0, c2.n_pulse, &o_lin, &lin));
        RS_CHECK_OK(rs_focus_backproject_opts(&c2, &g, 0, c2.n_pulse, &o_sinc, &sinc));

        size_t pl = 0, ps = 0;
        double al = 0.0, as = 0.0;
        for (size_t i = 0; i < g.n_x * g.n_y; i++) {
            const double ml = cabs(lin.data[i]), ms = cabs(sinc.data[i]);
            if (ml > al) { al = ml; pl = i; }
            if (ms > as) { as = ms; ps = i; }
        }
        double el = 0.0, es = 0.0;
        for (size_t i = 0; i < g.n_x * g.n_y; i++) {
            const double ml = cabs(lin.data[i]), ms = cabs(sinc.data[i]);
            el += ml * ml;
            es += ms * ms;
        }
        printf("    peak cell: linear %zu, sinc %zu; amplitude %.4g vs %.4g; "
               "total power %+.2f dB\n", pl, ps, al, as,
               10.0 * log10(es / el));

        /* Both kernels must agree about where the target is. */
        RS_CHECK(pl == ps);

        /* And the wider kernel must not throw energy away. Total power, not
         * peak amplitude: on this fixture the peak legitimately falls, for the
         * reason in the comment above. */
        RS_CHECK(es >= el * 0.99);

        rs_slc_free(&lin);
        rs_slc_free(&sinc);
        rs_cphd_free(&c2);
    }

    rs_slc_free(&sub);
    rs_slc_free(&img);
    rs_cphd_free(&cphd);
    RS_TEST_END();
}
