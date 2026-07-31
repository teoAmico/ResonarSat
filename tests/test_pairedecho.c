/* The simulator's vibrating target, checked against theory rather than against
 * the code that consumes it.
 *
 * WHY THIS EXISTS. `sim_cphd` and `tests/rs_sim.h` are the only place in this
 * project where the answer is known in advance, and two things now rest
 * entirely on them: `--null-static`, which after 2026-07-30 is the only null
 * test trusted for a phase observable, and the observation-ratio threshold in
 * `rs_validate()`, which was measured on this fixture and nowhere else. Nothing
 * checked the fixture itself. A fixture wrong in a self-consistent way is
 * exactly the failure a suite built on it cannot see.
 *
 * The paired-echo model closes that, because it predicts something about the
 * FOCUSED IMAGE from physics and not from any code in this repository. A target
 * whose line-of-sight displacement is A*sin(w*t) modulates the phase history by
 * exp(j*B*sin(w*t)) with B = 4*pi*A/lambda, and expanding that in Bessel
 * functions gives a train of images:
 *
 *     order k sits at   k * f_v * lambda*R / (2*v_p)
 *     with amplitude    J_k(B) relative to J_0(B)
 *
 * Offsets depend only on geometry and the injected frequency; amplitudes only
 * on the injected amplitude, through a function with no free parameters. Both
 * are sharply structured -- J_0 has a zero, and the first ghost can outshine the
 * target -- so passing by accident is implausible.
 *
 * See docs/SIMULATOR-PAIRED-ECHO-ORACLE.md for the derivation. */

#include "resonarsat/focus.h"
#include "rs_sim.h"
#include "rs_test.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* The fixture's geometry, from rs_sim_scene(). Repeated rather than shared
 * because the point of this file is to check that function against numbers
 * derived independently of it. */
#define PE_FC      9.6e9
#define PE_VPLAT   7500.0
#define PE_HEIGHT  500000.0
#define PE_OFFSET  350000.0

/* Bessel function of the first kind, by its power series.
 *
 * Adequate and accurate here: the arguments are under 3 and the order under 4,
 * where the series converges in a handful of terms. Written out rather than
 * taken from libm's jn() so the test does not depend on a platform's choice of
 * implementation for a value it is asserting against. */
static double rs_pe_besselj(int n, double x)
{
    double term = 1.0;
    for (int i = 1; i <= n; i++) term *= (x / 2.0) / (double)i;
    double sum = term;
    for (int m = 1; m < 40; m++) {
        term *= -(x / 2.0) * (x / 2.0) / ((double)m * (double)(m + n));
        sum += term;
        if (fabs(term) < 1e-18 * fabs(sum)) break;
    }
    return sum;
}

/* Peak magnitude over range for each along-track cell.
 *
 * Row is along-track and column across it, per rs_focus_backproject(). Taking
 * the maximum over range rather than one row makes the profile insensitive to
 * exactly which range cell the target lands in. */
static void rs_pe_profile(const rs_slc_t *img, double *prof)
{
    for (size_t r = 0; r < img->n_az; r++) {
        double best = 0.0;
        for (size_t c = 0; c < img->n_rg; c++) {
            const double m = cabs(img->data[r * img->n_rg + c]);
            if (m > best) best = m;
        }
        prof[r] = best;
    }
}

/* The largest value within 'search' metres of an expected azimuth offset, and
 * where it actually fell. */
static double rs_pe_peak_near(const double *prof, size_t n, double cell,
                              double x_want, double search, double *x_found)
{
    const double centre = 0.5 * (double)(n - 1);
    const long lo = (long)floor(centre + (x_want - search) / cell);
    const long hi = (long)ceil(centre + (x_want + search) / cell);
    double best = 0.0;
    long at = -1;
    for (long i = lo; i <= hi; i++) {
        if (i < 0 || i >= (long)n) continue;
        if (prof[i] > best) { best = prof[i]; at = i; }
    }
    if (x_found) *x_found = (at >= 0) ? ((double)at - centre) * cell : NAN;
    return best;
}

/* Coherent loss of a ghost displaced 'x' metres along track.
 *
 * THE BESSEL AMPLITUDES ARE NOT THE WHOLE PREDICTION, and the difference is a
 * property of backprojection rather than of the fixture. A sideband at Doppler
 * offset k*f_v places its image where the linear phase -(4*pi*x/lambda)*u/R
 * cancels it, which is the offset above. But the true range difference between
 * the target and a grid point displaced by x is -(x*u)/R(u), not -(x*u)/R:
 * the range to the platform grows along the aperture. The residual,
 *
 *     phi(u) = (4*pi*x/lambda) * u * (1/R - 1/R(u)),
 *
 * reaches 0.94, 1.88 and 2.81 radians at the aperture edge for orders 1, 2 and
 * 3 here, and partially decoheres each ghost in proportion. Ignoring it makes
 * the higher orders look 15-45% too faint against theory, which is what a first
 * run of this test showed.
 *
 * A residual of about 12% remains after this correction, from terms beyond
 * first order in x and from the cross-track geometry. It is not modelled, and
 * the tolerances below allow for it rather than pretending otherwise. */
static double rs_pe_ghost_loss(double x, double lambda, double R, double v,
                               double dwell)
{
    const size_t n = 4001;
    double re = 0.0, im = 0.0;
    for (size_t i = 0; i < n; i++) {
        const double u = (-0.5 + (double)i / (double)(n - 1)) * v * dwell;
        const double Ru = sqrt(R * R + u * u);
        const double phi = (4.0 * M_PI * x / lambda) * u * (1.0 / R - 1.0 / Ru);
        re += cos(phi); im += sin(phi);
    }
    return sqrt(re * re + im * im) / (double)n;
}

/* Focus one vibrating target over the full aperture and return its azimuth
 * profile. The caller owns 'prof', which must hold 'n_az' doubles. */
static resonarsat_status_t rs_pe_focus(double vib_hz, double vib_amp_m,
                                       size_t n_az, size_t n_rg, double cell,
                                       double *prof)
{
    const rs_sim_tgt_t tg = { .x = 0.0, .y = 0.0, .z = 0.0, .rcs = 1.0,
                              .vib_freq = vib_hz, .vib_amp = vib_amp_m };
    rs_cphd_t c;
    resonarsat_status_t st = rs_sim_scene(&c, &tg, 1, 20.0, 400.0, 256, 0.5);
    if (st != RS_OK) return st;

    rs_grid_t g = { .origin = {0,0,0}, .n_x = n_az, .n_y = n_rg,
                    .dx = cell, .dy = cell, .height = 0.0 };
    rs_slc_t img;
    if ((st = rs_slc_alloc(&img, g.n_x, g.n_y)) != RS_OK) {
        rs_cphd_free(&c); return st;
    }
    st = rs_focus_backproject(&c, &g, 0, c.n_pulse, &img);
    if (st == RS_OK) rs_pe_profile(&img, prof);
    rs_slc_free(&img);
    rs_cphd_free(&c);
    return st;
}

int main(void)
{
    const double lambda = RS_C_LIGHT / PE_FC;
    const double R = sqrt(PE_HEIGHT * PE_HEIGHT + PE_OFFSET * PE_OFFSET);
    const double proj = PE_HEIGHT / R;              /* vertical onto the LOS */
    const double per_hz = lambda * R / (2.0 * PE_VPLAT);  /* ghost offset per Hz */
    const double az_res = lambda * R / (2.0 * PE_VPLAT * 20.0);

    const size_t n_az = 512, n_rg = 16;
    const double cell = 0.05;
    double *prof = malloc(n_az * sizeof *prof);
    RS_CHECK(prof != NULL);

    printf("  fixture: lambda %.6f m, R %.1f m, projection %.4f\n", lambda, R, proj);
    printf("           ghost offset %.4f m per Hz, azimuth resolution %.4f m\n",
           per_hz, az_res);

    /* ------------------------------------------------------------------
     * The default target: 2 Hz, 5 mm. B is 1.65, where the FIRST GHOST IS
     * BRIGHTER THAN THE TARGET -- an assertion no scale error can fake.
     * ------------------------------------------------------------------ */
    RS_CASE("ghost positions and amplitudes match the Bessel series");
    {
        const double f_v = 2.0, amp = 0.005;
        const double B = (4.0 * M_PI / lambda) * amp * proj;
        const double spacing = f_v * per_hz;
        printf("    B = %.4f, ghost spacing %.4f m (%.0f resolution cells)\n",
               B, spacing, spacing / az_res);

        RS_CHECK_OK(rs_pe_focus(f_v, amp, n_az, n_rg, cell, prof));

        const double j0 = rs_pe_besselj(0, B);
        double x_at = 0.0;
        const double p0 = rs_pe_peak_near(prof, n_az, cell, 0.0, 0.3, &x_at);
        RS_CHECK(p0 > 0.0);

        /* Orders 1 to 3, both sides. Order 3 is included because the sinc
         * sidelobe envelope of the main lobe has fallen far below it by 7.6 m
         * -- the concern recorded in the design note was about level rather
         * than position and does not survive the arithmetic. */
        for (int k = 1; k <= 3; k++) {
            const double want_ratio = fabs(rs_pe_besselj(k, B) / j0);
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const double x_want = (double)(sgn * k) * spacing;
                double x_got = 0.0;
                const double pk = rs_pe_peak_near(prof, n_az, cell, x_want,
                                                  0.35, &x_got);
                const double ratio = pk / p0;
                printf("    k=%+d at %+7.3f m (want %+7.3f): ratio %.4f  "
                       "Bessel %.4f x defocus %.3f = %.4f\n",
                       sgn * k, x_got, x_want, ratio, want_ratio,
                       rs_pe_ghost_loss(fabs(x_want), lambda, R, PE_VPLAT, 20.0),
                       want_ratio * rs_pe_ghost_loss(fabs(x_want), lambda, R,
                                                     PE_VPLAT, 20.0));

                /* Position, to half a resolution cell. This is the exact
                 * part of the prediction and is asserted tightly. */
                RS_CHECK_NEAR(x_got, x_want, 0.5 * az_res + cell);

                /* Amplitude, against the Bessel ratio corrected for the
                 * displaced-ghost defocus above. 25% covers the unmodelled
                 * residual; without the correction the higher orders miss by
                 * far more than any honest tolerance would admit. */
                const double loss = rs_pe_ghost_loss(fabs(x_want), lambda, R,
                                                     PE_VPLAT, 20.0);
                RS_CHECK_REL(ratio, want_ratio * loss, 0.25);
            }
        }

        /* The headline: at this amplitude the first ghost outshines the target
         * itself, which is the property that makes this a real check. */
        RS_CHECK(fabs(rs_pe_besselj(1, B) / j0) > 1.0);
        double x_g = 0.0;
        RS_CHECK(rs_pe_peak_near(prof, n_az, cell, spacing, 0.35, &x_g) > p0);
    }

    /* ------------------------------------------------------------------
     * Spacing is linear in frequency. Catches a wrong R or v_platform in a
     * way one frequency cannot.
     * ------------------------------------------------------------------ */
    RS_CASE("ghost spacing doubles when the frequency doubles");
    {
        RS_CHECK_OK(rs_pe_focus(4.0, 0.005, n_az, n_rg, cell, prof));
        double x_got = 0.0;
        const double want = 4.0 * per_hz;
        rs_pe_peak_near(prof, n_az, cell, want, 0.35, &x_got);
        printf("    4 Hz: first ghost at %+.3f m (want %+.3f, and 2 Hz gave "
               "%+.3f)\n", x_got, want, 2.0 * per_hz);
        RS_CHECK_NEAR(x_got, want, 0.5 * az_res + cell);
    }

    /* ------------------------------------------------------------------
     * The strongest assertion available: at the first zero of J_0 the
     * target's own peak vanishes and only the ghosts remain. A fixture that
     * scales the injected amplitude wrongly nulls at a different amplitude;
     * one that ignores amplitude never nulls at all.
     * ------------------------------------------------------------------ */
    RS_CASE("at the first zero of J_0 the target's own peak disappears");
    {
        const double B_null = 2.404826;   /* first zero of J_0 */
        const double amp = B_null * lambda / (4.0 * M_PI * proj);
        const double f_v = 2.0;
        const double spacing = f_v * per_hz;
        printf("    amplitude for J_0 = 0 is %.4f mm\n", 1000.0 * amp);

        RS_CHECK_OK(rs_pe_focus(f_v, amp, n_az, n_rg, cell, prof));
        double xc = 0.0, xg = 0.0;
        const double centre = rs_pe_peak_near(prof, n_az, cell, 0.0, 0.3, &xc);
        const double ghost  = rs_pe_peak_near(prof, n_az, cell, spacing, 0.35, &xg);
        printf("    centre %.4g, first ghost %.4g, ratio %.4f\n",
               centre, ghost, centre / ghost);

        /* J_1(2.405) is 0.519 and J_0 is zero, so the centre must collapse
         * relative to the ghost rather than merely shrink. */
        RS_CHECK(ghost > 0.0);
        RS_CHECK(centre < 0.2 * ghost);
    }

    /* ------------------------------------------------------------------
     * Control: no vibration, no ghosts. Guards against the checks above
     * passing on focusing artefacts that happen to sit near the offsets.
     * ------------------------------------------------------------------ */
    RS_CASE("a stationary target produces no ghost train");
    {
        RS_CHECK_OK(rs_pe_focus(0.0, 0.0, n_az, n_rg, cell, prof));
        double x0 = 0.0, x1 = 0.0;
        const double p0 = rs_pe_peak_near(prof, n_az, cell, 0.0, 0.3, &x0);
        const double p1 = rs_pe_peak_near(prof, n_az, cell, 2.0 * per_hz,
                                          0.35, &x1);
        printf("    centre %.4g, where a 2 Hz ghost would be %.4g (%.1f dB down)\n",
               p0, p1, 20.0 * log10(p1 / p0));
        RS_CHECK(p1 < 0.05 * p0);
    }

    free(prof);
    RS_TEST_END();
}
