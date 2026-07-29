/* Tomographic focusing: parameter contract, model behaviour, and the
 * sensitivity that is the whole point of the exercise. */

#include "resonarsat/tomo.h"
#include "resonarsat/geom.h"
#include "rs_test.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Build a synthetic micro-motion result whose displacement series carry a
 * single known oscillation, so that a focusing model has something with
 * definite spectral content to work on. */
static resonarsat_status_t make_microm(rs_microm_t *m, size_t n_win, size_t n_looks,
                                       double dt, double freq)
{
    memset(m, 0, sizeof *m);
    m->disp_az  = calloc(n_win * n_looks, sizeof *m->disp_az);
    m->disp_rg  = calloc(n_win * n_looks, sizeof *m->disp_rg);
    m->vel_los  = calloc(n_win * n_looks, sizeof *m->vel_los);
    m->disp_los = calloc(n_win * n_looks, sizeof *m->disp_los);
    m->phase    = calloc(n_win * n_looks, sizeof *m->phase);
    m->quality  = calloc(n_win, sizeof *m->quality);
    if (!m->disp_az || !m->disp_rg || !m->vel_los || !m->disp_los ||
        !m->phase || !m->quality) {
        rs_microm_free(m);
        return RS_ERR_ALLOC;
    }
    m->n_win = n_win;
    m->n_win_az = n_win;
    m->n_win_rg = 1;
    m->n_looks = n_looks;
    m->dt = dt;
    m->f_max = 1.0 / (2.0 * dt);

    m->az_spacing_m = 1.0;
    m->rg_spacing_m = 1.0;

    for (size_t w = 0; w < n_win; w++) {
        m->quality[w] = 0.9;
        for (size_t k = 0; k < n_looks; k++) {
            const double t = (double)k * dt;

            /* The two coregistrator shift components in quadrature, which is
             * the form Eq. 20 of the source gives for the motion:
             * r(t) = [a*cos(w0*t), b*sin(w0*t)]. Assembled into the complex Y
             * of Eq. 21 they form a single rotating phasor. */
            m->disp_az[w * n_looks + k] = 0.005 * cos(2.0 * M_PI * freq * t);
            m->disp_rg[w * n_looks + k] = 0.005 * sin(2.0 * M_PI * freq * t);

            m->disp_los[w * n_looks + k] = 0.005 * sin(2.0 * M_PI * freq * t);
            /* Velocity is the derivative of that displacement, same frequency. */
            m->vel_los[w * n_looks + k] =
                0.005 * 2.0 * M_PI * freq * cos(2.0 * M_PI * freq * t);
        }
    }
    return RS_OK;
}

int main(void)
{
    /* The parameter contract is the most important thing this file tests. The
     * two assumed constants have no defaults precisely so that a depth axis
     * cannot be fabricated from library values. */
    RS_CASE("velocity and frequency are required, with no defaults");
    {
        rs_tomo_params_t p;
        rs_tomo_params_default(&p);
        RS_CHECK_ERR(rs_tomo_params_check(&p), RS_ERR_ARG);   /* both unset */

        p.velocity = 3000.0;
        RS_CHECK_ERR(rs_tomo_params_check(&p), RS_ERR_ARG);   /* frequency still unset */

        p.frequency = 12500.0;
        p.depth_cell = 1.0;
        p.depth_max = 60.0;
        RS_CHECK_OK(rs_tomo_params_check(&p));
    }

    /* A grid finer than the resolution the geometry supports must be refused.
     * Interpolating onto more cells makes a result look sharper without adding
     * information, which is exactly the impression this project must not
     * create. */
    RS_CASE("a depth grid finer than the resolution is refused");
    {
        rs_tomo_params_t p;
        rs_tomo_params_default(&p);
        p.velocity = 6000.0;
        p.frequency = 12500.0;
        p.slant_range = 650000.0;
        p.aperture = 84000.0;
        /* dT here is about 0.93 m. */
        p.depth_cell = 0.05;
        RS_CHECK_ERR(rs_tomo_params_check(&p), RS_ERR_RANGE);
        p.depth_cell = 1.0;
        RS_CHECK_OK(rs_tomo_params_check(&p));
    }

    RS_CASE("a depth extent beyond the unambiguous range is refused");
    {
        rs_tomo_params_t q;
        rs_tomo_params_default(&q);
        q.velocity = 3000.0;
        q.frequency = 12500.0;
        q.slant_range = 500000.0;
        q.aperture = 75000.0;
        q.incidence = 35.0 * M_PI / 180.0;

        /* The limit scales linearly with the look count, so more sub-apertures
         * buy depth range -- at the cost of azimuth resolution per look. */
        const double z16  = rs_tomo_max_depth(&q, 16);
        const double z128 = rs_tomo_max_depth(&q, 128);
        RS_CHECK(z16 > 0.0);
        RS_CHECK_NEAR(z128 / z16, 8.0, 1e-9);
    }

    RS_CASE("literal Eq. 22 scaling is explicit and limited to Model A LSTSQ");
    {
        rs_tomo_params_t q;
        rs_tomo_params_default(&q);
        q.velocity = 3000.0;
        q.frequency = 12500.0;
        q.slant_range = 500000.0;
        q.aperture = 75000.0;
        q.incidence = 35.0 * M_PI / 180.0;
        q.eq22_literal_t = 1.0;

        /* The DFT implementation assumes the conventional uniformly sampled
         * spatial-wavenumber grid and must not silently accept this scale. */
        RS_CHECK_ERR(rs_tomo_params_check(&q), RS_ERR_ARG);
        q.solver = RS_TOMO_SOLVER_LSTSQ;
        RS_CHECK_OK(rs_tomo_params_check(&q));

        const double literal = rs_tomo_max_depth(&q, 32);
        q.eq22_literal_t = 0.0;
        const double conventional = rs_tomo_max_depth(&q, 32);
        RS_CHECK_NEAR(literal * 2.0 * M_PI, conventional, 1e-9);
    }

    RS_CASE("implausible wave speeds are refused");
    {
        rs_tomo_params_t p;
        rs_tomo_params_default(&p);
        p.frequency = 12500.0;
        p.velocity = 5.0;
        RS_CHECK_ERR(rs_tomo_params_check(&p), RS_ERR_RANGE);
        p.velocity = 1e6;
        RS_CHECK_ERR(rs_tomo_params_check(&p), RS_ERR_RANGE);
    }

    /* Model A, both solvers. */
    rs_microm_t m;
    RS_CASE("build synthetic vibration observations");
    RS_CHECK_OK(make_microm(&m, 8, 32, 0.05, 2.0));

    rs_spectrum_t spec;
    RS_CHECK_OK(rs_spectrum_compute(&m, RS_SPEC_VELOCITY, &spec));

    RS_CASE("the spectrum recovers the injected frequency");
    {
        /* The series was built at 2 Hz; the estimator must find it within one
         * spectral bin. This is the check that the whole micro-motion chain
         * exists to satisfy. */
        RS_CHECK_NEAR(spec.dominant_freq[0], 2.0, spec.df * 1.5);
    }

    rs_tomo_params_t p;
    rs_tomo_params_default(&p);
    p.velocity = 3000.0;
    p.frequency = 12500.0;
    /* 32 sub-apertures at this geometry represent about 7.3 m without folding
     * (see rs_tomo_max_depth), so the extent must stay inside that. */
    p.depth_max = 6.0;
    p.depth_cell = 1.0;
    p.slant_range = 500000.0;
    p.aperture = 75000.0;

    RS_CASE("model A (DFT solver) produces a profile");
    {
        rs_tomo_t t;
        p.model = RS_TOMO_MODEL_A;
        p.solver = RS_TOMO_SOLVER_DFT;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t));
        RS_CHECK(t.n_depth == 7);
        RS_CHECK(t.lambda_ac > 0.0);
        RS_CHECK(t.dT > 0.0);

        double e = 0.0;
        for (size_t i = 0; i < t.n_win * t.n_depth; i++) e += t.profile[i];
        RS_CHECK(e > 0.0);
        rs_tomo_free(&t);
    }

    RS_CASE("model A (least-squares solver) produces a profile");
    {
        rs_tomo_t t;
        p.model = RS_TOMO_MODEL_A;
        p.solver = RS_TOMO_SOLVER_LSTSQ;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t));
        double e = 0.0;
        for (size_t i = 0; i < t.n_win * t.n_depth; i++) e += t.profile[i];
        RS_CHECK(e > 0.0);
        rs_tomo_free(&t);
    }

    RS_CASE("model B runs from spectra");
    {
        rs_tomo_t t;
        p.model = RS_TOMO_MODEL_B;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t));
        rs_tomo_free(&t);
    }

    /* Model C needs its own geometry, not Model A's. It inverts with the radar
     * wavelength over a real perpendicular baseline spread -- here 0 to 620 m,
     * a plausible repeat-pass tube -- so both its resolution cell and its
     * unambiguous extent differ from Model A's by nearly a factor of four.
     * Reusing Model A's 75 km "aperture" and 6 m extent would ask Model C for a
     * profile four times deeper than its geometry represents. */
    double *b = malloc(m.n_looks * sizeof *b);
    RS_CHECK(b != NULL);
    for (size_t i = 0; i < m.n_looks; i++) b[i] = (double)i * 20.0;

    rs_tomo_params_t pc = p;
    pc.model = RS_TOMO_MODEL_C;
    pc.aperture = b[m.n_looks - 1];
    pc.radar_wavelength = 0.031;
    pc.depth_cell = ceil(rs_tomo_resolution(pc.radar_wavelength, pc.slant_range,
                                            pc.aperture));
    pc.depth_max = floor(fmin(4.0 * pc.depth_cell,
                              rs_tomo_max_depth(&pc, m.n_looks)));

    RS_CASE("model C refuses to run without genuine baselines");
    {
        rs_tomo_t t;
        RS_CHECK_ERR(rs_tomo_focus(&m, &spec, &pc, NULL, &t), RS_ERR_ARG);
    }

    RS_CASE("model C runs when baselines are supplied");
    {
        rs_tomo_t t;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &pc, b, &t));
        rs_tomo_free(&t);
    }
    free(b);

    /* THE test for this stage. The depth axis is scaled by assumed constants,
     * so halving the assumed frequency must double every depth. If this ever
     * stops holding, either the chain has been broken or a measured quantity
     * has crept into the depth mapping -- both of which would be serious. */
    RS_CASE("depth scale tracks the assumed constants exactly");
    {
        p.model = RS_TOMO_MODEL_A;
        p.solver = RS_TOMO_SOLVER_DFT;

        rs_tomo_t t1, t2;
        p.frequency = 12500.0;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t1));
        p.frequency = 6250.0;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &t2));

        /* Halving f doubles lambda_ac and therefore doubles dT. */
        RS_CHECK_NEAR(t2.lambda_ac / t1.lambda_ac, 2.0, 1e-9);
        RS_CHECK_NEAR(t2.dT / t1.dT, 2.0, 1e-9);

        rs_tomo_free(&t1);
        rs_tomo_free(&t2);
    }

    /* The sharpest statement of the headline result, and the one a reader can
     * check without trusting a regression fit.
     *
     * Scale the depth grid in proportion to the assumed wavelength and the two
     * tomograms are not merely similar, they are the same numbers. Cell k holds
     * an identical value in both while being labelled at twice the depth. The
     * fitted slope of one is this fact measured indirectly; this is the fact.
     *
     * It also says exactly what a tomogram image is. Rendering slice k under
     * two assumptions produces pixel-identical pictures with captions differing
     * by the ratio of the assumptions, which is checkable with cmp(1). */
    RS_CASE("the same slice under two assumptions holds identical values");
    {
        p.model = RS_TOMO_MODEL_A;
        p.solver = RS_TOMO_SOLVER_DFT;

        rs_tomo_params_t q1 = p, q2 = p;
        q1.frequency = 12500.0;
        q1.depth_cell = 1.0;
        q1.depth_max = 6.0;

        /* Half the frequency doubles lambda_ac, so the grid doubles with it. */
        q2.frequency = 6250.0;
        q2.depth_cell = 2.0;
        q2.depth_max = 12.0;

        rs_tomo_t t1, t2;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &q1, NULL, &t1));
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &q2, NULL, &t2));

        RS_CHECK(t1.n_depth == t2.n_depth);
        RS_CHECK(t1.n_win == t2.n_win);
        RS_CHECK_NEAR(t2.lambda_ac / t1.lambda_ac, 2.0, 1e-12);

        /* Every cell of every window, not a summary statistic. */
        double worst = 0.0;
        for (size_t w = 0; w < t1.n_win; w++) {
            for (size_t k = 0; k < t1.n_depth; k++) {
                const size_t i = w * t1.n_depth + k;
                const double d = fabs(t1.profile[i] - t2.profile[i]);
                if (d > worst) worst = d;
            }
        }
        printf("    worst difference over %zu windows x %zu cells: %.3g\n",
               t1.n_win, t1.n_depth, worst);
        RS_CHECK(worst < 1e-9);

        /* And the labels differ by exactly the ratio of the assumptions. */
        for (size_t k = 1; k < t1.n_depth; k++) {
            RS_CHECK_NEAR(t2.depth[k] / t1.depth[k], 2.0, 1e-12);
        }
        printf("    slice %zu is labelled %.3f m and %.3f m\n",
               t1.n_depth / 2, t1.depth[t1.n_depth / 2], t2.depth[t1.n_depth / 2]);

        rs_tomo_free(&t1);
        rs_tomo_free(&t2);
    }

    RS_CASE("the two wavelength conventions differ by exactly two");
    {
        rs_tomo_t tp, tq;
        p.frequency = 12500.0;
        p.convention = RS_WAVELEN_PAPER;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &tp));
        p.convention = RS_WAVELEN_PATENT;
        RS_CHECK_OK(rs_tomo_focus(&m, &spec, &p, NULL, &tq));

        RS_CHECK_NEAR(tq.lambda_ac / tp.lambda_ac, 2.0, 1e-9);
        rs_tomo_free(&tp);
        rs_tomo_free(&tq);
    }

    rs_spectrum_free(&spec);
    rs_microm_free(&m);
        /* --------------------------------------------------------------------
     * Eq. 24 states h = A_dagger Y, the Moore-Penrose pseudoinverse. This
     * implementation solves the normal equations, which equals the
     * pseudoinverse only when A has full COLUMN rank -- that is, when there are
     * at least as many sub-looks as depth cells.
     *
     * That condition turns out to be guaranteed rather than assumed. Combining
     * the two guards this project already enforces:
     *
     *     depth cell   >= lambda*R/(2A)                    (resolution)
     *     depth extent <= 2*pi*k/Kz_max                    (unambiguous range)
     *
     * and substituting Kz_max = 4*pi*A/(lambda*R*sin(theta)) gives
     *
     *     F = extent/cell <= k * sin(theta) <= k
     *
     * So the steering matrix is never underdetermined in a valid run, and the
     * normal-equation solve at zero regularisation IS the exact pseudoinverse.
     * -------------------------------------------------------------------- */
    RS_CASE("the depth grid can never outnumber the sub-looks");
    {
        rs_tomo_params_t tp;
        rs_tomo_params_default(&tp);
        tp.velocity = 3000.0; tp.frequency = 200.0;
        tp.slant_range = 650000.0; tp.aperture = 42000.0;

        for (double inc_deg = 20.0; inc_deg <= 80.0; inc_deg += 20.0) {
            tp.incidence = inc_deg * M_PI / 180.0;
            const double lam = rs_acoustic_wavelength(tp.velocity, tp.frequency,
                                                      tp.convention);
            const double cell = rs_tomo_resolution(lam, tp.slant_range, tp.aperture);
            for (size_t k = 8; k <= 128; k *= 2) {
                const double zmax = rs_tomo_max_depth(&tp, k);
                const double F = zmax / cell;
                printf("    %2.0f deg, %3zu looks: at most %5.1f cells "
                       "(%.3f k)\n", inc_deg, k, F, F / (double)k);
                RS_CHECK(F <= (double)k + 1e-9);
            }
        }
    }

    RS_CASE("the pseudoinverse at zero regularisation solves an overdetermined grid");
    {
        rs_microm_t pm;
        rs_spectrum_t psp;
        memset(&pm, 0, sizeof pm);
        memset(&psp, 0, sizeof psp);

        const size_t k = 16, n_win = 4;
        pm.n_looks = k; pm.n_win = n_win; pm.n_win_az = 2; pm.n_win_rg = 2;
        pm.win_az = pm.win_rg = 8; pm.stride_az = pm.stride_rg = 8;
        pm.dt = 0.1; pm.az_spacing_m = 1.0; pm.rg_spacing_m = 1.0;
        pm.disp_az = calloc(n_win * k, sizeof *pm.disp_az);
        pm.disp_rg = calloc(n_win * k, sizeof *pm.disp_rg);
        pm.disp_los = calloc(n_win * k, sizeof *pm.disp_los);
        pm.vel_los = calloc(n_win * k, sizeof *pm.vel_los);
        pm.quality = calloc(n_win, sizeof *pm.quality);
        RS_CHECK(pm.disp_az && pm.disp_rg && pm.disp_los && pm.vel_los && pm.quality);
        for (size_t w = 0; w < n_win; w++) {
            pm.quality[w] = 1.0;
            for (size_t i = 0; i < k; i++) {
                pm.disp_az[w * k + i] = 0.01 * cos(0.4 * (double)i);
                pm.disp_rg[w * k + i] = 0.01 * sin(0.4 * (double)i);
            }
        }

        rs_tomo_params_t tp;
        rs_tomo_params_default(&tp);
        tp.velocity = 3000.0; tp.frequency = 200.0;
        tp.slant_range = 650000.0; tp.aperture = 42000.0;
        tp.incidence = 0.6; tp.solver = RS_TOMO_SOLVER_LSTSQ;
        tp.regularisation = 0.0;      /* exactly Eq. 24, no ridge */
        tp.window = 0;                /* and no taper the patent does not state */
        tp.remove_y_mean = 0;         /* and the raw Y of Eq. 21 */

        const double cell = rs_tomo_resolution(
            rs_acoustic_wavelength(tp.velocity, tp.frequency, tp.convention),
            tp.slant_range, tp.aperture);
        tp.depth_cell = cell; tp.depth_max = cell * 8.0;

        rs_tomo_t a;
        RS_CHECK_OK(rs_tomo_focus(&pm, &psp, &tp, NULL, &a));
        RS_CHECK(a.n_depth <= k);

        int finite = 1;
        for (size_t i = 0; i < a.n_win * a.n_depth; i++)
            if (!isfinite(a.profile[i])) { finite = 0; break; }
        printf("    %zu looks, %zu cells, mu = 0: all finite = %s\n",
               k, a.n_depth, finite ? "yes" : "NO");
        RS_CHECK(finite);

        rs_tomo_free(&a);
        free(pm.disp_az); free(pm.disp_rg); free(pm.disp_los);
        free(pm.vel_los); free(pm.quality);
    }

    RS_CASE("metadata patent-chain claim requires the full selected chain");
    {
        rs_tomo_t t;
        memset(&t, 0, sizeof t);
        rs_tomo_params_default(&t.params);
        t.params.model = RS_TOMO_MODEL_A;
        t.params.solver = RS_TOMO_SOLVER_LSTSQ;
        t.params.y_source = RS_TOMO_Y_SHIFTS;
        t.params.convention = RS_WAVELEN_PATENT;
        t.params.velocity = 3000.0;
        t.params.frequency = 1000.0;
        t.params.regularisation = 0.0;
        t.params.window = 0;
        t.params.remove_y_mean = 0;
        t.params.patent_exact = 1;
        t.params.subap_window = 0;
        t.params.coherence_min = 0.0;
        t.params.pair_reference = 1;

        char *buf = NULL;
        size_t len = 0;
        FILE *f = open_memstream(&buf, &len);
        RS_CHECK(f != NULL);
        rs_tomo_write_metadata(&t, f);
        fclose(f);
        RS_CHECK(strstr(buf, "PATENT-CHAIN INTERPRETATION") != NULL);
        RS_CHECK(strstr(buf, "NOT LITERAL EQ. 22") != NULL);
        free(buf);

        t.params.eq22_literal_t = 1.0;
        buf = NULL;
        len = 0;
        f = open_memstream(&buf, &len);
        RS_CHECK(f != NULL);
        rs_tomo_write_metadata(&t, f);
        fclose(f);
        RS_CHECK(strstr(buf, "EXPERIMENTAL LITERAL EQ. 22") != NULL);
        RS_CHECK(strstr(buf, "eq22_kz_scale") != NULL);
        free(buf);
        t.params.eq22_literal_t = 0.0;

        t.params.subap_window = 1;
        buf = NULL;
        len = 0;
        f = open_memstream(&buf, &len);
        RS_CHECK(f != NULL);
        rs_tomo_write_metadata(&t, f);
        fclose(f);
        RS_CHECK(strstr(buf, "NOT THE UNCONDITIONED PATENT CHAIN") != NULL);
        free(buf);

        /* The front end counts too. Eq. 21's Y is the offset of each look's
         * slave from its own master; a run that tracked everything against one
         * fixed look computed a different observable, and no exactness in
         * Eqs. 22-24 downstream makes that the published method. */
        t.params.subap_window = 0;
        t.params.pair_reference = 0;
        buf = NULL;
        len = 0;
        f = open_memstream(&buf, &len);
        RS_CHECK(f != NULL);
        rs_tomo_write_metadata(&t, f);
        fclose(f);
        RS_CHECK(strstr(buf, "NOT THE UNCONDITIONED PATENT CHAIN") != NULL);
        free(buf);
    }

    RS_TEST_END();
}
