/* Synthetic depth ground truth for Model A.
 *
 * WHAT THIS PROVES, AND WHAT IT DOES NOT. It validates the *inversion*: that the
 * steering matrix, the Kz mapping and the transform that inverts them are
 * mutually consistent, so that scatterers placed at known depths come back at
 * those depths. That is worth knowing and nothing else in the suite checks it.
 *
 * It does NOT validate the premise that along-orbit sub-aperture baselines carry
 * elevation information. It cannot: the fixture below synthesises its input FROM
 * that premise, by evaluating exactly the forward model the inversion assumes. A
 * reader who mistakes this for physical validation will over-trust every depth
 * this software reports.
 *
 * The distinction has teeth. tests/test_validation.c finds that Model A does not
 * distinguish a vibrating target from featureless ground on simulated radar
 * data. Two explanations were possible: the inversion is broken, or the geometry
 * carries no depth information. This test rules out the first. */

#include "resonarsat/tomo.h"
#include "rs_test.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Build a micro-motion result whose data vector Y is exactly what scatterers at
 * the given depths would produce under Eq. 22's forward model.
 *
 * For sub-aperture i the along-orbit baseline is B_i = A*(i/k), giving a
 * tomographic wavenumber Kz_i, and a scatterer of reflectivity 'amp' at depth
 * 'z' contributes amp*exp(j*Kz_i*z). The contributions sum, and the complex
 * total is split into the two coregistrator shift components the tracking stage
 * would have reported -- real part into azimuth, imaginary into range -- so that
 * rs_tomo_focus() reassembles precisely the Y that was synthesised.
 *
 * Pixel spacings are set to 1 m so shifts in pixels and metres coincide, which
 * keeps the fixture's arithmetic inspectable. */
static resonarsat_status_t make_depth_scene(rs_microm_t *m, size_t n_win, size_t n_looks,
                                            const rs_tomo_params_t *p,
                                            const double *depths, const double *amps,
                                            size_t n_scat)
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
    m->dt = 0.1;
    m->f_max = 5.0;
    m->az_spacing_m = 1.0;
    m->rg_spacing_m = 1.0;

    const double lambda_ac = rs_acoustic_wavelength(p->velocity, p->frequency,
                                                    p->convention);

    for (size_t w = 0; w < n_win; w++) {
        m->quality[w] = 0.9;
        for (size_t i = 0; i < n_looks; i++) {
            const double b_perp = p->aperture * ((double)i / (double)n_looks);
            const double kz = rs_tomo_wavenumber(b_perp, lambda_ac,
                                                 p->slant_range, p->incidence);

            double complex acc = 0.0;
            for (size_t sct = 0; sct < n_scat; sct++) {
                acc += amps[sct] * cexp(I * kz * depths[sct]);
            }

            const size_t idx = w * n_looks + i;
            m->disp_az[idx] = creal(acc);
            m->disp_rg[idx] = cimag(acc);
        }
    }
    return RS_OK;
}

/* Mean profile across windows, and the depth of its strongest cell above the
 * zero bin, refined parabolically. */
static double peak_depth(const rs_tomo_t *t)
{
    double *mean = calloc(t->n_depth, sizeof *mean);
    if (!mean) return -1.0;
    for (size_t j = 0; j < t->n_depth; j++) {
        for (size_t w = 0; w < t->n_win; w++) mean[j] += t->profile[w * t->n_depth + j];
        mean[j] /= (double)t->n_win;
    }
    size_t best = 1;
    for (size_t j = 2; j < t->n_depth; j++) if (mean[j] > mean[best]) best = j;

    double off = 0.0;
    if (best > 0 && best + 1 < t->n_depth) {
        const double y0 = mean[best - 1], y1 = mean[best], y2 = mean[best + 1];
        const double den = y0 - 2.0 * y1 + y2;
        if (den != 0.0) off = 0.5 * (y0 - y2) / den;
        if (off < -1.0) off = -1.0;
        if (off > 1.0) off = 1.0;
    }
    const double d = t->depth[best] + off * t->params.depth_cell;
    free(mean);
    return d;
}

int main(void)
{
    rs_tomo_params_t p;
    rs_tomo_params_default(&p);
    p.velocity    = 3000.0;
    p.frequency   = 12500.0;
    p.slant_range = 500000.0;
    p.aperture    = 75000.0;
    p.incidence   = 35.0 * M_PI / 180.0;
    p.depth_cell  = 0.5;
    p.depth_max   = 20.0;
    p.y_source    = RS_TOMO_Y_SHIFTS;

    const size_t n_looks = 128;
    const double dT = rs_tomo_resolution(
        rs_acoustic_wavelength(p.velocity, p.frequency, p.convention),
        p.slant_range, p.aperture);
    const double z_max = rs_tomo_max_depth(&p, n_looks);

    printf("  geometry: dT = %.3f m, unambiguous to %.2f m over %zu looks\n",
           dT, z_max, n_looks);
    RS_CHECK(p.depth_max <= z_max);

    /* ------------------------------------------------------------------
     * A single scatterer must come back where it was put.
     * ------------------------------------------------------------------ */
    RS_CASE("a single scatterer is recovered at its injected depth");
    {
        const double targets[] = { 3.0, 6.0, 9.0, 12.0 };
        for (size_t i = 0; i < sizeof targets / sizeof targets[0]; i++) {
            const double z = targets[i];
            const double amp = 1.0;

            rs_microm_t m;
            RS_CHECK_OK(make_depth_scene(&m, 4, n_looks, &p, &z, &amp, 1));

            rs_tomo_t t;
            RS_CHECK_OK(rs_tomo_focus(&m, NULL, &p, NULL, &t));

            const double got = peak_depth(&t);
            printf("    injected %.1f m -> recovered %.3f m (tolerance %.3f m)\n",
                   z, got, dT);
            RS_CHECK_NEAR(got, z, dT);

            rs_tomo_free(&t);
            rs_microm_free(&m);
        }
    }

    /* ------------------------------------------------------------------
     * The plan's own fixture: scatterers at 0, 3 and 6 m beneath a roof.
     * Separations here are several resolution cells, so all three should be
     * distinguishable rather than merging into one blur.
     * ------------------------------------------------------------------ */
    RS_CASE("three layered scatterers are separated");
    {
        const double depths[] = { 0.0, 3.0, 6.0 };
        const double amps[]   = { 1.0, 0.8, 0.6 };

        rs_microm_t m;
        RS_CHECK_OK(make_depth_scene(&m, 4, n_looks, &p, depths, amps, 3));

        rs_tomo_t t;
        RS_CHECK_OK(rs_tomo_focus(&m, NULL, &p, NULL, &t));

        /* Average profile, then check a local maximum sits near each injected
         * depth. Separation of 3 m against a resolution of dT means the layers
         * are resolved, not merged. */
        double *mean = calloc(t.n_depth, sizeof *mean);
        for (size_t j = 0; j < t.n_depth; j++) {
            for (size_t w = 0; w < t.n_win; w++) mean[j] += t.profile[w * t.n_depth + j];
            mean[j] /= (double)t.n_win;
        }

        for (size_t s = 1; s < 3; s++) {          /* skip z = 0, the zero bin */
            const size_t centre = (size_t)(depths[s] / p.depth_cell);
            const size_t span = (size_t)ceil(dT / p.depth_cell) + 1;
            double local = 0.0;
            for (size_t j = (centre > span ? centre - span : 0);
                 j < centre + span && j < t.n_depth; j++) {
                if (mean[j] > local) local = mean[j];
            }
            /* The layer must stand above the profile's own median level. */
            double sum = 0.0;
            for (size_t j = 1; j < t.n_depth; j++) sum += mean[j];
            const double avg = sum / (double)(t.n_depth - 1);
            printf("    layer at %.1f m: local peak %.4g vs profile mean %.4g\n",
                   depths[s], local, avg);
            RS_CHECK(local > avg);
        }

        free(mean);
        rs_tomo_free(&t);
        rs_microm_free(&m);
    }

    /* ------------------------------------------------------------------
     * The least-squares solver must agree with the DFT shortcut. The source
     * describes the steering matrix as a Fourier operator; if the two paths
     * disagree on a uniform depth grid, that description does not hold for
     * this geometry and the shortcut is invalid.
     * ------------------------------------------------------------------ */
    RS_CASE("the DFT and least-squares solvers agree");
    {
        const double z = 6.0, amp = 1.0;
        rs_microm_t m;
        RS_CHECK_OK(make_depth_scene(&m, 2, 32, &p, &z, &amp, 1));

        rs_tomo_params_t q = p;
        /* 32 looks reach about 7.3 m, a quarter of what 128 reach, so the grid
         * must shrink with them -- and must still contain the scatterer. */
        q.depth_max = 7.0;
        RS_CHECK(q.depth_max <= rs_tomo_max_depth(&q, 32));
        RS_CHECK(z < q.depth_max);

        rs_tomo_t td, tl;
        q.solver = RS_TOMO_SOLVER_DFT;
        RS_CHECK_OK(rs_tomo_focus(&m, NULL, &q, NULL, &td));
        q.solver = RS_TOMO_SOLVER_LSTSQ;
        RS_CHECK_OK(rs_tomo_focus(&m, NULL, &q, NULL, &tl));

        printf("    DFT peak %.3f m, least-squares peak %.3f m\n",
               peak_depth(&td), peak_depth(&tl));
        RS_CHECK_NEAR(peak_depth(&tl), peak_depth(&td), 2.0 * q.depth_cell);

        rs_tomo_free(&td);
        rs_tomo_free(&tl);
        rs_microm_free(&m);
    }

    /* ------------------------------------------------------------------
     * The recovered depth must be independent of how Y is scaled: doubling
     * every shift doubles the amplitude, not the depth. A depth that moved
     * with signal strength would mean the mapping had absorbed something it
     * should not have.
     * ------------------------------------------------------------------ */
    RS_CASE("recovered depth is independent of signal amplitude");
    {
        const double z = 6.0;
        const double a1 = 1.0, a2 = 25.0;

        rs_microm_t m1, m2;
        RS_CHECK_OK(make_depth_scene(&m1, 2, n_looks, &p, &z, &a1, 1));
        RS_CHECK_OK(make_depth_scene(&m2, 2, n_looks, &p, &z, &a2, 1));

        rs_tomo_t t1, t2;
        RS_CHECK_OK(rs_tomo_focus(&m1, NULL, &p, NULL, &t1));
        RS_CHECK_OK(rs_tomo_focus(&m2, NULL, &p, NULL, &t2));

        RS_CHECK_NEAR(peak_depth(&t2), peak_depth(&t1), 1e-6);

        rs_tomo_free(&t1); rs_tomo_free(&t2);
        rs_microm_free(&m1); rs_microm_free(&m2);
    }

    RS_TEST_END();
}
