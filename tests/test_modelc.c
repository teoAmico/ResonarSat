/* Model A against Model C: the disputed claim, quantified.
 *
 * Model C is classical multi-baseline SAR tomography — genuine perpendicular
 * baselines from a repeat-pass stack, the electromagnetic wavelength, physics
 * nobody disputes. Model A is the along-orbit sub-aperture construction of
 * Biondi & Malanga (2022), which substitutes an acoustic wavelength and treats
 * the arc flown during a single dwell as the tomographic baseline.
 *
 * Both are implemented here from the same interfaces, so the comparison is not
 * an argument. The plan (§Phase 6 item 6) asks for it as a number: if Model A
 * reports finer depth structure than Model C can support on comparable data,
 * that gap is the claim under dispute.
 *
 * The test does two things. It checks Model C recovers known heights from a
 * synthetic multi-baseline stack — establishing that the uncontested model works
 * in this implementation and is a fair reference. Then it computes both models'
 * depth resolution from their own stated geometry and reports the ratio. */

#include "resonarsat/tomo.h"
#include "rs_test.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Build a micro-motion result whose data vector is what scatterers at the given
 * heights produce under Model C's forward model: genuine perpendicular
 * baselines and the radar wavelength.
 *
 * This is the Model C counterpart of the fixture in test_depth.c, and it shares
 * that fixture's limitation — it synthesises its input from the model's own
 * premise, so it validates the inversion rather than the physics. The difference
 * is that for Model C the premise is uncontested: perpendicular baselines really
 * do carry elevation information, which is how conventional SAR tomography
 * works. */
static resonarsat_status_t make_baseline_scene(rs_microm_t *m, size_t n_win,
                                               const double *baselines, size_t n_base,
                                               const rs_tomo_params_t *p,
                                               const double *heights,
                                               const double *amps, size_t n_scat)
{
    memset(m, 0, sizeof *m);
    m->disp_az  = calloc(n_win * n_base, sizeof *m->disp_az);
    m->disp_rg  = calloc(n_win * n_base, sizeof *m->disp_rg);
    m->vel_los  = calloc(n_win * n_base, sizeof *m->vel_los);
    m->disp_los = calloc(n_win * n_base, sizeof *m->disp_los);
    m->phase    = calloc(n_win * n_base, sizeof *m->phase);
    m->quality  = calloc(n_win, sizeof *m->quality);
    if (!m->disp_az || !m->disp_rg || !m->vel_los || !m->disp_los ||
        !m->phase || !m->quality) {
        rs_microm_free(m);
        return RS_ERR_ALLOC;
    }

    m->n_win = n_win;
    m->n_win_az = n_win;
    m->n_win_rg = 1;
    m->n_looks = n_base;
    m->dt = 0.1;
    m->f_max = 5.0;
    m->az_spacing_m = 1.0;
    m->rg_spacing_m = 1.0;

    const double lambda = p->radar_wavelength;

    for (size_t w = 0; w < n_win; w++) {
        m->quality[w] = 0.9;
        for (size_t i = 0; i < n_base; i++) {
            const double kz = rs_tomo_wavenumber(baselines[i], lambda,
                                                 p->slant_range, p->incidence);
            double complex acc = 0.0;
            for (size_t s = 0; s < n_scat; s++) {
                acc += amps[s] * cexp(I * kz * heights[s]);
            }
            const size_t idx = w * n_base + i;
            m->disp_az[idx] = creal(acc);
            m->disp_rg[idx] = cimag(acc);
        }
    }
    return RS_OK;
}

/* Depth of the strongest cell above the zero bin, parabolically refined. */
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
    /* A repeat-pass stack with a realistic perpendicular baseline spread.
     * Sentinel-1 holds its orbital tube to a few hundred metres; TerraSAR-X and
     * COSMO-SkyMed stacks reach a kilometre or two. 1 km is generous. */
    enum { N_BASE = 24 };
    const double b_max = 1000.0;
    double baselines[N_BASE];
    for (size_t i = 0; i < N_BASE; i++) {
        baselines[i] = b_max * (double)i / (double)(N_BASE - 1);
    }

    rs_tomo_params_t pc;
    rs_tomo_params_default(&pc);
    pc.model = RS_TOMO_MODEL_C;
    pc.velocity = 3000.0;        /* unused by Model C, but required to validate */
    pc.frequency = 12500.0;
    pc.radar_wavelength = 0.031;
    pc.slant_range = 500000.0;
    pc.incidence = 35.0 * M_PI / 180.0;
    pc.aperture = b_max;         /* for Model C the baseline spread IS the aperture */

    const double dT_c = rs_tomo_resolution(pc.radar_wavelength, pc.slant_range, b_max);
    const double z_max_c = rs_tomo_max_depth(&pc, N_BASE);

    printf("  Model C: %d baselines to %.0f m, radar wavelength %.3f m\n",
           N_BASE, b_max, pc.radar_wavelength);
    printf("           resolution %.2f m, unambiguous to %.1f m\n", dT_c, z_max_c);

    pc.depth_cell = ceil(dT_c);
    pc.depth_max = fmin(60.0, floor(z_max_c));

    /* ------------------------------------------------------------------
     * Model C recovers known heights. Establishes it as a fair reference
     * rather than a straw man.
     * ------------------------------------------------------------------ */
    RS_CASE("model C recovers heights from genuine perpendicular baselines");
    {
        const double targets[] = { 10.0, 20.0, 30.0 };
        for (size_t i = 0; i < sizeof targets / sizeof targets[0]; i++) {
            if (targets[i] > pc.depth_max) continue;
            const double amp = 1.0;

            rs_microm_t m;
            RS_CHECK_OK(make_baseline_scene(&m, 4, baselines, N_BASE, &pc,
                                            &targets[i], &amp, 1));

            rs_tomo_t t;
            RS_CHECK_OK(rs_tomo_focus(&m, NULL, &pc, baselines, &t));

            const double got = peak_depth(&t);
            printf("    injected %.1f m -> recovered %.2f m (resolution %.2f m)\n",
                   targets[i], got, dT_c);
            RS_CHECK_NEAR(got, targets[i], 1.5 * dT_c);

            rs_tomo_free(&t);
            rs_microm_free(&m);
        }
    }

    /* ------------------------------------------------------------------
     * THE COMPARISON. Both models' depth resolution, each computed from the
     * geometry it claims for itself.
     * ------------------------------------------------------------------ */
    RS_CASE("how much finer does Model A claim to resolve than physics supports?");
    {
        rs_tomo_params_t pa;
        rs_tomo_params_default(&pa);
        pa.model = RS_TOMO_MODEL_A;
        pa.velocity = 3000.0;
        pa.frequency = 12500.0;
        pa.slant_range = 500000.0;
        pa.incidence = 35.0 * M_PI / 180.0;
        pa.aperture = 75000.0;     /* the along-orbit arc flown in a 10 s dwell */

        const double lambda_ac = rs_acoustic_wavelength(pa.velocity, pa.frequency,
                                                        pa.convention);
        const double dT_a = rs_tomo_resolution(lambda_ac, pa.slant_range, pa.aperture);

        printf("\n    %-34s %12s %12s\n", "", "Model A", "Model C");
        printf("    %-34s %12s %12s\n", "wavelength used",
               "acoustic", "radar");
        printf("    %-34s %12.4f %12.4f\n", "  ... in metres", lambda_ac,
               pc.radar_wavelength);
        printf("    %-34s %12.0f %12.0f\n", "baseline extent (m)",
               pa.aperture, b_max);
        printf("    %-34s %12.3f %12.3f\n", "depth resolution (m)", dT_a, dT_c);
        printf("\n    Model A claims %.1fx finer depth resolution.\n", dT_c / dT_a);

        /* Where the gain comes from. The acoustic substitution makes resolution
         * WORSE here, since the acoustic wavelength is longer than the radar
         * one; every bit of the advantage and more comes from the baseline. */
        const double from_wavelength = pc.radar_wavelength / lambda_ac;
        const double from_baseline = pa.aperture / b_max;
        printf("      from the wavelength substitution: %.2fx\n", from_wavelength);
        printf("      from treating a %.0f km along-track arc\n"
               "      as an elevation baseline:         %.0fx\n",
               pa.aperture / 1000.0, from_baseline);

        RS_CHECK_NEAR(from_wavelength * from_baseline, dT_c / dT_a, 1e-6);

        /* The claim under dispute, as a number. */
        RS_CHECK(dT_c / dT_a > 10.0);

        printf("\n    The entire advantage is the baseline term. Perpendicular\n"
               "    baselines carry elevation information because the two\n"
               "    observations see the target from different elevation angles;\n"
               "    an along-track arc does not, however long it is. The\n"
               "    sensitivity sweep (test_validation.c) measures the\n"
               "    consequence directly: Model A's depth axis tracks the assumed\n"
               "    acoustic wavelength with slope 1.045, carrying no independent\n"
               "    depth information at any resolution.\n\n");
    }

    /* ------------------------------------------------------------------
     * Model C's resolution responds to its baseline spread, as real
     * tomography must. A model whose resolution ignored its baselines would
     * not be a fair reference.
     * ------------------------------------------------------------------ */
    RS_CASE("model C resolution scales with baseline extent");
    {
        const double r1 = rs_tomo_resolution(pc.radar_wavelength, pc.slant_range, 500.0);
        const double r2 = rs_tomo_resolution(pc.radar_wavelength, pc.slant_range, 2000.0);
        printf("    500 m baselines -> %.2f m;  2000 m -> %.2f m\n", r1, r2);
        RS_CHECK_NEAR(r1 / r2, 4.0, 1e-9);
    }

    RS_TEST_END();
}
