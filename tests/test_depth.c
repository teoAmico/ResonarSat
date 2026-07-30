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

/* One look's contribution, for a signal whose SHAPE departs from the steering
 * assumption.
 *
 * 'i' is the look index and 'n' the look count, so a shape can vary along the
 * aperture -- which is what damping, chirp and quadratic phase all do. */
typedef double complex (*rs_shape_fn)(size_t i, size_t n, double kz, double z);

/* Build a micro-motion result from an arbitrary per-look signal shape.
 *
 * make_depth_scene() injects exactly the tone the matched filter is built to
 * recover, which tests conditioning and signal-to-noise but not robustness. This
 * injects whatever 'shape' returns, so the inversion can be asked what it does
 * with a signal the model does not describe. Same conventions otherwise. */
static resonarsat_status_t make_shaped_scene(rs_microm_t *m, size_t n_win,
                                             size_t n_looks,
                                             const rs_tomo_params_t *p,
                                             double z, rs_shape_fn shape,
                                             double noise)
{
    /* xorshift rather than rand(), so a fixture is identical on every
     * platform and a failure is reproducible from the source alone. */
    unsigned rs_ = 20260730u;
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
            double complex v = shape(i, n_looks, kz, z);
            if (noise > 0.0) {
                /* Two uniform draws summed to a rough normal, scaled to the
                 * requested amplitude. Added to both components, because a
                 * tracker's noise is on the shift it reports, not on a
                 * hypothetical analytic signal behind it. */
                double d[2];
                for (int k = 0; k < 2; k++) {
                    double s = 0.0;
                    for (int j = 0; j < 4; j++) {
                        rs_ ^= rs_ << 13; rs_ ^= rs_ >> 17; rs_ ^= rs_ << 5;
                        s += (double)(rs_ % 2000) / 1000.0 - 1.0;
                    }
                    d[k] = noise * s * 0.5;
                }
                v += d[0] + I * d[1];
            }
            const size_t idx = w * n_looks + i;
            m->disp_az[idx] = creal(v);
            m->disp_rg[idx] = cimag(v);
        }
    }
    return RS_OK;
}

/* The matched tone the steering matrix is built from: the positive control. */
static double complex shape_matched(size_t i, size_t n, double kz, double z)
{ (void)i; (void)n; return cexp(I * kz * z); }

/* A resonance losing amplitude along the aperture, moderately. */
static double complex shape_damped(size_t i, size_t n, double kz, double z)
{ return cexp(I * kz * z) * exp(-1.5 * (double)i / (double)n); }

/* The same, losing most of its amplitude within the aperture. */
static double complex shape_damped_heavy(size_t i, size_t n, double kz, double z)
{ return cexp(I * kz * z) * exp(-4.0 * (double)i / (double)n); }

/* A dispersive return whose apparent depth drifts along the aperture. */
static double complex shape_chirped(size_t i, size_t n, double kz, double z)
{ return cexp(I * kz * z * (1.0 + 0.6 * (double)i / (double)n)); }

/* A quadratic phase error across the aperture, as a defocus would give. */
static double complex shape_quadratic(size_t i, size_t n, double kz, double z)
{
    const double u = (double)i / (double)n;
    return cexp(I * (kz * z + 4.0 * u * u));
}

/* A feature too shallow for the geometry to resolve. */
static double complex shape_subres(size_t i, size_t n, double kz, double z)
{ (void)i; (void)n; return cexp(I * kz * z * 0.15); }

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

    /* ------------------------------------------------------------------
     * TWO STABILITY PROBES, and only one of them works.
     *
     * The unambiguous depth is 2*pi*n_looks/kz_max -- linear in the look
     * count -- so a reflector deeper than it folds back and appears shallow.
     * A folded peak is indistinguishable from a real one in a single
     * tomogram: same shape, same contrast, plausible depth. The question is
     * what perturbation separates them.
     *
     * Both probes come from the independent reproduction at
     * github.com/Hassanforeman/subsurface-sar-tomo, which found that one
     * discriminates and one does not, and kept the failure on purpose. Both
     * are reproduced here for the same reason: the negative is the more
     * useful of the two, because it is the one somebody would otherwise
     * reach for.
     * ------------------------------------------------------------------ */

    const double zmax32 = rs_tomo_max_depth(&p, 32);

    RS_CASE("look-count stability separates a real reflector from a folded one");
    {
        const size_t counts[] = { 32, 48, 64, 96, 128 };
        const size_t n_counts = sizeof counts / sizeof counts[0];

        /* Comfortably inside the smallest configuration's range, and well
         * beyond it -- 1.8x, so it wraps at 32 and 48 looks and is in range
         * from 64 up. */
        const double z_real   = 0.35 * zmax32;
        const double z_folded = 1.80 * zmax32;

        double rec_real[5], rec_folded[5];
        for (size_t i = 0; i < n_counts; i++) {
            rs_tomo_params_t q = p;
            /* Each configuration searches its OWN full unambiguous range, which
             * is what a caller following the depth guard would do. That is what
             * makes the folded case move: at 32 and 48 looks the reflector is
             * beyond the range and wraps to a different shallow depth each
             * time, and from 64 looks up the range finally covers it and it
             * appears where it really is. */
            q.depth_max = rs_tomo_max_depth(&p, counts[i]);

            const double amp = 1.0;
            rs_microm_t mr, mf;
            RS_CHECK_OK(make_depth_scene(&mr, 4, counts[i], &q, &z_real, &amp, 1));
            RS_CHECK_OK(make_depth_scene(&mf, 4, counts[i], &q, &z_folded, &amp, 1));

            rs_tomo_t tr, tf;
            RS_CHECK_OK(rs_tomo_focus(&mr, NULL, &q, NULL, &tr));
            RS_CHECK_OK(rs_tomo_focus(&mf, NULL, &q, NULL, &tf));
            rec_real[i] = peak_depth(&tr);
            rec_folded[i] = peak_depth(&tf);
            rs_tomo_free(&tr); rs_tomo_free(&tf);
            rs_microm_free(&mr); rs_microm_free(&mf);
        }

        double mr_ = 0.0, mf_ = 0.0;
        for (size_t i = 0; i < n_counts; i++) { mr_ += rec_real[i]; mf_ += rec_folded[i]; }
        mr_ /= (double)n_counts; mf_ /= (double)n_counts;
        double sr = 0.0, sf = 0.0;
        for (size_t i = 0; i < n_counts; i++) {
            sr += (rec_real[i] - mr_) * (rec_real[i] - mr_);
            sf += (rec_folded[i] - mf_) * (rec_folded[i] - mf_);
        }
        sr = sqrt(sr / (double)n_counts);
        sf = sqrt(sf / (double)n_counts);

        printf("    unambiguous depth at 32 looks: %.2f m\n", zmax32);
        printf("    real   (z = %.2f m):", z_real);
        for (size_t i = 0; i < n_counts; i++) printf(" %6.2f", rec_real[i]);
        printf("   sd %.3f m\n", sr);
        printf("    folded (z = %.2f m):", z_folded);
        for (size_t i = 0; i < n_counts; i++) printf(" %6.2f", rec_folded[i]);
        printf("   sd %.3f m\n", sf);

        /* The guard: a real reflector holds its depth across look counts. */
        RS_CHECK(sr < dT);
        /* And a folded one does not, by a wide margin. */
        RS_CHECK(sf > 4.0 * sr);
    }

    /* ------------------------------------------------------------------
     * SIGNATURE MODEL MISMATCH.
     *
     * Every positive control elsewhere in this file injects exp(j*Kz*z) --
     * exactly the tone the matched filter is built from. That tests
     * conditioning and noise, and says nothing about what happens when the
     * ground returns a signal of a different SHAPE. These six do.
     *
     * Two are documented blind spots and are asserted to FAIL, so that a later
     * change which fixes either one is noticed rather than absorbed silently.
     * From github.com/Hassanforeman/subsurface-sar-tomo's steering stress
     * tests; the depths and verdicts below are this project's own measurements.
     * ------------------------------------------------------------------ */
    RS_CASE("signals whose shape departs from the steering model");
    {
        const size_t nl = 64;
        rs_tomo_params_t q = p;
        q.depth_max = rs_tomo_max_depth(&p, nl);
        const double z = 0.4 * q.depth_max;

        const struct { const char *name; rs_shape_fn fn; } cases[] = {
            { "1. matched tone (positive control)", shape_matched      },
            { "2. moderately damped resonance",     shape_damped       },
            { "3. chirped / dispersive",            shape_chirped      },
            { "4. quadratic phase",                 shape_quadratic    },
            { "5. HEAVY sub-cycle damping",         shape_damped_heavy },
            { "6. very shallow / sub-resolution",   shape_subres       },
        };
        const size_t n_cases = sizeof cases / sizeof cases[0];
        double got[6];

        printf("    true depth %.2f m, %zu looks, resolution %.3f m\n", z, nl, dT);
        for (size_t i = 0; i < n_cases; i++) {
            rs_microm_t m;
            RS_CHECK_OK(make_shaped_scene(&m, 4, nl, &q, z, cases[i].fn, 0.0));
            rs_tomo_t t;
            RS_CHECK_OK(rs_tomo_focus(&m, NULL, &q, NULL, &t));

            double real_c = 0.0, nm = 0.0;
            RS_CHECK_OK(rs_tomo_alignment_null(&t, 100, 3u, &real_c, &nm,
                                               NULL, NULL, NULL));
            got[i] = peak_depth(&t);
            printf("    %-36s z_rec %6.2f m  (err %6.2f)  contrast %7.1f / %.1f\n",
                   cases[i].name, got[i], got[i] - z, real_c, nm);

            rs_tomo_free(&t);
            rs_microm_free(&m);
        }

        /* Recovered. Amplitude taper along the aperture -- cases 2 and 5 -- and
         * a quadratic phase error -- case 4 -- all leave the tone's frequency
         * intact, and frequency is what the matched filter reads. */
        RS_CHECK_NEAR(got[0], z, dT);
        RS_CHECK_NEAR(got[1], z, 2.0 * dT);
        RS_CHECK_NEAR(got[3], z, 2.0 * dT);
        RS_CHECK_NEAR(got[4], z, 2.0 * dT);

        /* BLIND SPOT: a chirped return has no single tone, and the inversion
         * reports a confident depth roughly 60% too deep rather than declining.
         * Asserted to miss, so that a change which fixes it is noticed. */
        RS_CHECK(fabs(got[2] - z) > 2.0 * dT);

        /* BLIND SPOT, and the worse of the two because its contrast is the
         * HIGHEST of all six: a sub-resolution feature is dragged to the
         * surface and reported with more confidence than the correct answer. */
        RS_CHECK(got[5] < 0.25 * z);
    }

    /* One case disagreed with the source, and the disagreement is informative
     * rather than a discrepancy to reconcile away. */
    RS_CASE("heavy damping does NOT reproduce as a blind spot here");
    {
        /* Case 5 recovers the depth correctly, where the source reports it
         * MISSED. Two explanations were proposed and both were measured and
         * rejected, which is why this case exists rather than a footnote.
         *
         * NOT the look count. The source uses 11 looks against 64 here, but 12
         * looks recovers it too, to 0.06 m.
         *
         * NOT the noise floor. exp(-4 i/n) leaves the last looks at 0.018 of
         * the first and the source adds noise of 0.02, so its damped tail is
         * plausibly buried -- but the sweep below recovers the depth at every
         * level up to 0.2, ten times the source's.
         *
         * The remaining structural difference is how the analytic signal is
         * formed. This fixture synthesises it directly, putting exp(j*Kz*z)
         * into the two shift components the tracker would have reported. The
         * source injects a real cosine and recovers the analytic signal by
         * Hilbert transform, and that step assumes an envelope varying slowly
         * against the carrier -- which is exactly what heavy damping violates.
         * If that is the cause, the blind spot belongs to the signal
         * construction rather than to the matched filter, and this project does
         * not inherit it. Stated as the surviving hypothesis, not as a result:
         * it has not been tested here, and testing it means building the
         * Hilbert path this code does not have. */
        const size_t nl = 64;
        rs_tomo_params_t q = p;
        q.depth_max = rs_tomo_max_depth(&p, nl);
        const double z = 0.4 * q.depth_max;
        const double levels[] = { 0.0, 0.002, 0.02, 0.2 };

        for (size_t i = 0; i < sizeof levels / sizeof levels[0]; i++) {
            rs_microm_t m;
            RS_CHECK_OK(make_shaped_scene(&m, 4, nl, &q, z,
                                          shape_damped_heavy, levels[i]));
            rs_tomo_t t;
            RS_CHECK_OK(rs_tomo_focus(&m, NULL, &q, NULL, &t));
            const double got = peak_depth(&t);
            printf("    noise %.3f (tail is %.3f): true %.2f m -> %.2f m, err %+.2f\n",
                   levels[i], exp(-4.0), z, got, got - z);
            RS_CHECK(isfinite(got));
            rs_tomo_free(&t);
            rs_microm_free(&m);
        }

        /* No threshold is asserted. Where the crossover falls depends on the
         * damping constant, the look count and the noise together, and pinning
         * one triple of numbers as the boundary would encode an arbitrary
         * choice as a requirement. The sweep is printed so a change is visible
         * in the ctest log. */
    }

    RS_CASE("grid extent does NOT separate them -- kept as a negative result");
    {
        /* A reflector 1.4x beyond the unambiguous depth at 64 looks folds to
         * 0.4x of it. That folded position is fixed by the wavenumber
         * sampling, not by how much depth the grid covers, so widening the
         * grid moves it not at all -- and the artefact looks exactly as
         * "stable" as a real reflector would. */
        const size_t n_looks_b = 64;
        const double zmax64 = rs_tomo_max_depth(&p, n_looks_b);
        const double z_folded = 1.4 * zmax64;
        const double widths[] = { 0.6, 0.8, 1.0 };
        const size_t n_w = sizeof widths / sizeof widths[0];

        double rec[3];
        for (size_t i = 0; i < n_w; i++) {
            rs_tomo_params_t q = p;
            q.depth_max = widths[i] * zmax64;

            const double amp = 1.0;
            rs_microm_t m;
            RS_CHECK_OK(make_depth_scene(&m, 4, n_looks_b, &q, &z_folded, &amp, 1));
            rs_tomo_t t;
            RS_CHECK_OK(rs_tomo_focus(&m, NULL, &q, NULL, &t));
            rec[i] = peak_depth(&t);
            rs_tomo_free(&t);
            rs_microm_free(&m);
        }

        double mean = 0.0;
        for (size_t i = 0; i < n_w; i++) mean += rec[i];
        mean /= (double)n_w;
        double sd = 0.0;
        for (size_t i = 0; i < n_w; i++) sd += (rec[i] - mean) * (rec[i] - mean);
        sd = sqrt(sd / (double)n_w);

        printf("    folded reflector (z = %.2f m, folds to %.2f m) across grid widths:",
               z_folded, z_folded - zmax64);
        for (size_t i = 0; i < n_w; i++) printf(" %6.2f", rec[i]);
        printf("   sd %.3f m\n", sd);

        /* THE ASSERTION IS THAT THE PROBE FAILS. If this ever starts
         * discriminating, the reason needs understanding before anyone treats
         * grid extent as a guard -- which is exactly the mistake this case
         * exists to prevent. */
        RS_CHECK(sd < dT);
    }

    RS_TEST_END();
}
