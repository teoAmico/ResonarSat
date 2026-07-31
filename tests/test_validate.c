/* The pre-flight validator, checked against configurations whose outcome is
 * already known from measurement.
 *
 * Every threshold in rs_validate() came from a run recorded in runs/, so the
 * cases below are not hypotheticals: each one asserts that the validator would
 * have called, in advance and in milliseconds, what a real run took twenty
 * minutes to discover. If a threshold is ever changed, these fail and the
 * evidence for the change has to be produced. */

#include "resonarsat/validate.h"
#include "rs_test.h"

#include <math.h>
#include <string.h>

/* Find one check's finding, or NULL. */
static const rs_validate_finding_t *find(const rs_validate_finding_t *f, size_t n,
                                         rs_validate_check_t which)
{
    for (size_t i = 0; i < n; i++) if (f[i].check == which) return &f[i];
    return NULL;
}

/* The Giza collect's geometry, as the reader derives it. */
static void giza(rs_validate_req_t *r)
{
    rs_validate_req_default(r);
    r->dwell_s       = 32.869;
    r->prf_hz        = 10196.35;
    r->lambda_m      = 0.0322;
    r->slant_range_m = 762800.0;
    r->v_platform_ms = 7263.0;
    r->incidence_rad = 39.5 * M_PI / 180.0;
    r->n_pulse       = 335141;
    r->n_rbin        = 4096;
    r->prf_min_hz    = 10121.0;
    r->prf_max_hz    = 10235.0;
    r->worst_gap_s   = 9.0 / 10196.35;
}

int main(void)
{
    rs_validate_finding_t f[RS_VALIDATE_N_CHECKS];
    size_t n = 0;

    /* ------------------------------------------------------------------
     * The configuration runs/giza/2026-07-30-validated-spot-khufu actually
     * used. It produced a complete, plausible result and was later shown to
     * be incapable of the measurement. The validator must refuse it.
     * ------------------------------------------------------------------ */
    RS_CASE("the Khufu configuration is flagged on the observation ratio");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 2.0;
        r.target_amp_m = 0.001;
        r.alpha = 0.05;
        r.overlap = 0.88;
        r.upsample = 40;
        r.cell_m = 1.0;
        r.grid_n = 512;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *eta = find(f, n, RS_VALIDATE_OBSERVATION_RATIO);
        /* WARN, not FAIL: the ghost-resolution mechanism is confirmed by
         * tests/test_pairedecho.c, but the threshold at which the tracker
         * actually breaks is not measured -- see rs_validate()'s comment. */
        RS_CHECK(eta != NULL && eta->level == RS_V_WARN);
        printf("    eta check: %s\n", eta ? eta->detail : "(missing)");

        /* And the aperture fraction, which the run was corrected TO, passes --
         * so the verdict is not just everything failing at once. */
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_APERTURE_FRACTION);
        RS_CHECK(a != NULL && a->level == RS_V_PASS);
    }

    /* ------------------------------------------------------------------
     * Lowering the frequency does NOT rescue that configuration, and this is
     * the case that says so. The wrap ceiling is set by the sub-look
     * resolution and the cell -- neither of which depends on the target -- so
     * at alpha 5% over a 32.9 s dwell the Giza collect has no admissible
     * amplitude at ANY frequency. What has to change is the sub-aperture.
     * ------------------------------------------------------------------ */
    RS_CASE("a lower frequency does not rescue an unusable sub-aperture");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 0.10;
        r.target_amp_m = 0.005;
        r.alpha = 0.05;
        r.overlap = 0.88;
        r.upsample = 40;
        r.cell_m = 1.0;
        r.grid_n = 512;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *eta = find(f, n, RS_VALIDATE_OBSERVATION_RATIO);
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_AMBIGUITY);
        RS_CHECK(eta != NULL && eta->level == RS_V_PASS);   /* eta is fine here */
        RS_CHECK(a != NULL && a->level == RS_V_FAIL);       /* and it does not matter */
    }

    /* Shortening the sub-aperture is what fixes it -- and doing so leaves the
     * aperture fraction BELOW the range the published campaigns validated.
     * Those campaigns run on far shorter dwells, so what transfers between
     * collects is t_sap in seconds, not alpha as a fraction. The two checks
     * disagree here on purpose; the tool reports both. */
    RS_CASE("a short sub-aperture opens a window, against the alpha rule");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 0.10;
        r.target_amp_m = 0.0;
        r.alpha = 0.008;           /* t_sap 0.26 s, well under the 4.5% floor */
        r.overlap = 0.0;
        r.upsample = 40;
        r.cell_m = 1.0;
        r.grid_n = 512;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_AMBIGUITY);
        const rs_validate_finding_t *al = find(f, n, RS_VALIDATE_APERTURE_FRACTION);
        RS_CHECK(a != NULL && a->level == RS_V_PASS);
        RS_CHECK(al != NULL && al->level != RS_V_PASS);
        printf("    ambiguity: %s\n", a->detail);
        printf("    alpha:     %s\n", al->detail);
    }

    /* ------------------------------------------------------------------
     * A target sitting on an integer observation ratio, where a
     * displacement-averaging observable has no response. This is how the
     * first positive control was accidentally built, at 1.0 Hz against a
     * 1.002 s sub-aperture.
     * ------------------------------------------------------------------ */
    RS_CASE("a target on an averaging null is caught");
    {
        rs_validate_req_t r;
        rs_validate_req_default(&r);
        r.dwell_s = 20.0;
        r.lambda_m = 0.031228;
        r.slant_range_m = 610328.0;
        r.v_platform_ms = 7500.0;
        r.incidence_rad = 35.0 * M_PI / 180.0;
        r.n_pulse = 8000;
        r.alpha = 0.0501;          /* t_sap = 1.002 s */
        r.overlap = 0.88;
        r.target_freq_hz = 1.0;    /* eta = 1.002, exactly a null */
        r.grid_n = 512;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *nul = find(f, n, RS_VALIDATE_AVERAGING_NULL);
        RS_CHECK(nul != NULL && nul->level == RS_V_FAIL);
        printf("    null check: %s\n", nul ? nul->detail : "(missing)");
    }

    /* ------------------------------------------------------------------
     * The two-sided excursion window, which is what actually disqualifies the
     * Khufu configuration. Both bounds were measured; see RS_VALIDATE_AMBIGUITY.
     * ------------------------------------------------------------------ */
    RS_CASE("a configuration with no admissible amplitude is refused");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 2.0;
        r.target_amp_m = 0.001;
        r.alpha = 0.05;            /* t_sap 1.64 s -> sub-look resolution 1.03 m */
        r.overlap = 0.88;
        r.upsample = 40;
        r.cell_m = 1.0;
        r.grid_n = 512;

        const rs_validate_level_t worst = rs_validate(&r, f, &n);
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_AMBIGUITY);
        RS_CHECK(a != NULL && a->level == RS_V_FAIL);
        RS_CHECK(worst == RS_V_FAIL);
        RS_CHECK(strstr(a->detail, "BELOW the floor") != NULL);
        printf("    %s\n", a->detail);
    }

    /* The shape that recovers an injected frequency 7 times out of 7, as the
     * fixed-displacement sweep in POSITIVE-CONTROL.md found: the simulator
     * geometry, pulse route, 128 looks, zero overlap. */
    RS_CASE("the operating point that works has a usable window");
    {
        rs_validate_req_t r;
        rs_validate_req_default(&r);
        r.dwell_s = 20.0;
        r.lambda_m = 0.031228;
        r.slant_range_m = 610328.0;
        r.v_platform_ms = 7500.0;
        r.incidence_rad = 34.99 * M_PI / 180.0;
        r.n_pulse = 8000;
        r.alpha = 1.0 / 128.0;     /* 128 looks, zero overlap */
        r.overlap = 0.0;
        r.target_freq_hz = 0.5;
        r.target_amp_m = 0.00955;  /* the amplitude that recovered 0.504 Hz */
        r.cell_m = 0.4;
        r.upsample = 10;
        r.grid_n = 320;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_AMBIGUITY);
        RS_CHECK(a != NULL && a->level == RS_V_PASS);
        printf("    %s\n", a->detail);

        /* And the amplitude that wrapped into the third harmonic is refused. */
        r.target_amp_m = 0.0382;   /* 4x, which reported 1.512 for 0.504 Hz */
        rs_validate(&r, f, &n);
        a = find(f, n, RS_VALIDATE_AMBIGUITY);
        RS_CHECK(a != NULL && a->level == RS_V_FAIL);
        RS_CHECK(strstr(a->detail, "would wrap") != NULL);
    }

    /* A marginally open window is not a usable one. This configuration has a
     * 4.8 px ceiling, which clears the 4 px floor a textured scene gives -- and
     * targets placed at 4.0, 4.4 and 4.8 px inside that nominal window all
     * missed. The conservative floor is what makes the validator say so. */
    RS_CASE("a marginally open window is still refused");
    {
        rs_validate_req_t r;
        rs_validate_req_default(&r);
        r.dwell_s = 20.0;
        r.lambda_m = 0.031228;
        r.slant_range_m = 610328.0;
        r.v_platform_ms = 7500.0;
        r.incidence_rad = 34.99 * M_PI / 180.0;
        r.n_pulse = 8000;
        r.alpha = 1.0 / 19.96;     /* 159 looks at 0.88 overlap */
        r.overlap = 0.88;
        r.target_freq_hz = 0.5;
        r.cell_m = 0.4;
        r.upsample = 10;
        r.grid_n = 320;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_AMBIGUITY);
        RS_CHECK(a != NULL && a->level == RS_V_FAIL);
        printf("    %s\n", a->detail);
    }

    /* The frequency a motionless scene produced, against the configuration
     * that produced it. The step-based f_max is 4.16 Hz and would call this
     * comfortably in band; the averaging response puts it out of reach. */
    RS_CASE("a frequency past the sub-aperture's averaging response is refused");
    {
        rs_validate_req_t r;
        rs_validate_req_default(&r);
        r.dwell_s = 20.0;
        r.lambda_m = 0.031228;
        r.slant_range_m = 610328.0;
        r.v_platform_ms = 7500.0;
        r.incidence_rad = 34.99 * M_PI / 180.0;
        r.n_pulse = 8000;
        r.alpha = 1.0 / 19.96;       /* t_sap 1.002 s -> band 0.499 Hz */
        r.overlap = 0.88;            /* dt 0.120 s -> step would say 4.16 Hz */
        r.cell_m = 0.4;
        r.upsample = 10;
        r.grid_n = 320;

        r.target_freq_hz = 1.569;    /* what the static scene reported */
        rs_validate(&r, f, &n);
        const rs_validate_finding_t *b = find(f, n, RS_VALIDATE_BAND);
        RS_CHECK(b != NULL && b->level == RS_V_FAIL);
        RS_CHECK(strstr(b->detail, "ABOVE the band") != NULL);
        printf("    %s\n", b->detail);

        /* eta = f*t_sap, so the band edge is exactly eta = 0.5. */
        r.target_freq_hz = 0.40;     /* eta 0.401 */
        rs_validate(&r, f, &n);
        b = find(f, n, RS_VALIDATE_BAND);
        RS_CHECK(b != NULL && b->level == RS_V_PASS);
    }

    /* And the configuration that recovers 7 injected frequencies out of 7 must
     * not have any of them refused -- the check has to be non-trivial. */
    RS_CASE("the working operating point keeps its whole recovered range");
    {
        rs_validate_req_t r;
        rs_validate_req_default(&r);
        r.dwell_s = 20.0;
        r.lambda_m = 0.031228;
        r.slant_range_m = 610328.0;
        r.v_platform_ms = 7500.0;
        r.incidence_rad = 34.99 * M_PI / 180.0;
        r.n_pulse = 8000;
        r.alpha = 1.0 / 128.0;       /* t_sap 0.156 s -> band 3.2 Hz */
        r.overlap = 0.0;
        r.cell_m = 0.4;
        r.upsample = 10;
        r.grid_n = 320;

        const double got[7] = { 0.101, 0.202, 0.302, 0.504, 0.706, 1.008, 1.411 };
        for (int i = 0; i < 7; i++) {
            r.target_freq_hz = got[i];
            rs_validate(&r, f, &n);
            const rs_validate_finding_t *b = find(f, n, RS_VALIDATE_BAND);
            RS_CHECK(b != NULL && b->level == RS_V_PASS);
        }
    }

    /* The coherence mask configuration A actually ran behind. At 0.99 overlap
     * almost every look pair shares band, so the estimator cannot report below
     * 0.574 whatever the scene does -- and the mask was 0.4. */
    RS_CASE("a coherence mask below its own floor is called vacuous");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 0.05;
        r.alpha = 1.0 / 2.27;        /* 128 looks at 0.99 overlap */
        r.overlap = 0.99;
        r.coherence_min = 0.4;
        r.win = 32;
        r.cell_m = 2.0;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *c = find(f, n, RS_VALIDATE_COHERENCE_GATE);
        RS_CHECK(c != NULL && c->level == RS_V_FAIL);
        RS_CHECK(strstr(c->detail, "vacuous") != NULL);
        printf("    %s\n", c->detail);

        /* Configuration B, same collect, 0.9 overlap over 2048 looks: only 1%
         * of pairs share band, and the same 0.4 mask is meaningful. */
        r.alpha = 1.0 / 205.7;
        r.overlap = 0.9;
        rs_validate(&r, f, &n);
        c = find(f, n, RS_VALIDATE_COHERENCE_GATE);
        RS_CHECK(c != NULL && c->level == RS_V_PASS);
        printf("    %s\n", c->detail);
    }

    /* The phase route's floor, which is the reason to prefer it. */
    RS_CASE("the phase floor is far below the correlation floor");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 0.10;
        r.alpha = 0.008;
        r.overlap = 0.0;
        r.coherence_min = 0.4;
        r.win = 32;
        r.cell_m = 1.0;

        rs_validate(&r, f, &n);
        const rs_validate_finding_t *p = find(f, n, RS_VALIDATE_PHASE_FLOOR);
        RS_CHECK(p != NULL && p->level == RS_V_PASS);
        printf("    %s\n", p->detail);
    }

    /* ------------------------------------------------------------------
     * Guards.
     * ------------------------------------------------------------------ */
    RS_CASE("a collect without geometry is refused immediately");
    {
        rs_validate_req_t r;
        rs_validate_req_default(&r);
        r.dwell_s = 0.0;
        const rs_validate_level_t worst = rs_validate(&r, f, &n);
        RS_CHECK(worst == RS_V_FAIL);
        RS_CHECK(n == 1);   /* nothing else is answerable */
    }

    RS_CASE("ground truth is always unknown, whatever else passes");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 0.05;
        r.alpha = 0.05;
        rs_validate(&r, f, &n);
        const rs_validate_finding_t *g = find(f, n, RS_VALIDATE_GROUND_TRUTH);
        RS_CHECK(g != NULL && g->level == RS_V_UNKNOWN);
        RS_CHECK(strstr(g->detail, "not a property of the collect") != NULL);
    }

    /* The sensitivity figure should scale as 1/f: the same displacement at a
     * lower frequency is a smaller velocity and so a smaller azimuth shift. */
    RS_CASE("the detection floor rises as the frequency falls");
    {
        rs_validate_req_t r;
        giza(&r);
        r.alpha = 0.05;
        r.cell_m = 1.0;
        r.upsample = 40;

        double amp[2];
        const double freqs[2] = { 0.05, 0.10 };
        for (int i = 0; i < 2; i++) {
            r.target_freq_hz = freqs[i];
            r.target_amp_m = 0.0;
            rs_validate(&r, f, &n);
            const rs_validate_finding_t *s = find(f, n, RS_VALIDATE_SENSITIVITY);
            RS_CHECK(s != NULL);
            /* Pull the millimetre figure back out of the message. */
            const char *p = s ? strstr(s->detail, "amplitude of ") : NULL;
            amp[i] = p ? atof(p + 13) : 0.0;
            printf("    %.2f Hz -> floor %.4f mm\n", freqs[i], amp[i]);
        }
        RS_CHECK(amp[0] > 0.0 && amp[1] > 0.0);
        RS_CHECK_REL(amp[0], 2.0 * amp[1], 0.02);
    }

    RS_TEST_END();
}
