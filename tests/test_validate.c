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
    RS_CASE("the Khufu configuration is refused, for the reason it failed");
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

        const rs_validate_level_t worst = rs_validate(&r, f, &n);
        const rs_validate_finding_t *eta = find(f, n, RS_VALIDATE_OBSERVATION_RATIO);
        RS_CHECK(worst == RS_V_FAIL);
        RS_CHECK(eta != NULL && eta->level == RS_V_FAIL);
        printf("    eta check: %s\n", eta ? eta->detail : "(missing)");

        /* And the aperture fraction, which the run was corrected TO, passes --
         * so the verdict is not just everything failing at once. */
        const rs_validate_finding_t *a = find(f, n, RS_VALIDATE_APERTURE_FRACTION);
        RS_CHECK(a != NULL && a->level == RS_V_PASS);
    }

    /* ------------------------------------------------------------------
     * The same collect asked for a frequency it can actually reach.
     * ------------------------------------------------------------------ */
    RS_CASE("a frequency inside the reachable band passes");
    {
        rs_validate_req_t r;
        giza(&r);
        r.target_freq_hz = 0.10;   /* below the 0.122 Hz eta limit at alpha 5% */
        r.target_amp_m = 0.005;
        r.alpha = 0.05;
        r.overlap = 0.88;
        r.upsample = 40;
        r.cell_m = 1.0;
        r.grid_n = 512;

        const rs_validate_level_t worst = rs_validate(&r, f, &n);
        const rs_validate_finding_t *eta = find(f, n, RS_VALIDATE_OBSERVATION_RATIO);
        RS_CHECK(eta != NULL && eta->level == RS_V_PASS);
        /* UNKNOWN from ground truth is not a failure; FAIL would be. */
        RS_CHECK(worst != RS_V_FAIL);
        printf("    eta check: %s\n", eta ? eta->detail : "(missing)");
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
