/* Pre-flight validation of a collect against a requested measurement.
 * See include/resonarsat/validate.h for what this is and is not. */

#include "resonarsat/validate.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/* The observation ratio at which the sub-look begins to resolve the target's
 * paired-echo train, measured rather than assumed.
 *
 * A vibrating target images as a train of echoes spaced f*lambda*R/(2*v_p), and
 * a sub-look resolves them when that exceeds its own azimuth resolution
 * lambda*R/(2*v_p*t_sap). The ratio of the two is exactly f*t_sap -- the
 * observation ratio. Past it the tracked feature fragments into comparably
 * bright spots and the correlation argmax hops between them.
 *
 * The bracket is from measurement: eta 0.078 and 0.141 recover an injected
 * frequency, eta 0.501 fails at any amplitude. Between 0.14 and 0.50 is
 * unmeasured, so WARN spans it rather than either bound being asserted.
 * Vattulainen et al. operate at 0.39-0.69 and report performance degrading
 * across that range, which is consistent with the upper bound being real and
 * with 0.5 being roughly where it bites. */
#define RS_ETA_GOOD 0.20
#define RS_ETA_BAD  0.50

/* The aperture fraction range the published validation works at, from Lotti et
 * al. (SHMII-13), whose three tests sit at 7.6%, 4.5% and 4.5%. Below it their
 * text warns the target "may no longer appear as a distinct feature". */
#define RS_ALPHA_LO 0.045
#define RS_ALPHA_HI 0.076

/* The quantisation floor of a correlation shift, in pixels, for a given
 * sub-pixel refinement factor. Matches what the spectrum stage reports:
 * 1/40 px gives 0.061225 and 1/10 px gives 0.2449. */
static double rs_validate_floor_px(size_t upsample)
{
    return (upsample > 0) ? 2.449 / (double)upsample : 1.0;
}

/* Record one finding, returning its level so the caller can track the worst. */
static rs_validate_level_t rs_v_add(rs_validate_finding_t *out, size_t *n,
                                    rs_validate_check_t check,
                                    rs_validate_level_t level,
                                    const char *name, const char *fmt, ...)
{
    va_list ap;
    rs_validate_finding_t *f = &out[*n];
    f->check = check;
    f->level = level;
    f->name = name;
    va_start(ap, fmt);
    vsnprintf(f->detail, sizeof f->detail, fmt, ap);
    va_end(ap);
    (*n)++;
    return level;
}

/* The one-word name of a level, for printing. */
const char *rs_validate_level_str(rs_validate_level_t level)
{
    switch (level) {
    case RS_V_PASS:    return "PASS";
    case RS_V_WARN:    return "WARN";
    case RS_V_FAIL:    return "FAIL";
    case RS_V_UNKNOWN: return "UNKNOWN";
    }
    return "?";
}

/* Fill a request with the published validation's operating point. */
void rs_validate_req_default(rs_validate_req_t *req)
{
    if (!req) return;
    memset(req, 0, sizeof *req);
    req->alpha = 0.05;
    req->overlap = 0.40;
    req->upsample = 40;
    req->cell_m = 1.0;
    req->grid_n = 512;
}

/* Run every check and write the findings. */
rs_validate_level_t rs_validate(const rs_validate_req_t *req,
                                rs_validate_finding_t *out, size_t *n_out)
{
    if (!req || !out || !n_out) return RS_V_FAIL;
    size_t n = 0;
    rs_validate_level_t worst = RS_V_PASS;
#define WORST(l) do { const rs_validate_level_t _l = (l); \
    if (_l == RS_V_FAIL || (_l == RS_V_WARN && worst != RS_V_FAIL)) worst = _l; } while (0)

    /* ---- what the pipeline cannot run without at all ---------------------- */
    if (!(req->dwell_s > 0.0) || !(req->lambda_m > 0.0) ||
        !(req->slant_range_m > 0.0) || !(req->v_platform_ms > 0.0)) {
        WORST(rs_v_add(out, &n, RS_VALIDATE_METADATA, RS_V_FAIL, "metadata",
              "missing or non-positive geometry: dwell %g s, lambda %g m, "
              "slant range %g m, platform speed %g m/s. Nothing downstream can "
              "run without all four.",
              req->dwell_s, req->lambda_m, req->slant_range_m, req->v_platform_ms));
        *n_out = n;
        return RS_V_FAIL;
    }
    WORST(rs_v_add(out, &n, RS_VALIDATE_METADATA, RS_V_PASS, "metadata",
          "dwell %.3f s, lambda %.4f m, slant range %.1f km, platform %.0f m/s, "
          "incidence %.1f deg, %zu pulses",
          req->dwell_s, req->lambda_m, req->slant_range_m / 1000.0,
          req->v_platform_ms, req->incidence_rad * 180.0 / M_PI, req->n_pulse));

    const double t_sap = req->alpha * req->dwell_s;
    const double f = req->target_freq_hz;
    const double df = 1.0 / req->dwell_s;

    /* ---- can the dwell resolve the frequency at all ----------------------- */
    if (f > 0.0) {
        const double bins = f / df;
        WORST(rs_v_add(out, &n, RS_VALIDATE_FREQ_RESOLUTION,
              (bins >= 3.0) ? RS_V_PASS : RS_V_FAIL, "frequency resolution",
              "%.4f Hz resolution from a %.3f s dwell; %.3f Hz is %.1f bins. "
              "%s", df, req->dwell_s, f, bins,
              (bins >= 3.0) ? "Resolvable."
                            : "Below three bins, so indistinguishable from drift."));
    }

    /* ---- Nyquist of the sub-aperture series ------------------------------- */
    {
        const double denom = (req->alpha > 0.0)
                           ? (1.0 / req->alpha) - req->overlap : 0.0;
        const double n_looks = (req->overlap < 1.0 && denom > 0.0)
                             ? denom / (1.0 - req->overlap) : 0.0;
        const double d = (n_looks > 1.0)
                       ? (req->dwell_s - t_sap) / (n_looks - 1.0) : req->dwell_s;
        const double f_max = 1.0 / (2.0 * d);
        WORST(rs_v_add(out, &n, RS_VALIDATE_BAND,
              (f <= 0.0 || f < f_max) ? RS_V_PASS : RS_V_FAIL, "observable band",
              "alpha %.3f%% and overlap %.2f give %.0f sub-apertures, dt %.4f s, "
              "so f_max is %.3f Hz.%s",
              100.0 * req->alpha, req->overlap, n_looks, d, f_max,
              (f > 0.0 && f >= f_max) ? " The target is above it and will alias."
                                      : ""));

        /* The spectral route splits a focused image and needs twice as many
         * azimuth lines as looks. The pulse route has no such constraint. */
        const size_t need = (size_t)(2.0 * n_looks + 0.5);
        WORST(rs_v_add(out, &n, RS_VALIDATE_GRID,
              (req->grid_n >= need) ? RS_V_PASS : RS_V_FAIL, "grid width",
              "the uniform spectral route needs n_az >= 2*n_looks = %zu; "
              "--size %zu %s. The pulse route is exempt.",
              need, req->grid_n,
              (req->grid_n >= need) ? "is enough" : "is NOT enough"));
    }

    /* ---- the observation ratio, which is the ghost-resolution limit ------- */
    if (f > 0.0) {
        const double eta = f * t_sap;
        const rs_validate_level_t lvl = (eta <= RS_ETA_GOOD) ? RS_V_PASS
                                      : (eta >= RS_ETA_BAD)  ? RS_V_FAIL
                                                             : RS_V_WARN;
        const double f_ok = (t_sap > 0.0) ? RS_ETA_GOOD / t_sap : 0.0;
        WORST(rs_v_add(out, &n, RS_VALIDATE_OBSERVATION_RATIO, lvl,
              "observation ratio",
              "t_sap %.4f s at %.3f Hz gives eta %.3f. %s "
              "At this aperture fraction eta stays under %.2f only below "
              "%.3f Hz.",
              t_sap, f, eta,
              (lvl == RS_V_PASS)
                ? "The sub-look does not resolve the target's paired echoes."
                : (lvl == RS_V_WARN)
                ? "Between the measured working range and the measured failure; "
                  "untested here and where the published validation sits."
                : "The sub-look RESOLVES the paired echoes, the tracked feature "
                  "fragments, and correlation tracking fails at any amplitude.",
              RS_ETA_GOOD, f_ok));
    }

    /* ---- displacement averaging vanishes at integer eta ------------------- */
    if (f > 0.0 && t_sap > 0.0) {
        const double eta = f * t_sap;
        const double k = floor(eta + 0.5);
        const double dist = (k >= 1.0) ? fabs(eta - k) : eta;
        const rs_validate_level_t lvl = (k >= 1.0 && dist < 0.10) ? RS_V_FAIL
                                      : (k >= 1.0 && dist < 0.20) ? RS_V_WARN
                                                                  : RS_V_PASS;
        WORST(rs_v_add(out, &n, RS_VALIDATE_AVERAGING_NULL, lvl,
              "averaging nulls",
              "the response vanishes at integer eta, so at %.3f, %.3f, %.3f Hz "
              "for this t_sap. The target sits %.3f from the nearest.%s",
              1.0 / t_sap, 2.0 / t_sap, 3.0 / t_sap, dist,
              (lvl == RS_V_FAIL)
                ? " That is on a null: a displacement-averaging observable has "
                  "no response there." : ""));
    }

    /* ---- aperture fraction against the range anyone has validated --------- */
    {
        const int inside = (req->alpha >= RS_ALPHA_LO && req->alpha <= RS_ALPHA_HI);
        WORST(rs_v_add(out, &n, RS_VALIDATE_APERTURE_FRACTION,
              inside ? RS_V_PASS : RS_V_WARN, "aperture fraction",
              "alpha %.3f%%; published validation sits at %.1f-%.1f%%. %s",
              100.0 * req->alpha, 100.0 * RS_ALPHA_LO, 100.0 * RS_ALPHA_HI,
              inside ? "Inside."
                     : (req->alpha < RS_ALPHA_LO)
                       ? "Below it: the target may stop being a distinct feature "
                         "in each sub-look."
                       : "Above it: longer sub-apertures integrate more of the "
                         "motion away."));
    }

    /* ---- the timing the sub-aperture stage assumes is uniform ------------- */
    if (req->prf_min_hz > 0.0 && req->prf_max_hz > 0.0) {
        const double spread = 100.0 * (req->prf_max_hz - req->prf_min_hz)
                            / req->prf_hz;
        const double gap_intervals = (req->prf_hz > 0.0)
                                   ? req->worst_gap_s * req->prf_hz : 0.0;
        const rs_validate_level_t lvl = (gap_intervals > 20.0 || spread > 5.0)
                                      ? RS_V_WARN : RS_V_PASS;
        WORST(rs_v_add(out, &n, RS_VALIDATE_PRF_STABILITY, lvl, "PRF stability",
              "instantaneous PRF spans %.2f%% of nominal; largest gap is %.1f "
              "pulse intervals. The sub-aperture stage lays centre times on a "
              "uniform grid, so this enters as a time-base distortion.%s",
              spread, gap_intervals,
              (lvl == RS_V_WARN) ? " Large enough to check before trusting a "
                                   "frequency." : ""));
    }

    /* ---- the smallest motion the tracker could see ------------------------ */
    if (f > 0.0 && req->cell_m > 0.0 && req->upsample > 0) {
        const double proj = cos(req->incidence_rad);
        const double floor_px = rs_validate_floor_px(req->upsample);
        /* A vertical displacement A at frequency f gives a radial velocity
         * 2*pi*f*A*proj, which a coregistrator sees as an azimuth shift of
         * R*v_r/V metres. Invert for the amplitude at the floor. */
        const double denom = req->slant_range_m * 2.0 * M_PI * f * proj;
        const double a_min = (denom > 0.0)
            ? floor_px * req->cell_m * req->v_platform_ms / denom : INFINITY;
        rs_validate_level_t lvl = RS_V_UNKNOWN;
        char verdict[220];
        if (req->target_amp_m > 0.0) {
            const double margin = req->target_amp_m / a_min;
            lvl = (margin >= 3.0) ? RS_V_PASS : (margin >= 1.0) ? RS_V_WARN
                                                                : RS_V_FAIL;
            snprintf(verdict, sizeof verdict,
                     "The %.2f mm asked for is %.1fx that. %s",
                     1000.0 * req->target_amp_m, margin,
                     (lvl == RS_V_PASS) ? "Comfortable."
                     : (lvl == RS_V_WARN) ? "Within a factor of three of the "
                                            "floor; marginal."
                                          : "Below the floor and undetectable.");
        } else {
            snprintf(verdict, sizeof verdict,
                     "No target amplitude given, so this is the threshold "
                     "rather than a verdict.");
        }
        WORST(rs_v_add(out, &n, RS_VALIDATE_SENSITIVITY, lvl, "sensitivity",
              "at 1/%zu px refinement the floor is %.4f px, which at %.2f m "
              "cells is a vertical amplitude of %.3f mm at %.3f Hz. %s",
              req->upsample, floor_px, req->cell_m, 1000.0 * a_min, f, verdict));
    }

    /* ---- what reading it will cost --------------------------------------- */
    if (req->n_pulse > 0 && req->n_rbin > 0) {
        const double gb = (double)req->n_pulse * (double)req->n_rbin * 8.0 / 1e9;
        WORST(rs_v_add(out, &n, RS_VALIDATE_MEMORY, RS_V_PASS, "memory",
              "%zu pulses x %zu range bins is %.1f GB as complex float. "
              "Machines under that will need a narrower --rbins.",
              req->n_pulse, req->n_rbin, gb));
    }

    /* ---- the question the file cannot answer ------------------------------ */
    rs_v_add(out, &n, RS_VALIDATE_GROUND_TRUTH, RS_V_UNKNOWN, "ground truth",
             "Whether anything in this scene moves is not a property of the "
             "collect, and no check above establishes it. Every validated "
             "result in this literature used a corner reflector on a shaker "
             "with a synchronous displacement sensor; docs/DATASETS.md records "
             "that no such collect is in any open archive. A pass here means "
             "the configuration is capable, not that the measurement is real.");

    *n_out = n;
    return worst;
#undef WORST
}
