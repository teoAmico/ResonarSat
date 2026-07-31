/* Pre-flight validation of a collect against a requested measurement.
 * See include/resonarsat/validate.h for what this is and is not. */

#include "resonarsat/validate.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/* The observation ratio at which the sub-look begins to resolve the target's
 * paired-echo train.
 *
 * THE MECHANISM IS REAL; THE THRESHOLD IS NOT YET MEASURED. A vibrating target
 * images as a train of echoes spaced f*lambda*R/(2*v_p) -- confirmed directly by
 * tests/test_pairedecho.c, which finds them at the predicted offsets with the
 * predicted Bessel amplitudes -- and a sub-look resolves that train when the
 * spacing exceeds its own azimuth resolution lambda*R/(2*v_p*t_sap). The ratio
 * of the two is exactly f*t_sap, the observation ratio. That much is arithmetic.
 *
 * What is NOT established is where the tracker actually breaks. The bracket
 * below was inferred by comparing a configuration that recovered an injected
 * frequency against one that did not, and those two differ in more than eta:
 * route, look count, overlap, sub-aperture length, injected displacement and
 * the modulation index B all move together, and eta, displacement and B cannot
 * be varied independently at all -- displacement goes as f*A and B as A, while
 * eta goes as f*t_sap.
 *
 * A ladder run afterwards found something worse. The uniform spectral-split
 * route reports 1.569 Hz at prominence 27.9 on a target that is not moving at
 * all, and reports the same 1.569 Hz for every injected frequency from 0.2 to
 * 1.4 Hz. A fixed spurious line of that strength swamps the comparison the
 * bracket came from, so the numbers below describe an expectation from physics
 * rather than a measured boundary.
 *
 * They are therefore a WARNING and not a refusal. See
 * runs/giza/2026-07-30-validated-spot-khufu/POSITIVE-CONTROL.md. */
#define RS_ETA_GOOD 0.20
#define RS_ETA_BAD  0.50

/* The aperture fraction range the published validation works at, from Lotti et
 * al. (SHMII-13), whose three tests sit at 7.6%, 4.5% and 4.5%. Below it their
 * text warns the target "may no longer appear as a distinct feature". */
/* The peak-to-peak azimuth excursion, in tracking-grid pixels, below which the
 * correlation tracker reports a fixed spurious frequency rather than the
 * target's. Measured on the pulse route at 128 looks and zero overlap over
 * 0.4 m cells, sweeping amplitude at a fixed 0.5 Hz: an isolated point misses
 * at 3, 4, 5 and 6 px and recovers from 7 px up.
 *
 * The floor DEPENDS ON THE SCENE. Repeating the sweep against a coherently
 * vibrating clutter patch -- which is what sim_cphd's help says to use when the
 * question is about the tracker rather than about focusing -- moves it to 4 px.
 * The conservative value is kept, for a measured reason: a configuration whose
 * window is only marginally open does not work. The uniform route at 159 looks
 * and 0.88 overlap has a 4.8 px ceiling, so against a 4 px textured floor it
 * has a nominal 4.0-4.8 px window; targets placed at 4.0, 4.4 and 4.8 px inside
 * that window all miss. A ceiling must clear the floor by a real margin, not
 * merely exceed it, and where that margin lies is not measured.
 *
 * HOW IT SCALES IS NOT KNOWN, AND THE PIXEL READING IS REFUTED. Expressing the
 * floor in pixels implies it shrinks in metres as the cell shrinks, so a finer
 * grid would buy sensitivity. It does not. Holding the physical scene and
 * window fixed at 128 m and 12.8 m and halving the cell from 0.4 m to 0.2 m,
 * a 1.2 m peak-to-peak excursion -- 6 px at the finer cell, comfortably above
 * a 4 px floor -- misses at BOTH cells. Two clean misses where the pixel
 * reading predicts a recovery.
 *
 * The obvious replacement, a fixed fraction of the sub-look resolution cell,
 * is refuted by different data: it makes the ceiling and floor both
 * proportional to that cell, so the window is always open at a fixed 7.5x
 * ratio, while the uniform route at 159 looks and 0.88 overlap was measured
 * CLOSED -- targets at 4.0, 4.4 and 4.8 px all missing inside the window it
 * predicts.
 *
 * So the number below holds at the geometry and cell it was measured at, and
 * the figures this file derives from it at other cells are NOT established.
 * Treat them as order-of-magnitude. Settling this needs a designed experiment
 * over cell, sub-look resolution and window size in pixels, which are confounded
 * in every run made so far. See runs/giza/2026-07-30-validated-spot-khufu/. */
#define RS_TRACK_FLOOR_PX 7.0

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
    req->win = 32;
    req->coherence_min = 0.40;   /* rs_microm_params_default() */
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

    /* ---- the band the sub-apertures can actually carry -------------------
     *
     * TWO DIFFERENT LIMITS, AND THE LOOSER ONE IS NOT THE REAL ONE. The step
     * between sub-apertures sets how finely the series is sampled, giving
     * 1/(2*dt). But each sub-aperture AVERAGES the target's motion over its own
     * duration, which is a lowpass whose first null is at 1/t_sap -- the same
     * null RS_VALIDATE_AVERAGING_NULL reports. Overlap shrinks dt without
     * shortening t_sap, so it buys finer sampling of the same band, not a wider
     * band. Quoting 1/(2*dt) at 88% overlap overstates the reach by a factor of
     * 1/(1-overlap), which is 8.3.
     *
     * This is not a theoretical worry. On a scene with NO motion at all the
     * uniform route at 159 looks and 0.88 overlap reports 1.569 Hz at
     * prominence 27.9, and reports the same 1.569 Hz for every injected
     * frequency from 0.2 to 1.4 Hz. Its dt-based f_max is 4.16 Hz, so that
     * number looks comfortably in band. Its t_sap is 1.002 s, so the first
     * averaging null is at 0.998 Hz and the artefact sits PAST it, where the
     * sub-aperture has essentially no response. Nothing measured can live
     * there.
     *
     * Note what this makes of the observation ratio: eta = f*t_sap, so
     * f < 1/(2*t_sap) is exactly eta < 0.5. The 0.5 bracket withdrawn from
     * RS_VALIDATE_OBSERVATION_RATIO was the right number reached by the wrong
     * argument -- it is not about resolving paired echoes, it is the
     * sub-aperture's own averaging response, and it is derived rather than
     * measured. */
    {
        const double denom = (req->alpha > 0.0)
                           ? (1.0 / req->alpha) - req->overlap : 0.0;
        const double n_looks = (req->overlap < 1.0 && denom > 0.0)
                             ? denom / (1.0 - req->overlap) : 0.0;
        const double d = (n_looks > 1.0)
                       ? (req->dwell_s - t_sap) / (n_looks - 1.0) : req->dwell_s;
        const double f_step = 1.0 / (2.0 * d);
        const double f_max  = (t_sap > 0.0) ? 1.0 / (2.0 * t_sap) : f_step;
        WORST(rs_v_add(out, &n, RS_VALIDATE_BAND,
              (f <= 0.0 || f < f_max) ? RS_V_PASS : RS_V_FAIL, "observable band",
              "alpha %.3f%% and overlap %.2f give %.0f sub-apertures. Each "
              "averages over %.4f s, so the band reaches %.3f Hz; the %.4f s "
              "step would suggest %.3f Hz, which overlap does not buy.%s",
              100.0 * req->alpha, req->overlap, n_looks, t_sap, f_max, d, f_step,
              (f > 0.0 && f >= f_max)
                ? " The target is ABOVE the band: past the sub-aperture's own "
                  "averaging response, where a reported peak cannot be signal."
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
        /* WARN rather than FAIL at both bounds: the mechanism is established
         * but the threshold is not, and a check that refuses a configuration
         * on an unmeasured boundary would be doing what this whole command
         * exists to prevent. */
        const rs_validate_level_t lvl = (eta <= RS_ETA_GOOD) ? RS_V_PASS
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
                : "The sub-look RESOLVES the target's own paired echoes, so what "
                  "the tracker follows is a train rather than a point. Expect "
                  "degradation; where it becomes fatal is NOT measured, and the "
                  "published validation operates at 0.39-0.69.",
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
              "uniform grid, so this enters as a time-base distortion -- a "
              "scale error on the frequency axis, not a phase artefact.%s",
              spread, gap_intervals,
              (lvl == RS_V_WARN) ? " Large enough to check before trusting a "
                                   "frequency." : ""));
    }

    /* A NOTE ON WHAT THIS CHECK IS AND IS NOT. Baehr (DGK Reihe C 719, KIT
     * 2013) Sect. 3.4.2 classes a biased PRF as a CLOCK error rather than a
     * timing error, and Sect. 3.4.3 states that "as long as coregistration is
     * implemented by amplitude cross-correlation, the interferometric phase
     * measurement is completely insensitive to errors in f_PRF and f_RSR".
     *
     * That does not make the check above redundant, because the mechanism here
     * is a different one. This project does not coregister two acquisitions; it
     * uses pulse times to place sub-aperture centres on the time axis of a
     * spectrum. A PRF error there distorts that axis directly, whatever the
     * phase does.
     *
     * What Baehr identifies as the error that DOES corrupt phase -- carrier
     * frequency drift -- this project has never examined. He cites short-term
     * ERS-1 drifts up to 82 Hz/s "that started and stopped abruptly and lasted
     * some tens of seconds", the same timescale as a long-dwell spotlight
     * collect, and notes such errors "were neither expected nor are they easy
     * to validate". Whether a drift within one dwell enters a single-pass
     * sub-aperture series, and what it would look like if it did, is an open
     * question and not checked anywhere here. */

    /* ---- the smallest motion the tracker could see ------------------------ */
    if (f > 0.0 && req->cell_m > 0.0 && req->upsample > 0) {
        const double proj = cos(req->incidence_rad);
        /* The floor is the excursion at which the tracker returns the target's
         * frequency rather than its own artefact, which is MEASURED and is not
         * the sub-pixel interpolation limit. Reporting the interpolation limit
         * here would overstate the sensitivity by the ratio printed below --
         * a factor of 57 at 1/40 px refinement -- and this project's whole
         * difficulty is that such a number looks like a demonstrated result. */
        const double floor_px = 0.5 * RS_TRACK_FLOOR_PX;
        const double interp_px = rs_validate_floor_px(req->upsample);
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
                     "The %.2f mm asked for is %.1fx that. %s Textured ground "
                     "lowers the floor to 4 px, so about 4/7 of the figure.",
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
              "the floor measured at 0.4 m cells is %.1f px p2p; carried to "
              "%.2f m cells -- which is NOT established, see rs_validate() -- "
              "that is a vertical amplitude of %.3f mm at %.3f Hz. %s Sub-pixel "
              "refinement at 1/%zu px reaches %.4f px, %.0fx finer -- that is "
              "an interpolation limit, not a demonstrated sensitivity.",
              RS_TRACK_FLOOR_PX, req->cell_m, 1000.0 * a_min, f, verdict,
              req->upsample, interp_px, floor_px / interp_px));
    }

    /* ---- is there any target amplitude this configuration can see? ------- */
    if (f > 0.0 && req->cell_m > 0.0 && t_sap > 0.0 && req->v_platform_ms > 0.0) {
        /* Two bounds squeeze the tracked azimuth excursion from both sides, and
         * a configuration is disqualified outright when they cross.
         *
         * CEILING. rs_microm_recommend_looks() requires the peak shift to stay
         * inside three quarters of a sub-look resolution cell; past that the
         * correlation argmax wraps. Measured directly: driving a target through
         * it made the tracker report 1.512 Hz for an injected 0.504 Hz -- the
         * third harmonic a wrapping sawtooth produces -- and at sixteen times
         * the amplitude the lowest bin.
         *
         * FLOOR. Below roughly 7 px peak-to-peak the tracker reports a fixed
         * spurious frequency instead of the target's. That artefact is not a
         * scene property: on a target with NO motion at all the uniform route
         * returns 1.569 Hz at prominence 27.9, and it returns the same 1.569 Hz
         * for every injected frequency from 0.2 to 1.4 Hz.
         *
         * The sub-look resolution is lambda*R/(2*v_p*t_sap), so a SHORT sub-
         * aperture raises the ceiling. Overlap shortens nothing -- it widens
         * each look's band, sharpening the sub-look and lowering the ceiling.
         * That is how a configuration ends up with no admissible amplitude at
         * all, which is what runs/giza/2026-07-30-validated-spot-khufu used.
         *
         * Both bounds are measured, on one geometry and one cell size, in
         * runs/.../POSITIVE-CONTROL.md. The crossing test below is the robust
         * part: it depends on their ratio, not on either absolute value. */
        const double res_sap  = req->lambda_m * req->slant_range_m
                              / (2.0 * req->v_platform_ms * t_sap);
        const double ceil_px  = 2.0 * 0.75 * res_sap / req->cell_m;  /* p2p */
        const double floor_px = RS_TRACK_FLOOR_PX;

        rs_validate_level_t lvl;
        char verdict[260];
        if (ceil_px <= floor_px) {
            lvl = RS_V_FAIL;
            snprintf(verdict, sizeof verdict,
                     "The ceiling is BELOW the floor, so no target amplitude "
                     "satisfies both: too small and the tracker reports its own "
                     "artefact, too large and the shift wraps. Shorten the "
                     "sub-aperture or reduce the overlap.");
        } else {
            const double proj = cos(req->incidence_rad);
            const double k = req->slant_range_m * 2.0 * M_PI * f * proj
                           / req->v_platform_ms;              /* m shift per m */
            const double lo = 0.5 * floor_px * req->cell_m / k;
            const double hi = 0.5 * ceil_px  * req->cell_m / k;
            lvl = RS_V_PASS;
            if (req->target_amp_m > 0.0)
                lvl = (req->target_amp_m >= lo && req->target_amp_m <= hi)
                    ? RS_V_PASS : RS_V_FAIL;
            snprintf(verdict, sizeof verdict,
                     "Admissible vertical amplitude at %.3f Hz is %.3f to "
                     "%.3f mm, a factor of %.1f. %s", f, 1000.0 * lo,
                     1000.0 * hi, ceil_px / floor_px,
                     (req->target_amp_m <= 0.0)
                       ? "No target amplitude given."
                       : (lvl == RS_V_PASS) ? "The target is inside it."
                       : (req->target_amp_m < lo)
                         ? "The target is BELOW it and would read as artefact."
                         : "The target is ABOVE it and the shift would wrap.");
        }
        WORST(rs_v_add(out, &n, RS_VALIDATE_AMBIGUITY, lvl, "ambiguity",
              "sub-look resolution %.2f m over %.4f s gives a wrap ceiling of "
              "%.1f px p2p against a %.1f px artefact floor. %s",
              res_sap, t_sap, ceil_px, floor_px, verdict));
    }

    /* ---- can the coherence mask reject anything at all? -------------------
     *
     * rs_splitband_shift() estimates coherence as the mean magnitude over every
     * pair of sub-looks. That is the standard estimator (ESA TM-19 Part C,
     * Eq. 1.14) applied to a stack that repeat-pass interferometry never has:
     * OVERLAPPING sub-apertures. Two looks whose bands overlap share spectral
     * content by construction, so their pairwise coherence is high whatever the
     * scene is doing -- for a white scene it is just the fraction of band they
     * share, about 1 - d*(1-overlap) for looks d apart.
     *
     * Averaged over all pairs this gives a FLOOR the estimator cannot report
     * below, and a mask set under that floor passes every window. TM-19's own
     * bias term (Eq. 1.15, sqrt(pi/4N) at true zero) is the same phenomenon
     * from sample count and is negligible here by comparison: at a 32x32 window
     * it is 0.03, where band overlap at 0.99 gives 0.57.
     *
     * This is not hypothetical. runs/giza/2026-07-30-uniform-phase-khufu
     * configuration A ran 128 looks at 0.99 overlap behind a 0.4 mask. Its
     * floor is 0.574. The mask admitted everything it was given. */
    if (req->alpha > 0.0 && req->overlap > 0.0 && req->overlap < 1.0) {
        const double dn = (1.0 / req->alpha) - req->overlap;
        const double nle = (dn > 0.0) ? dn / (1.0 - req->overlap) : 0.0;
        const size_t nl = (nle > 1.0) ? (size_t)(nle + 0.5) : 0;
        if (nl > 1) {
        double tot = 0.0, cnt = 0.0, ov_pairs = 0.0, all_pairs = 0.0;
        for (size_t d = 1; d < nl; d++) {
            const double frac = 1.0 - (double)d * (1.0 - req->overlap);
            const double wt = (double)(nl - d);
            tot += wt * ((frac > 0.0) ? frac : 0.0);
            cnt += wt;
            all_pairs += wt;
            if (frac > 0.0) ov_pairs += wt;
        }
        const double floor_g = (cnt > 0.0) ? tot / cnt : 0.0;
        const rs_validate_level_t lvl =
            (req->coherence_min <= 0.0)      ? RS_V_WARN
          : (req->coherence_min <= floor_g)  ? RS_V_FAIL
          : (req->coherence_min <= 1.5 * floor_g) ? RS_V_WARN : RS_V_PASS;
        WORST(rs_v_add(out, &n, RS_VALIDATE_COHERENCE_GATE, lvl, "coherence gate",
              "%.0f%% of look pairs share band, so an incoherent scene still "
              "reports |g| = %.3f. The mask is %.2f. %s",
              100.0 * (all_pairs > 0.0 ? ov_pairs / all_pairs : 0.0), floor_g,
              req->coherence_min,
              (req->coherence_min <= 0.0)
                ? "Masking is off, so nothing is rejected by design."
              : (req->coherence_min <= floor_g)
                ? "The mask is AT OR BELOW that floor and therefore vacuous: it "
                  "passes every window regardless of the scene."
              : (req->coherence_min <= 1.5 * floor_g)
                ? "Only just above the floor; it rejects little."
                : "Clear of the floor."));
        }
    }

    /* ---- the phase estimator's own noise floor ----------------------------
     *
     * The Cramer-Rao bound on interferometric phase (TM-19 Part C, Eq. 1.18;
     * Part A, Eq. 2.9) is sigma_phi = sqrt(1-g^2)/(g*sqrt(2N)), giving a line-
     * of-sight displacement noise of lambda*sigma_phi/(4*pi). This is the phase
     * route's equivalent of the excursion floor RS_VALIDATE_AMBIGUITY reports
     * for the correlation route, and it is roughly an order of magnitude
     * smaller -- which is the reason to prefer the phase estimator, and the
     * reason its null tests matter so much more.
     *
     * N is the number of INDEPENDENT samples in the window: cells only count
     * once per sub-look resolution cell, and a grid finer than the sub-look
     * merely oversamples. TM-19 notes the bound is reasonable only for N > 4,
     * and that at N of 1 or 2 the phase dispersion is unusable at any SNR.
     *
     * THIS CHECK IS ONLY AS GOOD AS THE COHERENCE FED TO IT. The bound falls
     * as coherence rises, so if RS_VALIDATE_COHERENCE_GATE reports that band
     * overlap inflates the estimator, the floor printed here is optimistic by
     * the same amount. Read the two together; they are not independent. */
    if (req->win > 0 && req->cell_m > 0.0 && t_sap > 0.0
        && req->coherence_min > 0.0 && req->lambda_m > 0.0) {
        const double res_sap = req->lambda_m * req->slant_range_m
                             / (2.0 * req->v_platform_ms * t_sap);
        const double per_edge = (double)req->win * req->cell_m / res_sap;
        const double n_ind = (per_edge < 1.0 ? 1.0 : per_edge) * (double)req->win;
        const double g = req->coherence_min;
        const double sphi = sqrt(1.0 - g * g) / (g * sqrt(2.0 * n_ind));
        const double d_m = req->lambda_m * sphi / (4.0 * M_PI);
        const rs_validate_level_t lvl = (n_ind <= 4.0) ? RS_V_WARN : RS_V_PASS;
        WORST(rs_v_add(out, &n, RS_VALIDATE_PHASE_FLOOR, lvl, "phase floor",
              "at the %.2f mask over %.0f independent samples the CRLB is "
              "%.3f rad, a line-of-sight noise of %.4f mm per look.%s",
              g, n_ind, sphi, 1000.0 * d_m,
              (n_ind <= 4.0)
                ? " Under 5 independent samples the bound does not hold and the"
                  " phase dispersion is unusable at any SNR."
                : ""));
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
