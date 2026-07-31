/* Minimal assertion harness for the ctest suite.
 *
 * Deliberately tiny and dependency-free: the point of these tests is to pin
 * numerical relations and algorithm behaviour, and a framework would add build
 * surface without adding confidence. Each test binary runs a handful of CHECK
 * macros and exits non-zero on the first failure, printing the file, line,
 * expected and observed values. */

#ifndef RS_TEST_H
#define RS_TEST_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int rs_test_failures = 0;
static const char *rs_test_current = "";

/* Announce the test currently running. The name is printed with each failure,
 * so a ctest log identifies the case without needing line numbers. */
#define RS_CASE(name) do { rs_test_current = (name); printf("  case: %s\n", (name)); } while (0)

/* Assert that a condition holds. */
#define RS_CHECK(cond) do { \
    if (!(cond)) { \
        printf("FAIL [%s] %s:%d: %s\n", rs_test_current, __FILE__, __LINE__, #cond); \
        rs_test_failures++; \
    } \
} while (0)

/* Assert that 'got' is within 'tol' of 'want', absolute tolerance. */
#define RS_CHECK_NEAR(got, want, tol) do { \
    const double _g = (double)(got), _w = (double)(want), _t = (double)(tol); \
    if (!(fabs(_g - _w) <= _t)) { \
        printf("FAIL [%s] %s:%d: %s = %.9g, expected %.9g +/- %.3g (off by %.3g)\n", \
               rs_test_current, __FILE__, __LINE__, #got, _g, _w, _t, fabs(_g - _w)); \
        rs_test_failures++; \
    } \
} while (0)

/* Assert that 'got' is within 'frac' relative tolerance of 'want'. */
#define RS_CHECK_REL(got, want, frac) do { \
    const double _g = (double)(got), _w = (double)(want), _f = (double)(frac); \
    const double _tol = fabs(_w) * _f; \
    if (!(fabs(_g - _w) <= _tol)) { \
        printf("FAIL [%s] %s:%d: %s = %.9g, expected %.9g +/- %.2f%% (off by %.3g)\n", \
               rs_test_current, __FILE__, __LINE__, #got, _g, _w, _f * 100.0, fabs(_g - _w)); \
        rs_test_failures++; \
    } \
} while (0)

/* Assert a call returned RS_OK, printing the status name and detail if not. */
#define RS_CHECK_OK(expr) do { \
    resonarsat_status_t _st = (expr); \
    if (_st != RS_OK) { \
        printf("FAIL [%s] %s:%d: %s -> %s (%s)\n", rs_test_current, __FILE__, __LINE__, \
               #expr, rs_status_str(_st), rs_last_error()); \
        rs_test_failures++; \
    } \
} while (0)

/* Assert a call returned a specific non-OK status. Used to pin the error
 * contract: malformed input must fail with a described, specific code. */
#define RS_CHECK_ERR(expr, want) do { \
    resonarsat_status_t _st = (expr); \
    if (_st != (want)) { \
        printf("FAIL [%s] %s:%d: %s -> %s, expected %s\n", rs_test_current, \
               __FILE__, __LINE__, #expr, rs_status_str(_st), rs_status_str(want)); \
        rs_test_failures++; \
    } \
} while (0)

/* Report the tally and produce the process exit code. */
/* Does a reported frequency TRACK an injected one, or is it fixed?
 *
 * WHY THIS EXISTS. The obvious criterion -- |reported - injected| < 2 bins at
 * one frequency -- is far too loose, and this project has drawn four wrong
 * conclusions from it in a single day. A tracker that emits a FIXED spurious
 * frequency passes it whenever that fixed value happens to fall near the
 * injected one, which at coarse bin spacing is often. Every one of those wrong
 * conclusions came from a single-point match; every one would have been caught
 * by sweeping the injected frequency and asking whether the answer moved.
 * See runs/giza/2026-07-30-validated-spot-khufu/POSITIVE-CONTROL.md.
 *
 * Two numbers are returned and both must be checked:
 *
 *   slope   least squares of reported against injected. A working chain gives
 *           1; a fixed artefact gives 0. This is the part a single point
 *           cannot test at all.
 *   rms     root mean square of (reported - injected), in Hz. Tightness. Half
 *           a bin is a reasonable bound and is what the working operating
 *           point achieves with room to spare.
 *
 * Neither alone is sufficient. Slope 1 with a large offset still tracks but is
 * biased; a small rms with slope 0 means the sweep was too narrow to tell.
 *
 * AND NEITHER IS SUFFICIENT ON ONE SCENE REALISATION. This criterion sweeps
 * frequency; it does not sweep the speckle. Measured on 2026-07-31: at 1.1 Hz
 * with everything else held fixed, zero overlap missed on one clutter seed and
 * recovered on two others, while 0.5 overlap recovered on two and reported
 * 6.275 Hz on the third. A five-frequency sweep on a single seed gave slope
 * 1.004 and rms 0.0030 Hz -- a clean pass -- for a configuration that fails
 * outright on a different realisation of the same scene.
 *
 * So a passing slope and rms establish that a chain tracked ON THAT SCENE.
 * Concluding anything about a CONFIGURATION needs the sweep repeated over
 * independent realisations (sim_cphd --seed) and the verdicts pooled. Single
 * realisations are how this project drew five wrong conclusions in one day.
 * Returns 0 and leaves the outputs at 0 when n < 3, since two points always
 * fit a line exactly and can never fail this. */
static inline int rs_track_fit(const double *injected, const double *reported,
                        size_t n, double *slope_out, double *rms_out)
{
    if (slope_out) *slope_out = 0.0;
    if (rms_out)   *rms_out = 0.0;
    if (!injected || !reported || n < 3) return 0;

    double mx = 0.0, my = 0.0;
    for (size_t i = 0; i < n; i++) { mx += injected[i]; my += reported[i]; }
    mx /= (double)n; my /= (double)n;

    double sxy = 0.0, sxx = 0.0, sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        sxy += (injected[i] - mx) * (reported[i] - my);
        sxx += (injected[i] - mx) * (injected[i] - mx);
        const double r = reported[i] - injected[i];
        sq += r * r;
    }
    if (slope_out) *slope_out = (sxx > 0.0) ? sxy / sxx : 0.0;
    if (rms_out)   *rms_out = sqrt(sq / (double)n);
    return 1;
}

#define RS_TEST_END() do { \
    if (rs_test_failures) { \
        printf("%d failure(s)\n", rs_test_failures); \
        return EXIT_FAILURE; \
    } \
    printf("ok\n"); \
    return EXIT_SUCCESS; \
} while (0)

#endif /* RS_TEST_H */
