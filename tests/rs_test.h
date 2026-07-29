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
#define RS_TEST_END() do { \
    if (rs_test_failures) { \
        printf("%d failure(s)\n", rs_test_failures); \
        return EXIT_FAILURE; \
    } \
    printf("ok\n"); \
    return EXIT_SUCCESS; \
} while (0)

#endif /* RS_TEST_H */
