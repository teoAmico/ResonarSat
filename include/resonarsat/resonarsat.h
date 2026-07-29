/* ResonarSat -- common status codes and error reporting.
 *
 * MIT licensed, see LICENSE at the repository root. */

#ifndef RESONARSAT_H
#define RESONARSAT_H

#include <stddef.h>

/* Status codes returned by every fallible function in the library.
 *
 * The contract, which readers in particular must honour: a malformed, truncated
 * or hostile input file produces one of these codes plus a descriptive message
 * retrievable with rs_last_error(). It never produces a crash, a partial write
 * into a caller's buffer, or a silently zero-filled result. */
typedef enum {
    RS_OK = 0,
    RS_ERR_ARG,           /* caller passed a NULL or out-of-range argument */
    RS_ERR_ALLOC,         /* out of memory */
    RS_ERR_IO,            /* open/read/write failed at the OS level */
    RS_ERR_FORMAT,        /* file is not the format it claims, or is truncated */
    RS_ERR_MISSING_META,  /* required metadata field absent (PRF, orbit, ...) */
    RS_ERR_UNSUPPORTED,   /* valid input this build cannot handle yet */
    RS_ERR_SINGULAR,      /* numerical solve failed to converge */
    RS_ERR_RANGE          /* computed quantity outside a physically sane range */
} resonarsat_status_t;

/* Return a short, stable, human-readable name for a status code, such as
 * "format error". The returned pointer is a string literal owned by the
 * library and must not be freed. Unknown codes yield "unknown status" rather
 * than NULL, so callers may always pass the result straight to printf. */
const char *rs_status_str(resonarsat_status_t st);

/* Record a descriptive message for the failure the caller is about to return.
 *
 * Takes a printf-style format string. The message is stored in a fixed-size
 * per-thread buffer, so it is safe to call from OpenMP parallel regions and
 * each thread sees its own most recent error. Messages longer than the buffer
 * are truncated rather than overflowing.
 *
 * The intended pattern is to call this immediately before returning a non-OK
 * status, at the point where the specific detail (file offset, field name,
 * observed value) is still in scope:
 *
 *     rs_set_error("annotation: azimuth spacing %g m outside [0.1, 1000]", v);
 *     return RS_ERR_RANGE;
 */
void rs_set_error(const char *fmt, ...);

/* Return the calling thread's most recent rs_set_error() message, or the empty
 * string if that thread has not recorded one. The pointer refers to thread-local
 * storage that the next rs_set_error() on the same thread overwrites, so copy
 * the string if it must outlive the immediate error path. Never returns NULL. */
const char *rs_last_error(void);

/* Print "context: <status name> (<detail>)" to stderr, where the detail is the
 * current rs_last_error() message when one is set. This is the standard way for
 * the CLI to surface a failure; library code should set the error and return a
 * status rather than printing anything itself. */
void rs_report_error(const char *context, resonarsat_status_t st);

#endif /* RESONARSAT_H */
