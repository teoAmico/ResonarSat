/* Malformed-input corpus for the remaining parsers.
 *
 * The CPHD 1.x reader got its corpus when it was written (test_cphd.c). These
 * are the two that did not: the UAVSAR reader, which parses an untrusted ASCII
 * annotation and a flat binary of unknown provenance, and rs_cphd_read(), which
 * parses this project's own interchange format -- "our own format" being no
 * protection at all, since the file on disk arrives from wherever the user got
 * it.
 *
 * Each case damages one thing in an otherwise valid pair and asserts the status.
 * The status code matters less than the two invariants behind it: nothing
 * crashes, and nothing is left allocated. A parser that returns an error while
 * leaving a half-filled struct behind is worse than one that crashes, because
 * the caller proceeds on data that looks real.
 *
 * The plan lists this as mandatory before the repository is published. These are
 * the inputs a public repository receives. */

#include "resonarsat/readers.h"
#include "resonarsat/focus.h"
#include "rs_test.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define ROWS 8u
#define COLS 6u

/* Write a UAVSAR annotation and its companion binary.
 *
 * 'damage' names one deliberate defect, or NULL for a well-formed pair. Returns
 * 0 on success. Values are plausible rather than realistic; the parser's job is
 * to reject what is malformed, not to judge what is likely. */
static int write_uavsar(const char *ann_path, const char *slc_path, const char *damage)
{
    FILE *a = fopen(ann_path, "w");
    if (!a) return -1;
    fprintf(a, "; a UAVSAR-style annotation\n");
    if (!(damage && !strcmp(damage, "no_rows"))) {
        fprintf(a, "Number of SLC Lines = %s\n",
                (damage && !strcmp(damage, "rows_nan"))     ? "not-a-number" :
                (damage && !strcmp(damage, "rows_zero"))    ? "0" :
                (damage && !strcmp(damage, "rows_huge"))    ? "1e12" : "8");
    }
    fprintf(a, "Number of SLC Range Bins = %u\n", COLS);
    fprintf(a, "Azimuth Spacing = 1.0\n");
    fprintf(a, "Range Spacing = 1.0\n");
    fclose(a);

    if (damage && !strcmp(damage, "no_slc")) {
        remove(slc_path);
        return 0;
    }

    FILE *f = fopen(slc_path, "wb");
    if (!f) return -1;
    size_t n = (size_t)ROWS * COLS;
    if (damage && !strcmp(damage, "short_slc")) n /= 2;      /* declared bigger than it is */
    for (size_t i = 0; i < n; i++) {
        const float iq[2] = { (float)i, (float)-i };
        fwrite(iq, sizeof iq, 1, f);
    }
    fclose(f);
    return 0;
}

/* Write a phase-history file in this project's interchange format, optionally
 * damaged. Mirrors rs_cphd_write()'s layout; see rs_cphd_read(). */
static int write_interchange(const char *path, const char *damage)
{
    struct {
        uint32_t magic, version;
        uint64_t n_pulse, n_rbin;
        double   r_near, dr, fc, lambda, prf;
    } h;
    memset(&h, 0, sizeof h);
    h.magic   = (damage && !strcmp(damage, "magic")) ? 0xDEADBEEFu : 0x52534348u;
    h.version = (damage && !strcmp(damage, "version")) ? 99u : 2u;
    h.n_pulse = (damage && !strcmp(damage, "dims_huge")) ? 1ull << 40 : 16;
    h.n_rbin  = (damage && !strcmp(damage, "dims_zero")) ? 0 : 8;
    h.r_near  = 500000.0;
    h.dr      = 1.0;
    h.fc      = 9.6e9;
    h.lambda  = 0.031;
    h.prf     = 1000.0;

    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    if (damage && !strcmp(damage, "stub")) {
        /* Too short to contain a header at all. */
        fwrite(&h, 4, 1, f);
        fclose(f);
        return 0;
    }

    fwrite(&h, sizeof h, 1, f);

    /* The payload rs_cphd_read() expects: positions, times, reference ranges,
     * then the samples. Truncated deliberately for the "truncated" case. */
    const size_t np = 16, nr = 8;
    const size_t n_sample = (damage && !strcmp(damage, "truncated")) ? np * nr / 4 : np * nr;
    for (size_t i = 0; i < np * 3; i++) { const double z = 0.0; fwrite(&z, sizeof z, 1, f); }
    for (size_t i = 0; i < np; i++)     { const double t = (double)i * 1e-3; fwrite(&t, sizeof t, 1, f); }
    for (size_t i = 0; i < np; i++)     { const double r = 500000.0; fwrite(&r, sizeof r, 1, f); }
    for (size_t i = 0; i < n_sample; i++) {
        const float iq[2] = { 1.0f, 0.0f };
        fwrite(iq, sizeof iq, 1, f);
    }
    fclose(f);
    return 0;
}

int main(void)
{
    const char *dir = getenv("TMPDIR");
    char ann[512], slc[512], cph[512];
    const int pid = (int)getpid();
    snprintf(ann, sizeof ann, "%srs_ann_%d.txt",  dir ? dir : "/tmp/", pid);
    snprintf(slc, sizeof slc, "%srs_slc_%d.bin",  dir ? dir : "/tmp/", pid);
    snprintf(cph, sizeof cph, "%srs_int_%d.cphd", dir ? dir : "/tmp/", pid);

    RS_CASE("a well-formed UAVSAR pair reads");
    {
        RS_CHECK(write_uavsar(ann, slc, NULL) == 0);
        rs_slc_t img;
        RS_CHECK_OK(rs_read_uavsar(slc, ann, &img));
        RS_CHECK(img.n_az == ROWS && img.n_rg == COLS);
        rs_slc_free(&img);
    }

    RS_CASE("malformed UAVSAR input is refused without crashing or leaking");
    {
        const struct { const char *damage; const char *why; } cases[] = {
            { "no_rows",   "annotation omits the line count" },
            { "rows_nan",  "line count is not a number" },
            { "rows_zero", "zero rows" },
            { "rows_huge", "implausible dimensions" },
            { "short_slc", "binary smaller than the declared dimensions" },
            { "no_slc",    "binary missing entirely" },
        };
        for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
            RS_CHECK(write_uavsar(ann, slc, cases[i].damage) == 0);
            rs_slc_t img;
            memset(&img, 0xA5, sizeof img);      /* poison: a partial write shows up */
            const resonarsat_status_t st = rs_read_uavsar(slc, ann, &img);
            printf("    %-10s -> %-22s (%s)\n", cases[i].damage,
                   rs_status_str(st), cases[i].why);
            RS_CHECK(st != RS_OK);
            /* The contract is that a failed read leaves no allocation, so this
             * must be safe on a struct the reader never completed. */
            rs_slc_free(&img);
        }
    }

    RS_CASE("null and missing paths are refused");
    {
        rs_slc_t img;
        RS_CHECK_ERR(rs_read_uavsar(NULL, ann, &img), RS_ERR_ARG);
        RS_CHECK_ERR(rs_read_uavsar(slc, NULL, &img), RS_ERR_ARG);
        RS_CHECK_ERR(rs_read_uavsar(slc, ann, NULL), RS_ERR_ARG);
        RS_CHECK_ERR(rs_read_uavsar("/nonexistent/x.bin", "/nonexistent/x.ann", &img),
                     RS_ERR_IO);
    }

    RS_CASE("a well-formed interchange file reads");
    {
        RS_CHECK(write_interchange(cph, NULL) == 0);
        rs_cphd_t c;
        RS_CHECK_OK(rs_cphd_read(&c, cph));
        RS_CHECK(c.n_pulse == 16 && c.n_rbin == 8);
        rs_cphd_free(&c);
    }

    RS_CASE("malformed interchange input is refused");
    {
        const struct { const char *damage; const char *why; } cases[] = {
            { "stub",      "shorter than its own header" },
            { "magic",     "wrong magic number" },
            { "version",   "future version" },
            { "dims_huge", "dimensions that would overflow an allocation" },
            { "dims_zero", "zero range bins" },
            { "truncated", "payload cut short" },
        };
        for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
            RS_CHECK(write_interchange(cph, cases[i].damage) == 0);
            rs_cphd_t c;
            memset(&c, 0x5A, sizeof c);
            const resonarsat_status_t st = rs_cphd_read(&c, cph);
            printf("    %-10s -> %-22s (%s)\n", cases[i].damage,
                   rs_status_str(st), cases[i].why);
            RS_CHECK(st != RS_OK);
            rs_cphd_free(&c);
        }
    }

    /* SICD damage cases. A real product is not synthesised here -- a valid NITF
     * needs several hundred correlated fields -- so these exercise the container
     * checks that fire before any of that matters: the magic, the declared
     * lengths, and the arithmetic tying the field walk to the header length. */
    RS_CASE("malformed SICD input is refused");
    {
        char sic[512];
        snprintf(sic, sizeof sic, "%srs_sicd_%d.nitf", dir ? dir : "/tmp/", pid);

        const struct { const char *bytes; size_t n; const char *why; } cases[] = {
            { "not a nitf file at all, not even close........................", 60,
              "wrong magic" },
            { "NITF02.10", 9, "correct magic but nothing after it" },
            { "NITF02.10" "0123456789", 19, "truncated inside the header" },
        };
        for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
            FILE *g = fopen(sic, "wb");
            RS_CHECK(g != NULL);
            fwrite(cases[i].bytes, 1, cases[i].n, g);
            fclose(g);

            rs_slc_t img;
            memset(&img, 0x3C, sizeof img);
            const resonarsat_status_t st = rs_read_sicd(sic, &img);
            printf("    %-24s -> %s\n", cases[i].why, rs_status_str(st));
            RS_CHECK(st != RS_OK);
            rs_slc_free(&img);

            /* The metadata-only path must refuse exactly what the full read
             * refuses. It skips the pixel pass, not the container checks, and a
             * screening command that accepts files the processing path rejects
             * would pass collects that cannot then be run. */
            rs_slc_t meta;
            memset(&meta, 0x3C, sizeof meta);
            const resonarsat_status_t mst = rs_read_sicd_meta(sic, &meta);
            RS_CHECK(mst == st);
            rs_slc_free(&meta);
        }
        RS_CHECK_ERR(rs_read_sicd(NULL, &(rs_slc_t){0}), RS_ERR_ARG);
        RS_CHECK_ERR(rs_read_sicd(sic, NULL), RS_ERR_ARG);
        RS_CHECK_ERR(rs_read_sicd_meta(NULL, &(rs_slc_t){0}), RS_ERR_ARG);
        RS_CHECK_ERR(rs_read_sicd_meta(sic, NULL), RS_ERR_ARG);
        rs_slc_t img2;
        RS_CHECK_ERR(rs_read_sicd("/nonexistent/x.nitf", &img2), RS_ERR_IO);
        rs_slc_t img3;
        RS_CHECK_ERR(rs_read_sicd_meta("/nonexistent/x.nitf", &img3), RS_ERR_IO);
        remove(sic);
    }

    RS_CASE("a directory and a null path are refused");
    {
        rs_cphd_t c;
        RS_CHECK_ERR(rs_cphd_read(&c, NULL), RS_ERR_ARG);
        RS_CHECK_ERR(rs_cphd_read(NULL, cph), RS_ERR_ARG);
        RS_CHECK(rs_cphd_read(&c, dir ? dir : "/tmp") != RS_OK);
    }

    remove(ann);
    remove(slc);
    remove(cph);
    RS_TEST_END();
}
