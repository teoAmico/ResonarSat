/* Raster output, checked by decoding what was written.
 *
 * The PNG encoder is hand-rolled to keep the project free of external
 * libraries, which means nothing else validates its output: an encoder that
 * emits a plausible header and a broken stream produces a file that 'file'
 * identifies happily and no viewer opens. So these cases decode the bytes back
 * -- chunk lengths, CRC-32 over type and payload, the zlib wrapper, the stored
 * DEFLATE blocks and the Adler-32 -- rather than checking that a file appeared.
 *
 * The strongest case is the last one: the greyscale PNG's pixel payload must be
 * byte-identical to the PGM written from the same data. Both come from the same
 * quantisation, so any disagreement is a fault in the container layer, and that
 * comparison needs no reference file to be kept in the repository. */

#include "resonarsat/raster.h"
#include "rs_test.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Read a whole file. Returns NULL and leaves '*len' zero when it cannot. */
static unsigned char *rs_slurp(const char *path, size_t *len)
{
    *len = 0;
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    const long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    rewind(f);
    unsigned char *b = malloc((size_t)n);
    if (!b) { fclose(f); return NULL; }
    const size_t got = fread(b, 1, (size_t)n, f);
    fclose(f);
    if (got != (size_t)n) { free(b); return NULL; }
    *len = got;
    return b;
}

/* CRC-32, independently written from the encoder's so a shared bug in the
 * polynomial or the initial value cannot cancel out. */
static unsigned long rs_ref_crc32(const unsigned char *b, size_t n)
{
    unsigned long c = 0xFFFFFFFFUL;
    for (size_t i = 0; i < n; i++) {
        c ^= b[i];
        for (int k = 0; k < 8; k++) {
            c = (c & 1) ? (0xEDB88320UL ^ (c >> 1)) : (c >> 1);
        }
    }
    return c ^ 0xFFFFFFFFUL;
}

static unsigned long rs_be32(const unsigned char *p)
{
    return ((unsigned long)p[0] << 24) | ((unsigned long)p[1] << 16) |
           ((unsigned long)p[2] << 8)  |  (unsigned long)p[3];
}

/* Walk the chunks, verify every CRC, and inflate the stored-block stream.
 * Returns the raw filtered bytes, with the geometry in the out parameters. */
static unsigned char *rs_png_decode(const unsigned char *d, size_t n,
                                    size_t *w, size_t *h, int *channels,
                                    size_t *raw_len)
{
    static const unsigned char sig[8] = { 137, 'P', 'N', 'G', '\r', '\n', 26, '\n' };
    if (n < 8 || memcmp(d, sig, 8) != 0) return NULL;

    size_t o = 8;
    unsigned char *idat = NULL;
    size_t idat_len = 0;
    int seen_end = 0;
    *w = *h = 0; *channels = 0;

    while (o + 12 <= n) {
        const size_t len = (size_t)rs_be32(d + o);
        if (o + 12 + len > n) { free(idat); return NULL; }
        const unsigned char *type = d + o + 4;
        const unsigned char *pay  = d + o + 8;

        /* CRC covers the type and the payload, not the length field. */
        unsigned char *tmp = malloc(4 + len);
        if (!tmp) { free(idat); return NULL; }
        memcpy(tmp, type, 4);
        memcpy(tmp + 4, pay, len);
        const unsigned long want = rs_be32(d + o + 8 + len);
        const unsigned long got  = rs_ref_crc32(tmp, 4 + len);
        free(tmp);
        if (want != got) { free(idat); return NULL; }

        if (memcmp(type, "IHDR", 4) == 0 && len == 13) {
            *w = (size_t)rs_be32(pay);
            *h = (size_t)rs_be32(pay + 4);
            if (pay[8] != 8) { free(idat); return NULL; }      /* bit depth */
            if (pay[9] == 0) *channels = 1;
            else if (pay[9] == 2) *channels = 3;
            else { free(idat); return NULL; }
            if (pay[10] || pay[11] || pay[12]) { free(idat); return NULL; }
        } else if (memcmp(type, "IDAT", 4) == 0) {
            unsigned char *g = realloc(idat, idat_len + len);
            if (!g) { free(idat); return NULL; }
            idat = g;
            memcpy(idat + idat_len, pay, len);
            idat_len += len;
        } else if (memcmp(type, "IEND", 4) == 0) {
            seen_end = 1;
        }
        o += 12 + len;
    }
    if (!seen_end || !idat || idat_len < 6) { free(idat); return NULL; }

    /* zlib header: the two bytes must form a multiple of 31, and the method
     * nibble must say deflate. */
    if ((idat[0] & 0x0F) != 8) { free(idat); return NULL; }
    if ((((unsigned)idat[0] << 8) | idat[1]) % 31u != 0) { free(idat); return NULL; }

    /* Inflate, accepting stored blocks only -- which is all the encoder emits,
     * so anything else here is a fault rather than an unsupported feature. */
    unsigned char *raw = NULL;
    size_t rn = 0;
    size_t p = 2;
    int final = 0;
    while (!final) {
        if (p + 5 > idat_len) { free(raw); free(idat); return NULL; }
        const unsigned char bh = idat[p];
        final = bh & 1;
        if ((bh >> 1) & 3) { free(raw); free(idat); return NULL; }   /* BTYPE != 00 */
        const size_t blen  = (size_t)idat[p + 1] | ((size_t)idat[p + 2] << 8);
        const size_t nlen  = (size_t)idat[p + 3] | ((size_t)idat[p + 4] << 8);
        if ((blen ^ 0xFFFFu) != nlen) { free(raw); free(idat); return NULL; }
        p += 5;
        if (p + blen > idat_len) { free(raw); free(idat); return NULL; }
        unsigned char *g = realloc(raw, rn + blen);
        if (!g) { free(raw); free(idat); return NULL; }
        raw = g;
        memcpy(raw + rn, idat + p, blen);
        rn += blen;
        p += blen;
    }

    /* Adler-32 of the uncompressed data closes the stream. */
    if (p + 4 > idat_len) { free(raw); free(idat); return NULL; }
    unsigned long a = 1, b = 0;
    for (size_t i = 0; i < rn; i++) { a = (a + raw[i]) % 65521UL; b = (b + a) % 65521UL; }
    if (rs_be32(idat + p) != ((b << 16) | a)) { free(raw); free(idat); return NULL; }
    free(idat);

    *raw_len = rn;
    return raw;
}

int main(void)
{
    const size_t R = 37, C = 53;                /* deliberately not round, and
                                                 * not a multiple of anything */
    double *map = malloc(R * C * sizeof *map);
    if (!map) return 1;
    for (size_t r = 0; r < R; r++) {
        for (size_t c = 0; c < C; c++) {
            map[r * C + c] = sin((double)c / 7.0) * cos((double)r / 5.0);
        }
    }

    RS_CASE("a greyscale PNG decodes, with every checksum intact");
    {
        RS_CHECK_OK(rs_raster_write_map(map, R, C, "t_raster_gray.png",
                                        0.0, 0.0, RS_PALETTE_GRAY));
        size_t n = 0;
        unsigned char *d = rs_slurp("t_raster_gray.png", &n);
        RS_CHECK(d != NULL);
        size_t w = 0, h = 0, rn = 0; int ch = 0;
        unsigned char *raw = rs_png_decode(d, n, &w, &h, &ch, &rn);
        RS_CHECK(raw != NULL);
        printf("    %zux%zu, %d channel(s), %zu raw bytes\n", w, h, ch, rn);
        RS_CHECK(w == C && h == R && ch == 1);
        RS_CHECK(rn == (C * 1 + 1) * R);
        free(raw); free(d);
    }

    RS_CASE("a viridis PNG decodes as three-channel colour");
    {
        RS_CHECK_OK(rs_raster_write_map(map, R, C, "t_raster_col.png",
                                        0.0, 0.0, RS_PALETTE_VIRIDIS));
        size_t n = 0;
        unsigned char *d = rs_slurp("t_raster_col.png", &n);
        RS_CHECK(d != NULL);
        size_t w = 0, h = 0, rn = 0; int ch = 0;
        unsigned char *raw = rs_png_decode(d, n, &w, &h, &ch, &rn);
        RS_CHECK(raw != NULL);
        RS_CHECK(w == C && h == R && ch == 3);
        RS_CHECK(rn == (C * 3 + 1) * R);

        /* The ramp must actually vary; a constant image would decode perfectly
         * and mean the palette never ran. */
        int varied = 0;
        for (size_t i = 1; i < rn && !varied; i++) if (raw[i] != raw[1]) varied = 1;
        RS_CHECK(varied);
        free(raw); free(d);
    }

    RS_CASE("the greyscale PNG carries the same pixels as the PGM");
    {
        RS_CHECK_OK(rs_raster_write_map(map, R, C, "t_raster.pgm",
                                        0.0, 0.0, RS_PALETTE_GRAY));
        size_t pn = 0;
        unsigned char *pgm = rs_slurp("t_raster.pgm", &pn);
        RS_CHECK(pgm != NULL);

        /* Skip the three PGM header lines. */
        size_t off = 0; int nl = 0;
        while (off < pn && nl < 3) if (pgm[off++] == '\n') nl++;
        RS_CHECK(pn - off == R * C);

        size_t n = 0;
        unsigned char *d = rs_slurp("t_raster_gray.png", &n);
        size_t w = 0, h = 0, rn = 0; int ch = 0;
        unsigned char *raw = rs_png_decode(d, n, &w, &h, &ch, &rn);
        RS_CHECK(raw != NULL);

        int same = 1;
        for (size_t r = 0; r < R && same; r++) {
            /* Every row must open with filter type 0, then match the PGM. */
            if (raw[r * (C + 1)] != 0) same = 0;
            else if (memcmp(raw + r * (C + 1) + 1, pgm + off + r * C, C) != 0) same = 0;
        }
        printf("    %zu bytes compared\n", R * C);
        RS_CHECK(same);
        free(raw); free(d); free(pgm);
    }

    RS_CASE("a palette request on a PGM path degrades instead of failing");
    {
        /* Asking for colour in a container that cannot carry it is cosmetic,
         * and refusing to write the image would be a worse answer than writing
         * it in grey. */
        RS_CHECK_OK(rs_raster_write_map(map, R, C, "t_raster_deg.pgm",
                                        0.0, 0.0, RS_PALETTE_VIRIDIS));
        size_t n = 0;
        unsigned char *d = rs_slurp("t_raster_deg.pgm", &n);
        RS_CHECK(d != NULL);
        RS_CHECK(n > 3 && d[0] == 'P' && d[1] == '5');
        free(d);
    }

    RS_CASE("a map that is entirely NaN still produces a readable image");
    {
        /* Autoscaling over NaN would previously give NaN limits and a blank
         * frame; a masked-out tomogram is a real result and must still render. */
        double *bad = malloc(R * C * sizeof *bad);
        RS_CHECK(bad != NULL);
        for (size_t i = 0; i < R * C; i++) bad[i] = NAN;
        RS_CHECK_OK(rs_raster_write_map(bad, R, C, "t_raster_nan.png",
                                        0.0, 0.0, RS_PALETTE_VIRIDIS));
        size_t n = 0;
        unsigned char *d = rs_slurp("t_raster_nan.png", &n);
        RS_CHECK(d != NULL);
        size_t w = 0, h = 0, rn = 0; int ch = 0;
        unsigned char *raw = rs_png_decode(d, n, &w, &h, &ch, &rn);
        RS_CHECK(raw != NULL);
        RS_CHECK(w == C && h == R);
        free(raw); free(d); free(bad);
    }

    RS_CASE("bad arguments are refused rather than written");
    {
        RS_CHECK_ERR(rs_raster_write_map(NULL, R, C, "t_x.png", 0, 0, RS_PALETTE_GRAY),
                     RS_ERR_ARG);
        RS_CHECK_ERR(rs_raster_write_map(map, 0, C, "t_x.png", 0, 0, RS_PALETTE_GRAY),
                     RS_ERR_ARG);
        RS_CHECK_ERR(rs_raster_write_map(map, R, 0, "t_x.png", 0, 0, RS_PALETTE_GRAY),
                     RS_ERR_ARG);
        RS_CHECK_ERR(rs_raster_write_map(map, R, C, NULL, 0, 0, RS_PALETTE_GRAY),
                     RS_ERR_ARG);
        /* An unwritable path must report IO, not crash or silently succeed. */
        RS_CHECK(rs_raster_write_map(map, R, C, "/nonexistent-dir/x.png",
                                     0, 0, RS_PALETTE_GRAY) != RS_OK);
    }

    free(map);
    remove("t_raster_gray.png"); remove("t_raster_col.png");
    remove("t_raster.pgm"); remove("t_raster_deg.pgm"); remove("t_raster_nan.png");

    RS_TEST_END();
}
