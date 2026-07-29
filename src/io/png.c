/* PNG encoding with stored DEFLATE blocks. See png.h for why. */

#include "png.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 * Checksums.
 *
 * PNG needs two, and they are different algorithms over different data: CRC-32
 * covers each chunk's type and payload, while Adler-32 covers the UNCOMPRESSED
 * bytes of the zlib stream. Mixing them up produces a file that looks correct
 * and fails to decode, so they are kept apart deliberately.
 * ------------------------------------------------------------------------- */

/* Lookup table for CRC-32, and whether it has been filled in yet. */
static unsigned long rs_crc_table[256];
static int rs_crc_table_ready = 0;

/* Fill the CRC-32 table, once, on first use.
 *
 * PNG specifies the reflected polynomial 0xEDB88320. Computing the 256 entries
 * here rather than pasting a 1 KB literal keeps the constant out of the source,
 * where a single mistyped digit would produce files that fail to decode in a way
 * no amount of reading the encoder would explain. */
static void rs_crc_build_table(void)
{
    for (unsigned long n = 0; n < 256; n++) {
        unsigned long c = n;
        for (int k = 0; k < 8; k++) {
            c = (c & 1) ? (0xEDB88320UL ^ (c >> 1)) : (c >> 1);
        }
        rs_crc_table[n] = c;
    }
    rs_crc_table_ready = 1;
}

/* Fold 'len' bytes into a running CRC-32 and return the new running value.
 *
 * Incremental rather than whole-buffer because a PNG chunk's checksum spans its
 * four-character type AND its payload, which are never contiguous in memory
 * here. The alternative is copying the entire image to concatenate them, which
 * for IDAT means doubling peak memory to compute a four-byte number.
 *
 * The value passed in and returned is the RAW running register, not a finished
 * checksum: start from 0xFFFFFFFF and complement the final result. Feeding a
 * complemented value back in silently produces a wrong checksum. */
static unsigned long rs_crc32_update(unsigned long crc,
                                     const unsigned char *buf, size_t len)
{
    if (!rs_crc_table_ready) rs_crc_build_table();
    for (size_t i = 0; i < len; i++) {
        crc = rs_crc_table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    }
    return crc;
}

/* Adler-32 over the raw stream. The modulus is the largest prime below 65536,
 * which is what makes the two running sums independent enough to catch the
 * byte-swap errors a plain sum would miss. */
static unsigned long rs_adler32(const unsigned char *buf, size_t len)
{
    unsigned long a = 1, b = 0;
    for (size_t i = 0; i < len; i++) {
        a = (a + buf[i]) % 65521UL;
        b = (b + a) % 65521UL;
    }
    return (b << 16) | a;
}

/* ---------------------------------------------------------------------------
 * Chunk writing.
 * ------------------------------------------------------------------------- */

/* Big-endian 32-bit store. PNG is big-endian throughout, on every field, which
 * is worth stating once here because it is the opposite of the byte order every
 * other binary this project writes uses. */
static void rs_put_u32(unsigned char *p, unsigned long v)
{
    p[0] = (unsigned char)((v >> 24) & 0xFF);
    p[1] = (unsigned char)((v >> 16) & 0xFF);
    p[2] = (unsigned char)((v >>  8) & 0xFF);
    p[3] = (unsigned char)( v        & 0xFF);
}

/* Emit one chunk: length, four-character type, payload, then the CRC computed
 * over the TYPE AND PAYLOAD but not over the length. Returns zero on a write
 * error so the caller can report the path. */
static int rs_png_chunk(FILE *f, const char *type,
                        const unsigned char *data, size_t len)
{
    unsigned char hdr[8];
    rs_put_u32(hdr, (unsigned long)len);
    memcpy(hdr + 4, type, 4);
    if (fwrite(hdr, 1, 8, f) != 8) return 0;
    if (len && fwrite(data, 1, len, f) != len) return 0;

    /* The CRC spans the type and the payload, which are not contiguous, so it
     * is accumulated across two calls and complemented once at the end. */
    unsigned long c = 0xFFFFFFFFUL;
    c = rs_crc32_update(c, (const unsigned char *)type, 4);
    c = rs_crc32_update(c, data, len);
    c ^= 0xFFFFFFFFUL;

    unsigned char crc[4];
    rs_put_u32(crc, c);
    return fwrite(crc, 1, 4, f) == 4;
}

/* ---------------------------------------------------------------------------
 * The encoder.
 * ------------------------------------------------------------------------- */

resonarsat_status_t rs_png_write(const char *path, const unsigned char *pixels,
                                 size_t n_col, size_t n_row, int channels)
{
    if (!path || !pixels || n_col == 0 || n_row == 0) return RS_ERR_ARG;
    if (channels != 1 && channels != 3) {
        rs_set_error("png: %d channels is not 1 (grey) or 3 (RGB)", channels);
        return RS_ERR_ARG;
    }

    /* Guard the row-stride and total-size arithmetic before allocating. An
     * image large enough to overflow size_t here would otherwise allocate a
     * small buffer and then be written past. */
    const size_t stride = n_col * (size_t)channels;
    if (stride / (size_t)channels != n_col) return RS_ERR_ARG;
    const size_t raw_len = (stride + 1) * n_row;          /* +1 filter byte/row */
    if ((stride + 1) != 0 && raw_len / (stride + 1) != n_row) return RS_ERR_ARG;

    /* The filtered raw stream: every row prefixed with filter type 0. */
    unsigned char *raw = malloc(raw_len);
    if (!raw) return RS_ERR_ALLOC;
    for (size_t r = 0; r < n_row; r++) {
        raw[r * (stride + 1)] = 0;
        memcpy(raw + r * (stride + 1) + 1, pixels + r * stride, stride);
    }

    FILE *f = fopen(path, "wb");
    if (!f) {
        free(raw);
        rs_set_error("png: cannot open %s for writing", path);
        return RS_ERR_IO;
    }

    int ok = 1;

    static const unsigned char sig[8] = { 137, 'P', 'N', 'G', '\r', '\n', 26, '\n' };
    if (fwrite(sig, 1, 8, f) != 8) ok = 0;

    /* IHDR. Colour type 0 is greyscale and 2 is truecolour; bit depth 8, no
     * interlacing, and compression method 0, which is the only one PNG defines. */
    if (ok) {
        unsigned char ihdr[13];
        rs_put_u32(ihdr,     (unsigned long)n_col);
        rs_put_u32(ihdr + 4, (unsigned long)n_row);
        ihdr[8]  = 8;
        ihdr[9]  = (channels == 1) ? 0 : 2;
        ihdr[10] = 0;
        ihdr[11] = 0;
        ihdr[12] = 0;
        if (!rs_png_chunk(f, "IHDR", ihdr, sizeof ihdr)) ok = 0;
    }

    /* IDAT. A zlib stream is a two-byte header, the DEFLATE data, and the
     * Adler-32 of the uncompressed bytes. 0x78 0x01 says deflate with a 32 KiB
     * window at the fastest setting; the pair is chosen so that the 16-bit value
     * is a multiple of 31, which is the check zlib readers apply. */
    if (ok) {
        const size_t n_blocks = (raw_len + 65534) / 65535;
        const size_t z_len = 2 + n_blocks * 5 + raw_len + 4;
        unsigned char *z = malloc(z_len);
        if (!z) {
            ok = 0;
        } else {
            size_t o = 0;
            z[o++] = 0x78;
            z[o++] = 0x01;

            size_t left = raw_len, off = 0;
            do {
                const size_t n = (left > 65535) ? 65535 : left;
                const int final = (left - n == 0);
                /* Stored block: one header byte whose low bit is BFINAL and
                 * whose next two bits are BTYPE 00, then the length and its
                 * one's complement, both little-endian -- the only place in the
                 * format that is not big-endian, because DEFLATE predates it. */
                z[o++] = (unsigned char)(final ? 1 : 0);
                z[o++] = (unsigned char)( n        & 0xFF);
                z[o++] = (unsigned char)((n >> 8)  & 0xFF);
                z[o++] = (unsigned char)((~n)       & 0xFF);
                z[o++] = (unsigned char)((~n >> 8)  & 0xFF);
                memcpy(z + o, raw + off, n);
                o += n; off += n; left -= n;
            } while (left > 0);

            rs_put_u32(z + o, rs_adler32(raw, raw_len));
            o += 4;

            if (!rs_png_chunk(f, "IDAT", z, o)) ok = 0;
            free(z);
        }
    }

    if (ok && !rs_png_chunk(f, "IEND", NULL, 0)) ok = 0;
    if (ferror(f)) ok = 0;

    if (fclose(f) != 0) ok = 0;
    free(raw);

    if (!ok) {
        rs_set_error("png: write failed for %s", path);
        return RS_ERR_IO;
    }
    return RS_OK;
}
