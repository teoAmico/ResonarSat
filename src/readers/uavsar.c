/* UAVSAR SLC reader: flat binary samples plus an ASCII annotation file.
 *
 * UAVSAR is the reader-development starter for this project because its format
 * needs no external library: the image is an uncompressed array of interleaved
 * float pairs, and everything needed to interpret it is in a plain-text .ann
 * file of "key = value ; comment" lines.
 *
 * The error contract from the implementation plan applies here in full. This
 * parser sees untrusted external files, so every field that determines an
 * allocation or an offset is validated before use, and a truncated or corrupt
 * file yields a described status rather than a crash or a silently short read. */

#include "resonarsat/readers.h"
#include "resonarsat/geom.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Strip leading and trailing whitespace from a string in place, returning a
 * pointer to the first non-space character. */
static char *rs_trim(char *s)
{
    while (*s && isspace((unsigned char)*s)) s++;
    if (!*s) return s;
    char *end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

/* Look up one key in a UAVSAR annotation file.
 *
 * Annotation lines take the form "Key Name (units) = value ; comment". Matching
 * is done on the portion of the key before any parenthesised unit, so callers
 * pass the bare name and need not reproduce the unit text exactly. The first
 * match wins; comments after ';' are discarded.
 *
 * Returns RS_ERR_MISSING_META, naming the key, if it does not appear. The
 * caller decides whether that is fatal: some fields are optional. */
static resonarsat_status_t rs_ann_lookup(const char *path, const char *key,
                                         char *value, size_t value_len)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        rs_set_error("uavsar: cannot open annotation %s", path);
        return RS_ERR_IO;
    }

    char line[1024];
    const size_t key_len = strlen(key);

    while (fgets(line, sizeof line, f)) {
        char *semi = strchr(line, ';');
        if (semi) *semi = '\0';

        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';

        char *name = rs_trim(line);
        /* Truncate at a parenthesised unit so "Range Spacing (m)" matches
         * a lookup for "Range Spacing". */
        char *paren = strchr(name, '(');
        if (paren) { *paren = '\0'; name = rs_trim(name); }

        if (strncasecmp(name, key, key_len) == 0 && strlen(name) == key_len) {
            char *val = rs_trim(eq + 1);
            snprintf(value, value_len, "%s", val);
            fclose(f);
            return RS_OK;
        }
    }

    fclose(f);
    rs_set_error("uavsar: annotation %s has no '%s' field", path, key);
    return RS_ERR_MISSING_META;
}

/* Look up a numeric annotation field. Returns RS_ERR_FORMAT if the value is
 * present but does not parse as a number, which is a different failure from the
 * field being absent and is worth distinguishing in the message. */
static resonarsat_status_t rs_ann_lookup_double(const char *path, const char *key, double *out)
{
    char buf[256];
    resonarsat_status_t st = rs_ann_lookup(path, key, buf, sizeof buf);
    if (st != RS_OK) return st;

    char *end = NULL;
    const double v = strtod(buf, &end);
    if (end == buf) {
        rs_set_error("uavsar: annotation field '%s' has non-numeric value '%s'", key, buf);
        return RS_ERR_FORMAT;
    }
    *out = v;
    return RS_OK;
}

/* Read a UAVSAR SLC and its annotation into an image. */
resonarsat_status_t rs_read_uavsar(const char *slc_path, const char *ann_path, rs_slc_t *img)
{
    if (!slc_path || !ann_path || !img) return RS_ERR_ARG;

    memset(img, 0, sizeof *img);

    /* Dimensions first: everything else depends on them and they bound the
     * allocation, so they are validated before a single sample is read. */
    double d_rows = 0.0, d_cols = 0.0;
    resonarsat_status_t st;
    if ((st = rs_ann_lookup_double(ann_path, "slc_1_1x1 Rows", &d_rows)) != RS_OK &&
        (st = rs_ann_lookup_double(ann_path, "Number of SLC Lines", &d_rows)) != RS_OK) {
        return st;
    }
    if ((st = rs_ann_lookup_double(ann_path, "slc_1_1x1 Columns", &d_cols)) != RS_OK &&
        (st = rs_ann_lookup_double(ann_path, "Number of SLC Range Bins", &d_cols)) != RS_OK) {
        return st;
    }

    if (!(d_rows >= 1.0 && d_rows <= 1e7) || !(d_cols >= 1.0 && d_cols <= 1e7)) {
        rs_set_error("uavsar: annotation declares implausible dimensions %g x %g",
                     d_rows, d_cols);
        return RS_ERR_FORMAT;
    }

    const size_t n_az = (size_t)d_rows, n_rg = (size_t)d_cols;

    /* Cross-check the declared dimensions against the actual file size before
     * allocating. A mismatch here is the single most common way a corrupt or
     * mislabelled product manifests, and catching it now turns a would-be
     * out-of-bounds read into a clear message. */
    FILE *f = fopen(slc_path, "rb");
    if (!f) {
        rs_set_error("uavsar: cannot open %s", slc_path);
        return RS_ERR_IO;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        rs_set_error("uavsar: cannot determine size of %s", slc_path);
        return RS_ERR_IO;
    }
    const long file_size = ftell(f);
    rewind(f);

    const size_t expect = n_az * n_rg * 2 * sizeof(float);
    if (file_size < 0 || (size_t)file_size < expect) {
        fclose(f);
        rs_set_error("uavsar: %s is %ld bytes, but %zux%zu complex samples need %zu",
                     slc_path, file_size, n_az, n_rg, expect);
        return RS_ERR_FORMAT;
    }

    if ((st = rs_slc_alloc(img, n_az, n_rg)) != RS_OK) { fclose(f); return st; }

    /* Read row by row rather than in one call, so a truncated file is caught at
     * the row that fails instead of after a multi-gigabyte read. */
    float *row = malloc(n_rg * 2 * sizeof *row);
    if (!row) { fclose(f); rs_slc_free(img); return RS_ERR_ALLOC; }

    for (size_t a = 0; a < n_az; a++) {
        if (fread(row, sizeof *row, n_rg * 2, f) != n_rg * 2) {
            free(row);
            fclose(f);
            rs_slc_free(img);
            rs_set_error("uavsar: %s truncated at azimuth line %zu of %zu",
                         slc_path, a, n_az);
            return RS_ERR_FORMAT;
        }
        for (size_t r = 0; r < n_rg; r++) {
            img->data[a * n_rg + r] = row[2 * r] + row[2 * r + 1] * I;
        }
    }
    free(row);
    fclose(f);

    /* Metadata. Azimuth timing comes from the pulse interval, never from a
     * field named "PRF"; see the comment on rs_slc_t. UAVSAR annotations give
     * the along-track spacing and platform speed, from which the line interval
     * follows directly and consistently. */
    double az_spacing = 0.0, rg_spacing = 0.0, v_platform = 0.0;
    double centre_freq = 0.0, r0 = 0.0, inc = 0.0;

    if ((st = rs_ann_lookup_double(ann_path, "1x1 SLC Azimuth Pixel Spacing",
                                   &az_spacing)) != RS_OK &&
        (st = rs_ann_lookup_double(ann_path, "Azimuth Spacing", &az_spacing)) != RS_OK) {
        rs_slc_free(img);
        return st;
    }
    if ((st = rs_ann_lookup_double(ann_path, "1x1 SLC Range Pixel Spacing",
                                   &rg_spacing)) != RS_OK &&
        (st = rs_ann_lookup_double(ann_path, "Range Spacing", &rg_spacing)) != RS_OK) {
        rs_slc_free(img);
        return st;
    }

    /* These three are optional: absent ones leave sensible fallbacks rather
     * than failing, since a product missing platform speed is still usable for
     * anything that does not need absolute geometry. */
    if (rs_ann_lookup_double(ann_path, "Average Altitude", &r0) != RS_OK) r0 = 0.0;
    if (rs_ann_lookup_double(ann_path, "Global Average Terrain Height", &inc) != RS_OK) inc = 0.0;
    if (rs_ann_lookup_double(ann_path, "Average Along Track Velocity", &v_platform) != RS_OK)
        v_platform = 0.0;
    if (rs_ann_lookup_double(ann_path, "Center Frequency", &centre_freq) != RS_OK)
        centre_freq = 1.2575e9;   /* UAVSAR L-band nominal */

    /* UAVSAR annotations state the centre frequency in MHz. */
    if (centre_freq < 1e6) centre_freq *= 1e6;

    img->fc = centre_freq;
    img->az_spacing_m = az_spacing;
    img->rg_spacing_m = rg_spacing;
    img->v_platform = v_platform;
    img->r0 = (r0 > 0.0) ? r0 : 0.0;

    /* Derive the azimuth line interval from spacing and speed. If the platform
     * speed is unavailable, fall back to a nominal UAVSAR value so the image is
     * still processable, and say so in the source string so it is traceable. */
    double v_used = v_platform;
    int assumed_speed = 0;
    if (!(v_used > 0.0)) { v_used = 220.0; assumed_speed = 1; }

    if (!(az_spacing > 0.0)) {
        rs_slc_free(img);
        rs_set_error("uavsar: azimuth spacing %g m is not positive", az_spacing);
        return RS_ERR_FORMAT;
    }

    img->azimuth_time_interval = az_spacing / v_used;
    img->t_dwell = (double)n_az * img->azimuth_time_interval;

    if ((st = rs_slc_finalise_metadata(img)) != RS_OK) { rs_slc_free(img); return st; }

    snprintf(img->source, sizeof img->source, "UAVSAR%s",
             assumed_speed ? " (assumed platform speed)" : "");

    if ((st = rs_slc_validate(img)) != RS_OK) {
        /* Validation failure is reported but not fatal for a reader: the
         * message tells the operator what looks wrong, and the CLI's info
         * command surfaces it. Refusing to load would make diagnosing a
         * questionable product harder, not easier. */
        rs_slc_free(img);
        return st;
    }

    return RS_OK;
}
