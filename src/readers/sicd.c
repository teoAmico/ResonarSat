/* SICD reader: a focused complex image in a NITF 2.1 container.
 *
 * This is the format the published method actually starts from. The source's
 * own flowchart runs SLC image -> two-dimensional transform -> band-pass filter
 * -> inverse transform -> pixel tracking, so a focused product feeding
 * rs_subaperture_split() reproduces the described processing directly. The CPHD
 * route in this project, which focuses shifted pulse windows from phase history,
 * is a refinement of that rather than a reproduction of it.
 *
 * NITF is fixed-width ASCII fields in a rigid order, with no key names and no
 * delimiters: every offset is implied by the widths of everything before it.
 * That makes it unforgiving of a miscount and easy to validate, since the
 * declared header length must equal the number of bytes the fields consume. The
 * reader checks exactly that, which is how the field table below was verified
 * against a real product rather than against the specification alone.
 *
 * Layout, all offsets derived rather than assumed:
 *
 *     [0, HL)                     file header
 *     [HL, HL+LISH)               image subheader
 *     [HL+LISH, +LI)              pixel data
 *     [HL+LISH+LI, +LDSH)         data extension subheader
 *     [.. +LD)                    the SICD XML
 *
 * The XML sits at the END of the file, after the pixel data. On a multi-gigabyte
 * product that means the metadata cannot be read without seeking past the image,
 * which is why the geometry is parsed before a single sample is loaded.
 *
 * One transposition matters. SICD's Row axis is range and its Col axis is
 * azimuth, and the pixels are stored row-major, so the file is ordered
 * [range][azimuth] while rs_slc_t is [azimuth][range]. Loading it verbatim would
 * put cross-track data on the azimuth axis, which is a defect this project has
 * already made once and which hides well: the image still looks like a scene. */

#include "resonarsat/geocode.h"
#include "resonarsat/readers.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RS_NITF_HDR_MAX  (1u << 20)   /* file header cannot plausibly exceed this */
#define RS_SICD_XML_MAX  (16u << 20)

/* Read 'n' ASCII digits at 'p' as an unsigned value.
 *
 * NITF numeric fields are zero-padded decimal with no sign and no terminator, so
 * strtoul on the surrounding buffer would run past the field. Returns 0 and sets
 * '*ok' to zero if any character is not a digit, which is how a miscounted
 * offset announces itself rather than silently yielding a plausible number. */
static unsigned long rs_nitf_num(const unsigned char *p, size_t n, int *ok)
{
    unsigned long v = 0;
    for (size_t i = 0; i < n; i++) {
        if (p[i] < '0' || p[i] > '9') { *ok = 0; return 0; }
        v = v * 10u + (unsigned long)(p[i] - '0');
    }
    return v;
}

/* Decode a big-endian IEEE float. NITF stores numeric pixel data most
 * significant byte first regardless of host order. */
static float rs_be_float(const unsigned char *p, int swap)
{
    uint32_t u;
    float f;
    memcpy(&u, p, sizeof u);
    if (swap) {
        u = ((u & 0x000000ffU) << 24) | ((u & 0x0000ff00U) << 8) |
            ((u & 0x00ff0000U) >> 8)  | ((u & 0xff000000U) >> 24);
    }
    memcpy(&f, &u, sizeof f);
    return f;
}

/* Decode a big-endian signed 16-bit sample, as RE16I_IM16I stores I and Q.
 *
 * The integers are returned unscaled. SICD gives no scale factor for this pixel
 * type -- unlike AMP8I_PHS8I, which carries an amplitude table -- so there is no
 * radiometric constant to apply and inventing one would put a fabricated number
 * where a measured one belongs. Nothing downstream is harmed by that: this
 * project treats amplitudes as qualitative throughout, and every quantity it
 * reports is a frequency, a displacement or a correlation, all of which are
 * invariant to a constant scaling of the image. */
static float rs_be_i16(const unsigned char *p, int swap)
{
    uint16_t u;
    memcpy(&u, p, sizeof u);
    if (swap) u = (uint16_t)(((u & 0x00ffU) << 8) | ((u & 0xff00U) >> 8));
    return (float)(int16_t)u;
}

/* Text of the first <tag> in 'xml', or 0 if absent. See the note in cphd.c on
 * why a tag scan rather than a parser: the same handful of leaf elements, and
 * the same reluctance to take on an XML dependency for them. */
static int rs_sicd_tag(const char *xml, const char *tag, double *out)
{
    char open[64];
    if ((size_t)snprintf(open, sizeof open, "<%s>", tag) >= sizeof open) return 0;
    const char *p = strstr(xml, open);
    if (!p) return 0;
    p += strlen(open);
    char buf[64];
    size_t n = 0;
    while (n < sizeof buf - 1 && p[n] && p[n] != '<') { buf[n] = p[n]; n++; }
    buf[n] = '\0';
    char *end = NULL;
    const double v = strtod(buf, &end);
    if (end == buf) return 0;
    *out = v;
    return 1;
}

/* As rs_sicd_tag(), but for the n-th occurrence, counting from zero. SICD reuses
 * element names across the Row and Col grid descriptions, so position is the
 * only way to tell a range quantity from an azimuth one. */
static int rs_sicd_tag_nth(const char *xml, const char *tag, size_t nth, double *out)
{
    char open[64];
    if ((size_t)snprintf(open, sizeof open, "<%s>", tag) >= sizeof open) return 0;
    const char *p = xml;
    for (size_t i = 0; i <= nth; i++) {
        p = strstr(p, open);
        if (!p) return 0;
        p += strlen(open);
    }
    char buf[64];
    size_t n = 0;
    while (n < sizeof buf - 1 && p[n] && p[n] != '<') { buf[n] = p[n]; n++; }
    buf[n] = '\0';
    char *end = NULL;
    const double v = strtod(buf, &end);
    if (end == buf) return 0;
    *out = v;
    return 1;
}

/* Read a SICD product, optionally stopping before the pixels.
 *
 * 'meta_only' skips the allocation and the pixel pass, leaving 'img->data' NULL
 * while every metadata field is filled exactly as a full read would fill it.
 * The dimensions are still reported, because a caller asking only about
 * geometry still needs to know how large the image is.
 *
 * The split exists because a spotlight SICD runs to tens of thousands of samples
 * on a side -- 40000 x 40000 is 12.8 GB as complex float -- and answering a
 * question about the collect's geometry should not cost that. The pixel pass is
 * the only part that does, so it is the only part guarded. */
static resonarsat_status_t rs_sicd_read(const char *path, rs_slc_t *img,
                                        int meta_only)
{
    if (!path || !img) return RS_ERR_ARG;
    memset(img, 0, sizeof *img);

    const uint16_t probe = 1;
    const int swap = (*(const unsigned char *)&probe == 1);

    resonarsat_status_t st = RS_ERR_FORMAT;
    char *xml = NULL;
    unsigned char *row = NULL;

    FILE *f = fopen(path, "rb");
    if (!f) {
        rs_set_error("sicd: cannot open %s", path);
        return RS_ERR_IO;
    }
    if (fseek(f, 0, SEEK_END) != 0) { rs_set_error("sicd: %s is not seekable", path); goto done; }
    const long fsize = ftell(f);
    rewind(f);

    unsigned char hdr[1024];
    if (fread(hdr, 1, sizeof hdr, f) != sizeof hdr) {
        rs_set_error("sicd: %s is too short to hold a NITF header", path);
        goto done;
    }
    if (memcmp(hdr, "NITF02.10", 9) != 0) {
        rs_set_error("sicd: %s is not a NITF 2.1 file", path);
        goto done;
    }

    /* Fixed-width walk. The widths up to FL are constant for NITF 2.1; rather
     * than list forty security fields, the total is taken as the constant it is
     * and then checked against the declared header length below. */
    int ok = 1;
    const size_t OFF_FL = 342;
    const unsigned long fl = rs_nitf_num(hdr + OFF_FL, 12, &ok);
    const unsigned long hl = rs_nitf_num(hdr + OFF_FL + 12, 6, &ok);
    size_t off = OFF_FL + 18;
    const unsigned long numi = rs_nitf_num(hdr + off, 3, &ok);
    off += 3;
    if (!ok || hl == 0 || hl > RS_NITF_HDR_MAX || numi == 0) {
        rs_set_error("sicd: implausible NITF header (length %lu, %lu images)", hl, numi);
        goto done;
    }
    if (fl != (unsigned long)fsize) {
        rs_set_error("sicd: header declares %lu bytes but %s is %ld", fl, path, fsize);
        goto done;
    }

    /* Only the first image segment is read. SICD products carry one; a file with
     * several is a valid NITF this build does not interpret, and saying so beats
     * silently reading a fraction of it. */
    const unsigned long lish = rs_nitf_num(hdr + off, 6, &ok);
    const unsigned long li   = rs_nitf_num(hdr + off + 6, 10, &ok);
    off += 16u * (size_t)numi;
    if (numi > 1) {
        rs_set_error("sicd: %lu image segments; this build reads single-segment products",
                     numi);
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }
    off += 9;                                   /* NUMS, NUMX, NUMT */
    const unsigned long numdes = rs_nitf_num(hdr + off, 3, &ok);
    off += 3;
    unsigned long ldsh = 0, ld = 0;
    if (numdes >= 1) {
        ldsh = rs_nitf_num(hdr + off, 4, &ok);
        ld   = rs_nitf_num(hdr + off + 4, 9, &ok);
    }
    off += 13u * (size_t)numdes;
    off += 13;                                  /* NUMRES, UDHDL, XHDL */

    if (!ok) { rs_set_error("sicd: non-numeric field in the NITF header"); goto done; }
    /* The arithmetic check that validates the whole field table: the walk must
     * land exactly on the declared header length. A miscount anywhere shows up
     * here rather than as a wrong image later. */
    if (off != (size_t)hl) {
        rs_set_error("sicd: header field walk consumed %zu bytes, header declares %lu",
                     off, hl);
        goto done;
    }
    if (numdes < 1 || ld == 0 || ld > RS_SICD_XML_MAX) {
        rs_set_error("sicd: no usable data extension segment; SICD metadata is required");
        st = RS_ERR_MISSING_META;
        goto done;
    }

    /* ---- image subheader ---- */
    unsigned char ish[1024];
    if (lish < 400 || lish > sizeof ish ||
        fseek(f, (long)hl, SEEK_SET) != 0 ||
        fread(ish, 1, lish, f) != lish) {
        rs_set_error("sicd: cannot read the %lu-byte image subheader", lish);
        goto done;
    }
    /* Offsets within the image subheader, again fixed-width. */
    const size_t OFF_NROWS = 333;
    const unsigned long nrows  = rs_nitf_num(ish + OFF_NROWS,      8, &ok);
    const unsigned long ncols  = rs_nitf_num(ish + OFF_NROWS + 8,  8, &ok);
    const char pvtype = (char)ish[OFF_NROWS + 16];
    const unsigned long abpp = rs_nitf_num(ish + OFF_NROWS + 35, 2, &ok);

    if (!ok || nrows == 0 || ncols == 0) {
        rs_set_error("sicd: implausible image dimensions %lu x %lu", nrows, ncols);
        goto done;
    }
    /* Two sample layouts, distinguished by what the product declares rather than
     * by who made it.
     *
     * RE32F_IM32F is NITF PVTYPE 'R' at 32 bits and is what ICEYE ships.
     * RE16I_IM16I is PVTYPE 'SI' at 16 bits and is what Capella ships, so the
     * entire Capella SICD catalogue was unreadable while only the first was
     * accepted.
     *
     * Keyed on the declared type, never on CollectorName. The container says
     * which layout it holds and is right about it; a vendor-string test would
     * be wrong for a new vendor and would go on being wrong silently. That is
     * the failure mode item 3 of docs/FOLLOW-UPS.md records for the CPHD SGN
     * override, and it is not worth reproducing here to save one comparison. */
    size_t bps;                       /* bytes per complex sample */
    if (pvtype == 'R' && abpp == 32)       bps = 8;
    else if (pvtype == 'S' && abpp == 16)  bps = 4;
    else {
        rs_set_error("sicd: pixel type '%c' at %lu bits is not supported; this build "
                     "reads 32-bit IEEE float and 16-bit signed integer I/Q pairs",
                     pvtype, abpp);
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }
    /* Two bands of the declared width is the only layout consistent with the
     * declared segment size; checking it catches a compressed or blocked product
     * without having to interpret the compression fields. */
    if ((unsigned long long)nrows * ncols * (unsigned long long)bps
        != (unsigned long long)li) {
        rs_set_error("sicd: %lu x %lu complex samples need %llu bytes, segment holds %lu "
                     "-- the product is probably blocked or compressed",
                     nrows, ncols,
                     (unsigned long long)nrows * ncols * (unsigned long long)bps, li);
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }

    /* ---- SICD XML, which lives after the pixel data ---- */
    const long des_xml = (long)hl + (long)lish + (long)li + (long)ldsh;
    xml = malloc(ld + 1);
    if (!xml) { st = RS_ERR_ALLOC; goto done; }
    if (fseek(f, des_xml, SEEK_SET) != 0 || fread(xml, 1, ld, f) != ld) {
        rs_set_error("sicd: cannot read the %lu-byte SICD XML at offset %ld", ld, des_xml);
        goto done;
    }
    xml[ld] = '\0';

    /* Row is range and Col is azimuth, so the file's row-major order is
     * [range][azimuth] and the two spacings belong to opposite axes from the
     * ones their element names suggest. */
    double ss_rg = 0.0, ss_az = 0.0, fx_min = 0.0, fx_max = 0.0;
    double dur = 0.0, inc_deg = 0.0, slant = 0.0;
    (void)rs_sicd_tag_nth(xml, "SS", 0, &ss_rg);
    (void)rs_sicd_tag_nth(xml, "SS", 1, &ss_az);
    (void)rs_sicd_tag(xml, "Min", &fx_min);
    (void)rs_sicd_tag(xml, "Max", &fx_max);
    (void)rs_sicd_tag(xml, "CollectDuration", &dur);
    (void)rs_sicd_tag(xml, "IncidenceAng", &inc_deg);
    (void)rs_sicd_tag(xml, "SlantRange", &slant);

    /* The image plane. SICD gives the scene centre point in ECF and the Row and
     * Col directions as ECF unit vectors, which is exactly a plane definition.
     * Row is range and Col is azimuth, matching the transpose applied below, so
     * u_x is the azimuth axis and u_y the range axis of the stored image.
     *
     * This plane lies in the SLANT plane, unlike a CPHD scene plane which lies
     * in the ground plane. Recording the vectors rather than assuming an
     * orientation is what allows the two to be reconciled: the angle between
     * their range axes is the grazing angle, and a comparison that ignores it
     * measures the projection instead of the scene. */
    const char *scp = strstr(xml, "<SCP>");
    const char *grid = strstr(xml, "<Grid>");
    if (scp && grid) {
        double ox = 0.0, oy = 0.0, oz = 0.0;
        const char *ecf_el = strstr(scp, "<ECF>");
        if (ecf_el) {
            (void)rs_sicd_tag(ecf_el, "X", &ox);
            (void)rs_sicd_tag(ecf_el, "Y", &oy);
            (void)rs_sicd_tag(ecf_el, "Z", &oz);
        }
        const char *row_el = strstr(grid, "<Row>");
        const char *col_el = strstr(grid, "<Col>");
        double rx=0, ry=0, rz=0, cx2=0, cy2=0, cz2=0;
        if (row_el) {
            const char *u = strstr(row_el, "<UVectECF>");
            if (u) { (void)rs_sicd_tag(u,"X",&rx); (void)rs_sicd_tag(u,"Y",&ry);
                     (void)rs_sicd_tag(u,"Z",&rz); }
        }
        if (col_el) {
            const char *u = strstr(col_el, "<UVectECF>");
            if (u) { (void)rs_sicd_tag(u,"X",&cx2); (void)rs_sicd_tag(u,"Y",&cy2);
                     (void)rs_sicd_tag(u,"Z",&cz2); }
        }
        if ((ox != 0.0 || oy != 0.0 || oz != 0.0) &&
            (rx != 0.0 || ry != 0.0 || rz != 0.0) &&
            (cx2 != 0.0 || cy2 != 0.0 || cz2 != 0.0)) {
            img->plane.origin[0]=ox; img->plane.origin[1]=oy; img->plane.origin[2]=oz;
            img->plane.u_x[0]=cx2;  img->plane.u_x[1]=cy2;  img->plane.u_x[2]=cz2;
            img->plane.u_y[0]=rx;   img->plane.u_y[1]=ry;   img->plane.u_y[2]=rz;
            /* Shift the origin from the scene centre point to image sample
             * (0,0), so a caller can geolocate a pixel from its index times its
             * spacing without also needing to know where the SCP sits. SICD puts
             * the SCP near but not exactly at the image centre. */
            double scp_row = 0.0, scp_col = 0.0;
            const char *sp = strstr(xml, "<SCPPixel>");
            if (sp) {
                (void)rs_sicd_tag(sp, "Row", &scp_row);
                (void)rs_sicd_tag(sp, "Col", &scp_col);
            }
            double ss_r = 0.0, ss_c = 0.0;
            (void)rs_sicd_tag_nth(xml, "SS", 0, &ss_r);
            (void)rs_sicd_tag_nth(xml, "SS", 1, &ss_c);
            for (int k = 0; k < 3; k++) {
                img->plane.origin[k] -= scp_col * ss_c * img->plane.u_x[k]
                                      + scp_row * ss_r * img->plane.u_y[k];
            }
            /* SICD's grid is the slant plane, and the ARP at aperture centre
             * is what a projection to the ground has to look from. */
            const char *arp = strstr(xml, "<ARPPos>");
            if (arp) {
                (void)rs_sicd_tag(arp, "X", &img->plane.sensor[0]);
                (void)rs_sicd_tag(arp, "Y", &img->plane.sensor[1]);
                (void)rs_sicd_tag(arp, "Z", &img->plane.sensor[2]);
            }
            const char *arv = strstr(xml, "<ARPVel>");
            if (arv) {
                (void)rs_sicd_tag(arv, "X", &img->plane.sensor_vel[0]);
                (void)rs_sicd_tag(arv, "Y", &img->plane.sensor_vel[1]);
                (void)rs_sicd_tag(arv, "Z", &img->plane.sensor_vel[2]);
            }
            double scp_ecf[3] = { ox, oy, oz }, scp_h = 0.0;
            (void)rs_geo_ecf_to_llh(scp_ecf, NULL, NULL, &scp_h);
            img->plane.ref_hae = scp_h;
            img->plane.is_slant = 1;
            img->plane.valid = 1;
        }
    }

    double vx = 0.0, vy = 0.0, vz = 0.0;
    const char *vel = strstr(xml, "<ARPVel>");
    if (vel) {
        (void)rs_sicd_tag(vel, "X", &vx);
        (void)rs_sicd_tag(vel, "Y", &vy);
        (void)rs_sicd_tag(vel, "Z", &vz);
    }

    if (!(fx_min > 0.0 && fx_max > fx_min)) {
        rs_set_error("sicd: RadarCollection/TxFrequency is absent or degenerate");
        st = RS_ERR_MISSING_META;
        goto done;
    }

    /* ---- pixels, transposed on the way in ---- */
    if (meta_only) {
        /* rs_slc_alloc() would have set these. A caller asking about geometry
         * still needs them -- the memory and grid-width checks are stated in
         * samples -- so report the dimensions without buying the buffer. */
        img->n_az = (size_t)ncols;
        img->n_rg = (size_t)nrows;
    } else {
        st = rs_slc_alloc(img, (size_t)ncols, (size_t)nrows); /* n_az = Col, n_rg = Row */
        if (st != RS_OK) goto done;

        row = malloc((size_t)ncols * bps);
        if (!row) { st = RS_ERR_ALLOC; goto fail; }
        if (fseek(f, (long)hl + (long)lish, SEEK_SET) != 0) {
            rs_set_error("sicd: cannot seek to the pixel data");
            st = RS_ERR_IO;
            goto fail;
        }
        const size_t half = bps / 2u;
        for (size_t r = 0; r < (size_t)nrows; r++) {
            if (fread(row, bps, (size_t)ncols, f) != (size_t)ncols) {
                rs_set_error("sicd: %s truncated at range line %zu of %lu",
                             path, r, nrows);
                st = RS_ERR_FORMAT;
                goto fail;
            }
            for (size_t c = 0; c < (size_t)ncols; c++) {
                const unsigned char *s = row + c * bps;
                const float re = (bps == 8) ? rs_be_float(s, swap)
                                            : rs_be_i16(s, swap);
                const float im = (bps == 8) ? rs_be_float(s + half, swap)
                                            : rs_be_i16(s + half, swap);
                img->data[c * (size_t)nrows + r] = re + im * I;   /* [az][rg] */
            }
        }
    }

    img->fc = 0.5 * (fx_min + fx_max);
    img->rg_spacing_m = (ss_rg > 0.0) ? ss_rg : 1.0;
    img->az_spacing_m = (ss_az > 0.0) ? ss_az : 1.0;
    img->r0 = slant;
    img->incidence = inc_deg * M_PI / 180.0;
    img->v_platform = sqrt(vx * vx + vy * vy + vz * vz);
    img->t_dwell = dur;

    /* Azimuth timing for a FOCUSED product, which is not what it is for a raw
     * one. An earlier version set this to CollectDuration/NumCols, treating the
     * image columns as successive acquisition times. That is wrong: in a focused
     * image the azimuth axis is spatial, and in spotlight mode every column is
     * synthesised from the whole aperture, so no column corresponds to an
     * instant. Here the platform travels 7.7 km while the scene spans 4.9 km in
     * azimuth -- aperture length and scene extent are simply different
     * quantities.
     *
     * The meaningful rate for a spatially sampled grid is the equivalent PRF,
     * v/az_spacing. rs_slc_validate() cross-checks exactly that relation and
     * rejected the earlier derivation with a 1.56x disagreement, which is the
     * check earning its place. */
    if (img->v_platform > 0.0 && img->az_spacing_m > 0.0) {
        img->azimuth_time_interval = img->az_spacing_m / img->v_platform;
    }
    /* A focused product records no transmit PRF, and inventing one would put a
     * fabricated number where the pipeline reaches for a measured one. */
    img->pulse_prf = 0.0;

    /* Let the shared derivation set fs_az and lambda. A reader must never fill
     * those itself: doing so is how a product's "PRF" field ends up on the
     * azimuth sampling axis, which is the mistake slc.h exists to prevent. */
    st = rs_slc_finalise_metadata(img);
    if (st != RS_OK) {
        rs_set_error("sicd: %s lacks the azimuth timing or carrier frequency the "
                     "pipeline needs", path);
        goto fail;
    }
    snprintf(img->source, sizeof img->source, "SICD %lux%lu", ncols, nrows);

    st = RS_OK;
    goto done;

fail:
    rs_slc_free(img);
done:
    free(row);
    free(xml);
    fclose(f);
    return st;
}

/* Read a SICD product. See readers.h. */
resonarsat_status_t rs_read_sicd(const char *path, rs_slc_t *img)
{
    return rs_sicd_read(path, img, 0);
}

/* Read a SICD product's metadata without its pixels. See readers.h. */
resonarsat_status_t rs_read_sicd_meta(const char *path, rs_slc_t *img)
{
    return rs_sicd_read(path, img, 1);
}
