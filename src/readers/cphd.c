/* CPHD 1.x phase-history reader.
 *
 * CPHD is the format the sub-aperture strategy needs, because it is the only
 * one that still contains pulses. A focused product has already integrated the
 * whole aperture; recovering sub-looks from it means filtering after the fact,
 * which approximates the source's section 3.1 rather than implementing it. See
 * focus.h.
 *
 * The container is refreshingly plain for a defence-standard format: an ASCII
 * header of "KEY := VALUE" lines, terminated by a form feed, whose values are
 * byte offsets and sizes for three blocks.
 *
 *     CPHD/1.1.0
 *     XML_BLOCK_SIZE := 15887
 *     XML_BLOCK_BYTE_OFFSET := 1024
 *     PVP_BLOCK_SIZE := 1385520
 *     PVP_BLOCK_BYTE_OFFSET := 340224
 *     SIGNAL_BLOCK_SIZE := 173143816
 *     SIGNAL_BLOCK_BYTE_OFFSET := 1725760
 *
 * The XML block describes the collect; the PVP block carries per-pulse platform
 * state; the signal block carries the samples. This reader extracts only what
 * the pipeline uses, and validates every extent against the file's real size
 * before allocating anything, per the error contract in readers.h.
 *
 * The one substantive piece of signal processing here is range compression. A
 * CPHD signal array is normally in the FX domain -- samples are a function of
 * transmitted frequency, not of delay -- so the pulses are not range compressed
 * and rs_cphd_t is documented to require that they are. Transforming each pulse
 * along its sample axis produces the profile. Doing it in the reader rather
 * than downstream keeps chirp and sampling conventions where they belong, as a
 * format concern. */

#include "resonarsat/readers.h"
#include "resonarsat/fft.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define RS_CPHD_HDR_MAX   (64u * 1024u)   /* header must terminate within this */
#define RS_CPHD_XML_MAX   (16u * 1024u * 1024u)
#define RS_CPHD_PVP_WORDS 256u            /* sanity cap on NumBytesPVP/8 */

/* Speed of light used to convert delay to slant range, m/s. */
static const double RS_C = 299792458.0;

/* ------------------------------------------------------------------ */
/* Byte order                                                          */
/* ------------------------------------------------------------------ */

/* Return non-zero if this machine stores integers little-endian.
 *
 * CPHD binary data are big-endian without exception, so every multi-byte value
 * read from the PVP and signal blocks is swapped on the platforms this project
 * targets. Detecting at runtime rather than with preprocessor guesswork keeps
 * the reader correct on any host without a configure step. */
static int rs_host_is_le(void)
{
    const uint16_t probe = 1;
    return *(const unsigned char *)&probe == 1;
}

/* Reverse the bytes of a 64-bit quantity. */
static uint64_t rs_bswap64(uint64_t v)
{
    return ((v & 0x00000000000000ffULL) << 56) | ((v & 0x000000000000ff00ULL) << 40) |
           ((v & 0x0000000000ff0000ULL) << 24) | ((v & 0x00000000ff000000ULL) <<  8) |
           ((v & 0x000000ff00000000ULL) >>  8) | ((v & 0x0000ff0000000000ULL) >> 24) |
           ((v & 0x00ff000000000000ULL) >> 40) | ((v & 0xff00000000000000ULL) >> 56);
}

/* Reverse the bytes of a 32-bit quantity. */
static uint32_t rs_bswap32(uint32_t v)
{
    return ((v & 0x000000ffU) << 24) | ((v & 0x0000ff00U) << 8) |
           ((v & 0x00ff0000U) >> 8)  | ((v & 0xff000000U) >> 24);
}

/* Decode a big-endian IEEE double from an eight-byte buffer.
 *
 * Goes through memcpy rather than a cast because a double and a uint64_t may
 * not be accessed through each other's type, and because the PVP block offers
 * no alignment guarantee. */
static double rs_be_f8(const unsigned char *p, int swap)
{
    uint64_t u;
    double d;
    memcpy(&u, p, sizeof u);
    if (swap) u = rs_bswap64(u);
    memcpy(&d, &u, sizeof d);
    return d;
}

/* Big-endian signed 16-bit, for the CI4 sample format.
 *
 * CI4 stores each sample as two signed 16-bit integers rather than two floats,
 * which halves the file at the cost of dynamic range. It is what most vendors
 * ship; supporting only CF8 fitted this reader to the one collector it was
 * developed against. The values are left unscaled: CPHD carries no scale factor
 * for CI4, every stage downstream is relative, and inventing a normalisation
 * would put a fabricated constant in front of the data. */
static float rs_be_i2(const unsigned char *p, int swap)
{
    const unsigned lo = swap ? p[1] : p[0];
    const unsigned hi = swap ? p[0] : p[1];
    const int16_t v = (int16_t)((hi << 8) | lo);
    return (float)v;
}

/* Decode a big-endian IEEE float from a four-byte buffer, for CF8 samples. */
static float rs_be_f4(const unsigned char *p, int swap)
{
    uint32_t u;
    float f;
    memcpy(&u, p, sizeof u);
    if (swap) u = rs_bswap32(u);
    memcpy(&f, &u, sizeof f);
    return f;
}

/* ------------------------------------------------------------------ */
/* Header and XML scanning                                             */
/* ------------------------------------------------------------------ */

/* Every declared block extent from the ASCII header, plus the values this
 * reader needs from the XML. Kept in one struct so validation can consider the
 * whole picture before anything is allocated. */
typedef struct {
    unsigned long xml_off, xml_size;
    unsigned long pvp_off, pvp_size;
    unsigned long sig_off, sig_size;
} rs_cphd_blocks_t;

/* Read the value of a "KEY := VALUE" line from the ASCII header text.
 *
 * Returns 1 and writes the parsed unsigned long on success, 0 if the key is
 * absent. The header is small, bounded and already in memory, so a linear scan
 * per key is simpler than tokenising and costs nothing measurable. */
static int rs_hdr_ulong(const char *hdr, const char *key, unsigned long *out)
{
    const char *p = strstr(hdr, key);
    if (!p) return 0;
    p += strlen(key);
    while (*p == ' ' || *p == '\t') p++;
    if (p[0] != ':' || p[1] != '=') return 0;
    p += 2;
    char *end = NULL;
    const unsigned long v = strtoul(p, &end, 10);
    if (end == p) return 0;
    *out = v;
    return 1;
}

/* Locate the text content of the first <tag> in 'xml', searching from 'from'.
 *
 * Writes the bounds of the content into 'begin'/'end' and returns 1, or returns
 * 0 if the tag is absent. Namespace prefixes are not handled because CPHD
 * instances in the wild do not use them on these elements; a prefixed document
 * fails the missing-metadata check rather than being silently misread.
 *
 * This is a tag scan, not an XML parser. It is sufficient because the fifteen
 * values this pipeline needs are all simple leaf elements with unique names
 * within their scope, and because a real parser would be a dependency the
 * project does not otherwise carry. */
static int rs_xml_find(const char *xml, const char *tag, const char *from,
                       const char **begin, const char **end)
{
    char open[64];
    if ((size_t)snprintf(open, sizeof open, "<%s>", tag) >= sizeof open) return 0;
    const char *p = strstr(from ? from : xml, open);
    if (!p) return 0;
    p += strlen(open);
    const char *q = strchr(p, '<');
    if (!q) return 0;
    *begin = p;
    *end = q;
    return 1;
}

/* Extract a double-valued XML element. Returns 1 on success, 0 if absent. */
static int rs_xml_double(const char *xml, const char *tag, const char *from, double *out)
{
    const char *b, *e;
    if (!rs_xml_find(xml, tag, from, &b, &e)) return 0;
    char buf[64];
    const size_t n = (size_t)(e - b);
    if (n == 0 || n >= sizeof buf) return 0;
    memcpy(buf, b, n);
    buf[n] = '\0';
    char *stop = NULL;
    const double v = strtod(buf, &stop);
    if (stop == buf) return 0;
    *out = v;
    return 1;
}

/* Extract an unsigned XML element. Returns 1 on success, 0 if absent. */
static int rs_xml_ulong(const char *xml, const char *tag, const char *from, unsigned long *out)
{
    double d;
    if (!rs_xml_double(xml, tag, from, &d)) return 0;
    if (d < 0.0 || d > 4.0e9) return 0;
    *out = (unsigned long)d;
    return 1;
}

/* Compare an XML element's text against an expected literal.
 *
 * Returns 1 on match, 0 on mismatch, -1 if the element is absent, so that a
 * caller can tell "wrong value" from "no value" and report the right status. */
static int rs_xml_is(const char *xml, const char *tag, const char *from, const char *want)
{
    const char *b, *e;
    if (!rs_xml_find(xml, tag, from, &b, &e)) return -1;
    const size_t n = (size_t)(e - b), w = strlen(want);
    return (n == w && memcmp(b, want, w) == 0) ? 1 : 0;
}

/* Read the word offset of a named per-vector parameter.
 *
 * PVP layout is declared as a list of <Name><Offset>k</Offset>... entries where
 * k counts eight-byte words from the start of the vector. Returns 1 and writes
 * the offset on success, 0 if the parameter is not declared. */
static int rs_pvp_offset(const char *xml, const char *name, unsigned long *out)
{
    const char *b, *e;
    char open[64];
    if ((size_t)snprintf(open, sizeof open, "<%s>", name) >= sizeof open) return 0;
    const char *p = strstr(xml, open);
    if (!p) return 0;
    if (!rs_xml_find(xml, "Offset", p, &b, &e)) return 0;
    return rs_xml_ulong(xml, "Offset", p, out);
}

/* ------------------------------------------------------------------ */
/* Geometry                                                            */
/* ------------------------------------------------------------------ */

/* Read an ECF triple written as <X>..</X><Y>..</Y><Z>..</Z> following 'from'.
 * Returns 1 on success, 0 if any component is missing. */
static int rs_xml_ecf(const char *xml, const char *from, double v[3])
{
    return rs_xml_double(xml, "X", from, &v[0]) &&
           rs_xml_double(xml, "Y", from, &v[1]) &&
           rs_xml_double(xml, "Z", from, &v[2]);
}

/* Normalise a 3-vector in place, returning its original length. */
static double rs_norm3(double v[3])
{
    const double n = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n > 0.0) { v[0] /= n; v[1] /= n; v[2] /= n; }
    return n;
}

/* Cross product, c = a x b. */
static void rs_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/* Euclidean distance between two points. */
static double rs_dist3(const double a[3], const double b[3])
{
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/* ------------------------------------------------------------------ */
/* Reader                                                              */
/* ------------------------------------------------------------------ */

/* Read a CPHD 1.x file into the focusing container. See readers.h. */
resonarsat_status_t rs_read_cphd(const char *path, const rs_cphd_read_opts_t *opts,
                                 rs_cphd_t *cphd)
{
    if (!path || !cphd) return RS_ERR_ARG;
    memset(cphd, 0, sizeof *cphd);

    const int swap = rs_host_is_le();
    resonarsat_status_t st = RS_ERR_FORMAT;
    char *xml = NULL;
    unsigned char *pvp = NULL, *raw = NULL;
    unsigned long *valid = NULL;
    float complex *line = NULL;
    rs_fft_plan *plan = NULL;

    FILE *f = fopen(path, "rb");
    if (!f) {
        rs_set_error("cphd: cannot open %s", path);
        return RS_ERR_IO;
    }

    /* File size first: every declared extent is checked against it. */
    if (fseek(f, 0, SEEK_END) != 0) { rs_set_error("cphd: %s is not seekable", path); goto io; }
    const long fsize = ftell(f);
    if (fsize <= 0) { rs_set_error("cphd: %s is empty", path); goto done; }
    rewind(f);

    /* ---- ASCII header ---- */
    char hdr[RS_CPHD_HDR_MAX];
    const size_t hread = fread(hdr, 1, sizeof hdr - 1, f);
    hdr[hread] = '\0';
    if (hread < 16 || memcmp(hdr, "CPHD/", 5) != 0) {
        rs_set_error("cphd: %s does not begin with a CPHD version line", path);
        goto done;
    }
    /* The header terminates at a form feed, and everything before it is ASCII
     * by definition. Both halves of that matter. Truncating at the form feed
     * stops a value appearing later in the XML from being read as a header key;
     * stopping at the first byte that cannot appear in an ASCII header stops a
     * file with no terminator at all from being accepted because some byte deep
     * in the signal block happened to be 0x0C. A plain memchr over the read
     * buffer does the first but not the second, and a file whose header simply
     * ends would then parse -- the signal block is float data, so it supplies a
     * stray form feed roughly one time in every few hundred bytes. */
    size_t hend = 0;
    while (hend < hread && hdr[hend] != '\f') {
        const unsigned char c = (unsigned char)hdr[hend];
        if (c != '\n' && c != '\r' && c != '\t' && (c < 0x20 || c > 0x7e)) break;
        hend++;
    }
    if (hend >= hread || hdr[hend] != '\f') {
        rs_set_error("cphd: %s has no ASCII header terminator (stopped at byte %zu)",
                     path, hend);
        goto done;
    }
    hdr[hend] = '\0';

    rs_cphd_blocks_t b;
    memset(&b, 0, sizeof b);
    if (!rs_hdr_ulong(hdr, "XML_BLOCK_BYTE_OFFSET", &b.xml_off) ||
        !rs_hdr_ulong(hdr, "XML_BLOCK_SIZE", &b.xml_size) ||
        !rs_hdr_ulong(hdr, "PVP_BLOCK_BYTE_OFFSET", &b.pvp_off) ||
        !rs_hdr_ulong(hdr, "PVP_BLOCK_SIZE", &b.pvp_size) ||
        !rs_hdr_ulong(hdr, "SIGNAL_BLOCK_BYTE_OFFSET", &b.sig_off) ||
        !rs_hdr_ulong(hdr, "SIGNAL_BLOCK_SIZE", &b.sig_size)) {
        rs_set_error("cphd: %s header is missing a block offset or size", path);
        goto done;
    }
    if (b.xml_size == 0 || b.xml_size > RS_CPHD_XML_MAX) {
        rs_set_error("cphd: XML block of %lu bytes is implausible", b.xml_size);
        goto done;
    }
    /* Each block must lie wholly inside the file. Written as a subtraction so
     * that a hostile offset near ULONG_MAX cannot wrap an addition. */
    const unsigned long fsz = (unsigned long)fsize;
    if (b.xml_off > fsz || b.xml_size > fsz - b.xml_off ||
        b.pvp_off > fsz || b.pvp_size > fsz - b.pvp_off ||
        b.sig_off > fsz || b.sig_size > fsz - b.sig_off) {
        rs_set_error("cphd: a declared block extends past the end of %s (%ld bytes)",
                     path, fsize);
        goto done;
    }

    /* ---- XML block ---- */
    xml = malloc(b.xml_size + 1);
    if (!xml) { st = RS_ERR_ALLOC; goto done; }
    if (fseek(f, (long)b.xml_off, SEEK_SET) != 0 ||
        fread(xml, 1, b.xml_size, f) != b.xml_size) {
        rs_set_error("cphd: short read of the XML block in %s", path);
        goto done;
    }
    xml[b.xml_size] = '\0';

    /* Sample format. CF8 is two big-endian float32; CI4 is two big-endian
     * int16. Both are complex pairs and differ only in width and type, so the
     * decode below branches once per sample rather than per format. Anything
     * else is refused rather than guessed at. */
    unsigned samp_bytes = 0;
    if (rs_xml_is(xml, "SignalArrayFormat", NULL, "CF8") == 1) {
        samp_bytes = 8u;
    } else if (rs_xml_is(xml, "SignalArrayFormat", NULL, "CI4") == 1) {
        samp_bytes = 4u;
    } else {
        rs_set_error("cphd: sample format is neither CF8 nor CI4");
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }
    const int domain_fx = rs_xml_is(xml, "DomainType", NULL, "FX");
    if (domain_fx < 0) {
        rs_set_error("cphd: Global/DomainType is absent");
        st = RS_ERR_MISSING_META;
        goto done;
    }
    if (rs_xml_is(xml, "CollectType", NULL, "MONOSTATIC") == 0) {
        rs_set_error("cphd: only monostatic collects are supported");
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }

    unsigned long n_chan = 1, n_vec = 0, n_samp = 0, pvp_bytes = 0;
    (void)rs_xml_ulong(xml, "NumCPHDChannels", NULL, &n_chan);
    if (n_chan != 1) {
        rs_set_error("cphd: %lu channels; only single-channel files are supported", n_chan);
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }
    /* NumVectors and NumSamples must be read from inside <Channel>, not from
     * the document at large. Both names recur elsewhere -- support arrays and
     * the image grid description carry their own -- and in a real product the
     * first occurrence is not the channel's. Taking the first match silently
     * yields dimensions that disagree with the signal block by an order of
     * magnitude, which the block-extent check then reports as a corrupt file. */
    const char *chan = strstr(xml, "<Channel>");
    if (!chan) {
        rs_set_error("cphd: Data/Channel is absent");
        st = RS_ERR_MISSING_META;
        goto done;
    }
    if (!rs_xml_ulong(xml, "NumVectors", chan, &n_vec) ||
        !rs_xml_ulong(xml, "NumSamples", chan, &n_samp) ||
        !rs_xml_ulong(xml, "NumBytesPVP", NULL, &pvp_bytes)) {
        rs_set_error("cphd: NumVectors, NumSamples or NumBytesPVP is absent");
        st = RS_ERR_MISSING_META;
        goto done;
    }
    if (n_vec == 0 || n_samp < 2 || pvp_bytes == 0 || pvp_bytes % 8 != 0 ||
        pvp_bytes / 8 > RS_CPHD_PVP_WORDS) {
        rs_set_error("cphd: implausible dimensions (%lu vectors, %lu samples, %lu PVP bytes)",
                     n_vec, n_samp, pvp_bytes);
        goto done;
    }
    /* The declared dimensions must account for exactly the declared blocks.
     * This is the check that catches a truncated download, which is otherwise
     * indistinguishable from a valid file until the last pulse reads short. */
    if ((unsigned long)pvp_bytes * n_vec > b.pvp_size) {
        rs_set_error("cphd: PVP block holds %lu bytes, %lu vectors need %lu",
                     b.pvp_size, n_vec, pvp_bytes * n_vec);
        goto done;
    }
    if (n_samp > b.sig_size / samp_bytes / n_vec) {
        rs_set_error("cphd: signal block holds %lu bytes, %lu x %lu %s samples need %llu",
                     b.sig_size, n_vec, n_samp,
                     samp_bytes == 8u ? "CF8" : "CI4",
                     (unsigned long long)n_vec * n_samp * samp_bytes);
        goto done;
    }

    /* Carrier. FxBand brackets the transmitted sweep; its centre is the
     * frequency the backprojector's phase term must remove. */
    double fx_min = 0.0, fx_max = 0.0;
    if (!rs_xml_double(xml, "FxMin", NULL, &fx_min) ||
        !rs_xml_double(xml, "FxMax", NULL, &fx_max) || fx_max <= fx_min) {
        rs_set_error("cphd: Global/FxBand is absent or degenerate");
        st = RS_ERR_MISSING_META;
        goto done;
    }
    /* The phase reference is the first sample's frequency, not the band centre.
     *
     * Transforming S[k] = A*exp(-j*2*pi*(SC0 + k*SCSS)*tau) along k puts the
     * peak at bin N*SCSS*tau, but leaves the amplitude multiplied by
     * exp(-j*2*pi*SC0*tau). That residual is what the backprojector's phase term
     * has to undo, so the carrier this container reports must be SC0 -- using
     * the band centre instead leaves an error of 4*pi*dR*(fc-SC0)/c, about
     * 3.4 radians for every metre a target sits from the reference point, which
     * defocuses the image completely while looking entirely reasonable. */
    double fc = 0.0;   /* set from SC0 once the PVP block is read, below */
    if (fx_max < 1.0e8 || fx_max > 1.0e11) {
        rs_set_error("cphd: transmit band %.4g-%.4g Hz is outside any radar band",
                     fx_min, fx_max);
        st = RS_ERR_RANGE;
        goto done;
    }

    /* Sign convention for the FX transform. */
    double sgn = -1.0;
    (void)rs_xml_double(xml, "SGN", NULL, &sgn);

    /* Scene frame: the file's own planar image-area axes, origin at the image
     * area reference point. Grid coordinates in focus.h are expressed here. */
    const char *iarp = strstr(xml, "<IARP>");
    const char *uiax = strstr(xml, "<uIAX>");
    const char *uiay = strstr(xml, "<uIAY>");
    double origin[3], ex[3], ey[3], ez[3];
    if (!iarp || !uiax || !uiay || !rs_xml_ecf(xml, iarp, origin) ||
        !rs_xml_ecf(xml, uiax, ex) || !rs_xml_ecf(xml, uiay, ey)) {
        rs_set_error("cphd: SceneCoordinates IARP or planar reference axes are absent");
        st = RS_ERR_MISSING_META;
        goto done;
    }
    rs_norm3(ex);
    rs_norm3(ey);
    rs_cross3(ex, ey, ez);
    if (rs_norm3(ez) < 1.0e-6) {
        rs_set_error("cphd: planar reference axes are parallel");
        st = RS_ERR_RANGE;
        goto done;
    }

    /* A NOTE ON HANDEDNESS, AND A CORRECTION THAT WAS TRIED AND WITHDRAWN.
     *
     * Capella products declare uIAX and uIAY LEFT-handed about local up:
     * (uIAX x uIAY) . up = -1.000 on every one of fourteen collects sampled
     * across six satellites, three modes and two months. It is a fixed
     * convention, and it is NOT a defect.
     *
     * It is tempting to "fix" it by rebuilding the cross-track axis to be
     * right-handed, and that was done here for a few hours before the evidence
     * caught up. It is wrong, for reasons worth leaving behind so that nobody
     * repeats it. Capella images BOTH sides, and the handedness does not change
     * with the side. Of those fourteen collects, nine are right-looking
     * and five left-looking -- computed from the platform track and the SRP,
     * with no file axis involved -- and all fourteen declare the same
     * left-handed pair. So handedness carries no information about geometry,
     * a correction keyed on it fires identically on every product, and it is
     * therefore equivalent to mirroring every Capella geolocation
     * unconditionally.
     *
     * The metadata is also self-consistent as declared, in three separate
     * ways. Projecting a file's own ImageAreaCornerPoints with its own axes
     * reproduces its own ImageArea exactly, in standard corner order. The XML
     * carries <SideOfTrack>, and on the reference collect it reads L, agreeing
     * with the look side computed independently from the platform track and
     * the SRP. And uIAY's bearing matches the declared <AzimuthAngle> to four
     * decimal places -- 231.6099 against 231.6145 -- which is the documented
     * convention: uIAY points from the scene toward the sensor, which is why
     * the platform comes out at positive Y.
     *
     * WHAT REMAINS UNEXPLAINED is the Giza imagery. There, three structures
     * 447 m apart -- matching the pyramids' 463 and 476 m spacing -- include
     * two differing by a factor of two in size, with chevron arms of 100-125 m
     * and 190-230 m against Menkaure's 105 m base and Khufu's 230 m. Projected
     * with the declared axes, the SMALL structure lands nearest Khufu and the
     * large one nearest Menkaure, which is backwards.
     *
     * That reading is itself suspect, and the suspicion is recorded so the
     * next person does not inherit it as fact: the structures were measured on
     * frames centred at the MEASURED positions rather than the projected ones,
     * so the projected Khufu sat at a frame edge and the chevron measured in
     * that frame was some 290 m away from it. It may simply be a different
     * structure. Settling this needs frames centred on the PROJECTED positions,
     * which has not been done. See data/README.md. */

    /* PVP word offsets. */
    unsigned long o_txtime = 0, o_txpos = 0, o_rcvpos = 0, o_srppos = 0, o_sc0 = 0, o_scss = 0;
    if (!rs_pvp_offset(xml, "TxTime", &o_txtime) || !rs_pvp_offset(xml, "TxPos", &o_txpos) ||
        !rs_pvp_offset(xml, "RcvPos", &o_rcvpos) || !rs_pvp_offset(xml, "SRPPos", &o_srppos) ||
        !rs_pvp_offset(xml, "SC0", &o_sc0) || !rs_pvp_offset(xml, "SCSS", &o_scss)) {
        rs_set_error("cphd: a required PVP parameter (TxTime/TxPos/RcvPos/SRPPos/SC0/SCSS) "
                     "is not declared");
        st = RS_ERR_MISSING_META;
        goto done;
    }
    /* SIGNAL is optional: when present it marks each vector as valid signal
     * data or not. Files in the wild really do carry invalid vectors, and they
     * carry them with NaN timestamps and positions, so ignoring the flag does
     * not merely process a little junk -- one NaN poisons the dwell span, the
     * PRF and every geometry derived from them. */
    unsigned long o_signal = 0;
    const int have_signal = rs_pvp_offset(xml, "SIGNAL", &o_signal);

    const unsigned long words = pvp_bytes / 8;
    if (have_signal && o_signal >= words) {
        rs_set_error("cphd: SIGNAL offset lies outside the %lu-word vector", words);
        goto done;
    }
    if (o_txtime >= words || o_txpos + 3 > words || o_rcvpos + 3 > words ||
        o_srppos + 3 > words || o_sc0 >= words || o_scss >= words) {
        rs_set_error("cphd: a PVP offset lies outside the %lu-word vector", words);
        goto done;
    }

    /* ---- PVP block ---- */
    pvp = malloc((size_t)pvp_bytes * n_vec);
    if (!pvp) { st = RS_ERR_ALLOC; goto done; }
    if (fseek(f, (long)b.pvp_off, SEEK_SET) != 0 ||
        fread(pvp, (size_t)pvp_bytes, n_vec, f) != n_vec) {
        rs_set_error("cphd: short read of the PVP block in %s", path);
        goto done;
    }

    /* Sample spacing in frequency is declared per pulse, and on real hardware it
     * genuinely wobbles: the reference collect varies by 2.8 parts per million
     * across its 63291 pulses. Backprojection needs one range grid for the whole
     * aperture, so the question is not whether the spacing is constant but
     * whether treating it as constant moves anything by an amount that matters.
     *
     * Measure the spread here; the tolerance is applied once the retained window
     * is known, because the error grows with distance from the profile centre
     * and a narrow window tolerates far more drift than a wide one. */
    const double scss0 = rs_be_f8(pvp + (size_t)o_scss * 8, swap);
    if (!(scss0 > 0.0)) {
        rs_set_error("cphd: first pulse declares a non-positive sample spacing");
        st = RS_ERR_RANGE;
        goto done;
    }
    double scss_lo = scss0, scss_hi = scss0;
    for (unsigned long p = 1; p < n_vec; p++) {
        const double s = rs_be_f8(pvp + (size_t)p * pvp_bytes + (size_t)o_scss * 8, swap);
        if (!isfinite(s) || s <= 0.0) continue;   /* invalid vectors carry NaN */
        if (s < scss_lo) scss_lo = s;
        if (s > scss_hi) scss_hi = s;
    }
    const double scss_spread = (scss_hi - scss_lo) / scss_lo;

    /* Capella declares SGN = +1 while shipping data the standard would call
     * SGN = -1; see the note at the transform below. The override is keyed on
     * the collector rather than applied globally so that a conformant CPHD
     * from anyone else still reads correctly. */
    double sgn_eff = sgn;
    {
        const char *b2, *e2;
        if (rs_xml_find(xml, "CollectorName", NULL, &b2, &e2)) {
            const size_t n2 = (size_t)(e2 - b2);
            for (size_t i = 0; i + 7 <= n2; i++) {
                if (strncasecmp(b2 + i, "capella", 7) == 0) {
                    sgn_eff = -sgn;
                    fprintf(stderr,
                        "note: Capella product -- its signal does not follow its "
                        "declared SGN=%+.0f, so the\n      FX-to-delay transform "
                        "is inverted to match Capella's own reference reader.\n"
                        "      Without this the image is mirrored in range.\n", sgn);
                    break;
                }
            }
        }
    }

    fc = rs_be_f8(pvp + (size_t)o_sc0 * 8, swap);
    if (fc < 1.0e8 || fc > 1.0e11) {
        rs_set_error("cphd: SC0 of %.4g Hz is outside any radar band", fc);
        st = RS_ERR_RANGE;
        goto done;
    }

    /* Range bin spacing follows from the transmitted bandwidth: a profile of
     * n_samp bins spans c/(2*SCSS) in slant range. */
    const double dr = RS_C / (2.0 * (double)n_samp * scss0);
    if (!(dr > 0.0) || dr > 1.0e4) {
        rs_set_error("cphd: derived range bin spacing %.4g m is not sane", dr);
        st = RS_ERR_RANGE;
        goto done;
    }

    /* ---- Valid vectors ----
     *
     * A vector is usable when the file does not mark it invalid and when every
     * quantity focusing needs is finite. Both halves are necessary: the SIGNAL
     * flag is optional and some files omit it, while a vector flagged invalid
     * carries NaN rather than stale numbers. Screening here, once, keeps every
     * later stage free of the question. */
    valid = malloc((size_t)n_vec * sizeof *valid);
    if (!valid) { st = RS_ERR_ALLOC; goto done; }
    size_t n_valid = 0, n_flagged = 0, n_nonfinite = 0;
    for (unsigned long p = 0; p < n_vec; p++) {
        const unsigned char *v = pvp + (size_t)p * pvp_bytes;
        if (have_signal) {
            int64_t flag;
            uint64_t u;
            memcpy(&u, v + (size_t)o_signal * 8, sizeof u);
            if (swap) u = rs_bswap64(u);
            memcpy(&flag, &u, sizeof flag);
            if (flag == 0) { n_flagged++; continue; }
        }
        int ok = isfinite(rs_be_f8(v + (size_t)o_txtime * 8, swap));
        for (int k = 0; k < 3 && ok; k++) {
            ok = isfinite(rs_be_f8(v + ((size_t)o_txpos  + (size_t)k) * 8, swap)) &&
                 isfinite(rs_be_f8(v + ((size_t)o_rcvpos + (size_t)k) * 8, swap)) &&
                 isfinite(rs_be_f8(v + ((size_t)o_srppos + (size_t)k) * 8, swap));
        }
        if (!ok) { n_nonfinite++; continue; }
        valid[n_valid++] = p;
    }
    if (n_flagged || n_nonfinite) {
        fprintf(stderr, "cphd: skipped %zu vector(s) flagged invalid and %zu with "
                        "non-finite geometry, of %lu\n", n_flagged, n_nonfinite, n_vec);
    }
    if (n_valid < 2) {
        rs_set_error("cphd: only %zu usable vector(s) of %lu", n_valid, n_vec);
        st = RS_ERR_FORMAT;
        goto done;
    }

    /* ---- Pulse and bin selection ---- */
    const size_t stride = (opts && opts->pulse_stride > 1) ? opts->pulse_stride : 1;
    size_t n_pulse = (n_valid + stride - 1) / stride;
    if (opts && opts->max_pulses && n_pulse > opts->max_pulses) n_pulse = opts->max_pulses;

    size_t n_rbin = (size_t)n_samp;
    if (opts && opts->rbin_window && opts->rbin_window < n_rbin) n_rbin = opts->rbin_window;
    if (n_rbin < 2 || n_pulse < 2) {
        rs_set_error("cphd: selection leaves %zu pulses of %zu bins, too few to focus",
                     n_pulse, n_rbin);
        st = RS_ERR_ARG;
        goto done;
    }

    /* Now the sample-spacing tolerance can be judged by its consequence. A
     * relative spread of 'scss_spread' displaces the outermost retained bin by
     * that fraction of its distance from the profile centre, so the worst-case
     * error is (n_rbin/2) * dr * spread metres. A tenth of a bin is the limit:
     * below it the common range grid is indistinguishable from per-pulse grids
     * after the linear interpolation backprojection already applies. */
    const double edge_err_bins = 0.5 * (double)n_rbin * scss_spread;
    if (edge_err_bins > 0.1) {
        rs_set_error("cphd: sample spacing varies by %.3g across pulses, displacing the "
                     "edge of a %zu-bin window by %.2f bins (%.3f m); a common range grid "
                     "cannot represent this collect. Try a narrower --rbins.",
                     scss_spread, n_rbin, edge_err_bins, edge_err_bins * dr);
        st = RS_ERR_UNSUPPORTED;
        goto done;
    }

    /* Refuse a container the machine cannot hold, before asking for it.
     *
     * A focused product is read one pulse at a time, but the container it lands
     * in is the whole collect: pulses times range bins times sixteen bytes. On
     * a large spotlight collect that is tens of gigabytes, and the failure mode
     * is not a NULL from malloc. Systems that overcommit hand back a pointer,
     * then kill the process when it touches the pages -- so the caller sees a
     * silent SIGKILL with no message, no output and no indication that a range
     * window was all that was needed.
     *
     * The check names the figure and the option that fixes it, because the
     * quantity that has to shrink is not obvious from the outside. */
    {
        const double want = (double)n_pulse * (double)n_rbin *
                            (double)sizeof(float complex);
        double have = 0.0;
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
        const long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGESIZE);
        if (pages > 0 && psz > 0) have = (double)pages * (double)psz;
#endif
        /* Half of physical memory: the caller still needs room for the focused
         * image, the sub-look stack and the working buffers. */
        if (have > 0.0 && want > 0.5 * have) {
            const size_t fits = (size_t)(0.5 * have /
                                ((double)n_pulse * (double)sizeof(float complex)));
            rs_set_error("cphd: %lu pulses x %zu range bins needs %.1f GB, and this "
                         "machine has %.1f GB. Pass --rbins %zu or fewer to read a "
                         "range window, or --pulse-stride to thin the pulses.",
                         n_pulse, n_rbin, want / 1e9, have / 1e9,
                         fits > 0 ? fits : 1);
            st = RS_ERR_ALLOC;
            goto done;
        }
    }

    st = rs_cphd_alloc(cphd, n_pulse, n_rbin);
    if (st != RS_OK) goto done;

    /* ---- Signal block, one pulse at a time ---- */
    raw  = malloc((size_t)n_samp * samp_bytes);   /* CF8 is 8, CI4 is 4 */
    line = malloc((size_t)n_samp * sizeof *line);
    if (!raw || !line) { st = RS_ERR_ALLOC; goto fail; }
    if (domain_fx == 1) {
        st = rs_fft_plan_create((size_t)n_samp, &plan);
        if (st != RS_OK) goto fail;
    }

    /* Zero delay is the scene reference point, because spotlight phase history
     * is motion compensated: the receive window follows the scene rather than
     * sitting at a fixed absolute range. The transform leaves zero delay at bin
     * 0 with negative delays wrapped into the top of the profile, so the window
     * is taken about bin 0 and rotated to put it at the centre -- which is the
     * position rs_cphd_t documents for r_ref. */
    const size_t half = n_rbin / 2u;

    for (size_t i = 0; i < n_pulse; i++) {
        const unsigned long p = (unsigned long)valid[i * stride];
        const unsigned char *v = pvp + (size_t)p * pvp_bytes;

        double tx[3], rcv[3], srp[3];
        for (int k = 0; k < 3; k++) {
            tx[k]  = rs_be_f8(v + ((size_t)o_txpos  + (size_t)k) * 8, swap);
            rcv[k] = rs_be_f8(v + ((size_t)o_rcvpos + (size_t)k) * 8, swap);
            srp[k] = rs_be_f8(v + ((size_t)o_srppos + (size_t)k) * 8, swap);
        }

        /* Monostatic phase centre, and the reference range the receive window
         * is centred on. */
        double ap[3];
        for (int k = 0; k < 3; k++) ap[k] = 0.5 * (tx[k] + rcv[k]);
        cphd->r_ref[i] = 0.5 * (rs_dist3(tx, srp) + rs_dist3(rcv, srp));

        /* Into the scene's planar frame. */
        const double d[3] = { ap[0] - origin[0], ap[1] - origin[1], ap[2] - origin[2] };
        cphd->pos[3 * i + 0] = d[0] * ex[0] + d[1] * ex[1] + d[2] * ex[2];
        cphd->pos[3 * i + 1] = d[0] * ey[0] + d[1] * ey[1] + d[2] * ey[2];
        cphd->pos[3 * i + 2] = d[0] * ez[0] + d[1] * ez[1] + d[2] * ez[2];
        cphd->t[i] = rs_be_f8(v + (size_t)o_txtime * 8, swap);

        const long off = (long)b.sig_off +
                         (long)((size_t)p * (size_t)n_samp * samp_bytes);
        if (fseek(f, off, SEEK_SET) != 0 ||
            fread(raw, samp_bytes, (size_t)n_samp, f) != (size_t)n_samp) {
            rs_set_error("cphd: short read of pulse %lu in %s", p, path);
            st = RS_ERR_FORMAT;
            goto fail;
        }
        if (samp_bytes == 8u) {
            for (size_t k = 0; k < (size_t)n_samp; k++) {
                line[k] = rs_be_f4(raw + k * 8u, swap)
                        + I * rs_be_f4(raw + k * 8u + 4u, swap);
            }
        } else {
            for (size_t k = 0; k < (size_t)n_samp; k++) {
                line[k] = rs_be_i2(raw + k * 4u, swap)
                        + I * rs_be_i2(raw + k * 4u + 2u, swap);
            }
        }

        /* Range compression. FX-domain samples are a function of transmitted
         * frequency, so one transform along the sample axis yields the delay
         * profile. SGN records the collector's convention: a negative sign
         * means the stored data are already a forward transform of the profile,
         * which the inverse undoes. */
        /* THE FX-TO-DELAY TRANSFORM DIRECTION, AND A VENDOR THAT MISLABELS IT.
         *
         * The standard is unambiguous. SARPy, NGA's own reference
         * implementation, documents Global/SGN as
         *
         *     Phase(fx) = SGN * fx * dTOA_TGT,   phase in cycles
         *
         * so the signal is exp(+j*2*pi*SGN*fx*dTau) and recovering dTau takes
         * the opposite-sign kernel: SGN = +1 wants the FORWARD transform. That
         * is what this branch does, and it is correct for conformant data.
         *
         * CAPELLA'S PRODUCTS DO NOT MATCH THEIR OWN DECLARED SGN. They declare
         * SGN = +1, but their signal is exp(-j*2*pi*fx*dTau), which the
         * standard would label SGN = -1. Their own reference notebook
         * (CPHD_by_Example.ipynb, capellaspace/jupyter-notebooks)
         * range-compresses with the INVERSE transform and never reads SGN at
         * all:
         *
         *     signal_time = np.fft.fftshift(np.fft.ifft(signal_chip, axis=1),
         *                                   axes=1)
         *
         * Read to the standard, their data comes back MIRRORED IN RANGE about
         * the reference point -- and mirrored coherently, so it focuses
         * sharply and looks entirely correct. On the Giza collect that put the
         * Great Pyramid roughly a kilometre from its own coordinates, with
         * open desert where it should be. Focusing at the three pyramids'
         * declared positions and measuring the bright structure at each:
         *
         *              true base   as standard        with the override
         *   Khufu        230 m     81 m @ 229 m off   271 m @ 117 m off
         *   Khafre       215 m    267 m @   7 m off   284 m @ 134 m off
         *   Menkaure     105 m    174 m @ 237 m off   130 m @  74 m off
         *
         * Khafre is the trap in all of this: at Y = -88 it sits almost on the
         * mirror axis, so it barely moves and appeared to confirm every wrong
         * hypothesis tried -- including two attempts to blame the image-plane
         * axes, which were correct throughout.
         *
         * So the standard behaviour is kept and the override is confined to
         * the vendor that needs it, keyed on CollectorName. A conformant CPHD
         * is unaffected. */
        if (domain_fx == 1) {
            st = (sgn_eff < 0.0) ? rs_fft_inverse(plan, line) : rs_fft_forward(plan, line);
            if (st != RS_OK) goto fail;
        }

        float complex *dst = cphd->signal + i * n_rbin;
        for (size_t k = 0; k < n_rbin; k++) {
            dst[k] = line[(k + (size_t)n_samp - half) % (size_t)n_samp];
        }
    }

    /* The scene plane, so a focused grid cell can be put on the earth. These
     * are the same vectors the platform positions were rotated into above, so
     * grid coordinates and geolocation cannot drift apart. */
    for (int i = 0; i < 3; i++) {
        cphd->plane.origin[i] = origin[i];
        cphd->plane.u_x[i] = ex[i];
        cphd->plane.u_y[i] = ey[i];
    }
    {
        double h = 0.0;
        (void)rs_geo_ecf_to_llh(origin, NULL, NULL, &h);
        cphd->plane.ref_hae = h;
    }
    cphd->plane.is_slant = 0;   /* the image area reference surface is the ground */
    cphd->plane.valid = 1;

    /* A real collect is demodulated against the reference point as it is
     * received, so the reference delay is already out of the data. */
    cphd->phase_ref_srp = 1;

    cphd->dr = dr;
    cphd->fc = fc;
    cphd->lambda = RS_C / fc;
    cphd->r_near = cphd->r_ref[0] - (double)half * dr;

    /* PRF from the pulse times actually kept, so that a strided read reports
     * the rate the samples were taken at rather than the rate they were
     * transmitted at -- the sub-aperture stage derives its frequency reach from
     * this and would otherwise overstate it. */
    const double span = cphd->t[n_pulse - 1] - cphd->t[0];
    if (!(span > 0.0)) {
        rs_set_error("cphd: pulse times do not increase; cannot derive PRF");
        st = RS_ERR_RANGE;
        goto fail;
    }
    cphd->prf = (double)(n_pulse - 1) / span;
    snprintf(cphd->source, sizeof cphd->source, "cphd:%lux%lu:%.1fs",
             n_vec, n_samp, span);

    st = RS_OK;
    goto done;

fail:
    rs_cphd_free(cphd);
done:
    rs_fft_plan_destroy(plan);
    free(line);
    free(raw);
    free(valid);
    free(pvp);
    free(xml);
    fclose(f);
    return st;

io:
    fclose(f);
    free(xml);
    return RS_ERR_IO;
}
