/* Product readers.
 *
 * Every function here parses an untrusted external file, and every one obeys
 * the same contract: a truncated, corrupt or hostile input produces a status
 * code and a descriptive message, never a crash, a partial write into the
 * caller's struct, or a silently zero-filled image. Fields that determine an
 * allocation or a file offset are validated before they are used.
 *
 * UAVSAR, CPHD 1.x and SICD are implemented. CSK (HDF5) and Sentinel-1 SAFE
 * (GeoTIFF plus annotation XML, and requiring TOPS deburst and azimuth deramp
 * before any sub-aperture work) are the planned additions.
 *
 * CPHD was written before SICD because it is much the easier container -- an
 * ASCII key-value header carrying explicit byte offsets, against NITF's forty
 * unnamed fixed-width fields at implied offsets -- and because focusing shifted
 * pulse windows from phase history is a more faithful sub-aperture operation
 * than filtering an already focused image.
 *
 * SICD earns its place regardless: it is the format the published method
 * actually starts from, so a product read through it and passed to
 * rs_subaperture_split() follows the source's own flowchart, and it is the
 * format most commercial X-band data ships in. */

#ifndef RESONARSAT_READERS_H
#define RESONARSAT_READERS_H

#include "resonarsat/focus.h"
#include "resonarsat/resonarsat.h"
#include "resonarsat/slc.h"

/* Read a UAVSAR single-look complex image and its ASCII annotation.
 *
 * 'slc_path' is the flat binary of interleaved float pairs; 'ann_path' is the
 * accompanying .ann file. On success 'img' owns a sample buffer the caller must
 * release with rs_slc_free(); on any failure 'img' is left with no allocation.
 *
 * The declared dimensions are cross-checked against the actual file size before
 * anything is allocated, and samples are read row by row so that a truncated
 * file fails at the offending line rather than after a multi-gigabyte read.
 *
 * Azimuth timing is derived from the annotated along-track spacing and platform
 * speed, not from any field named "PRF". Where the annotation omits platform
 * speed a nominal value is assumed and recorded in the image's source string,
 * so a derived timing is always traceable to whether it was measured or
 * assumed.
 *
 * Returns RS_ERR_MISSING_META naming the field if a required annotation entry
 * is absent, RS_ERR_FORMAT if the file is truncated or the dimensions are
 * implausible, and RS_ERR_RANGE if the assembled metadata fails
 * rs_slc_validate(). */
resonarsat_status_t rs_read_uavsar(const char *slc_path, const char *ann_path, rs_slc_t *img);

/* What to keep when reading a CPHD, and why keeping everything is not an
 * option.
 *
 * A long-dwell collect is large in the one dimension that does not matter. The
 * reference collect for this reader carries 63291 pulses of 7499 range samples
 * -- 3.8 GB of signal -- across a 7 km swath, while the sub-aperture stage
 * consumes pulses and discards nearly all of that swath. Range compression has
 * to see every sample of a pulse to form its profile at full resolution, but
 * once formed only the bins around the target need to be retained. Compressing
 * pulse by pulse and keeping a window turns 3.8 GB into 259 MB for a 512-bin
 * window, which is the difference between the file fitting in memory and not.
 *
 * 'rbin_window' is the number of range bins to retain, centred on the scene
 * reference point; 0 keeps every bin. 'pulse_stride' keeps every nth pulse and
 * 0 or 1 keeps all of them -- note that striding lowers the effective PRF and
 * so lowers the vibration frequency the sub-aperture stage can reach without
 * aliasing, which is why it defaults off. 'max_pulses' caps the count after
 * striding, for smoke tests; 0 means no cap. */
typedef struct rs_cphd_read_opts {
    size_t rbin_window;
    size_t pulse_stride;
    size_t max_pulses;
} rs_cphd_read_opts_t;

/* Read a CPHD 1.x phase-history file into the focusing container.
 *
 * 'opts' may be NULL, which keeps every pulse and every range bin. On success
 * 'cphd' owns buffers the caller must release with rs_cphd_free(); on any
 * failure 'cphd' is left with no allocation.
 *
 * The file is an ASCII header of "KEY := VALUE" lines terminated by a form
 * feed, giving byte offsets for an XML metadata block, a per-vector parameter
 * (PVP) block and a signal block. Only the handful of XML fields this pipeline
 * needs are extracted, by tag scan rather than a general parser -- a full XML
 * implementation would be a dependency and an attack surface out of proportion
 * to reading fifteen scalars.
 *
 * Two format details drive the implementation. All CPHD binary data are
 * big-endian, so every double and float is byte-swapped on the little-endian
 * machines this runs on. And the signal array is normally in the FX (frequency)
 * domain, meaning it is *not* range compressed: this reader transforms each
 * pulse along the sample axis to produce the range-compressed profiles that
 * rs_cphd_t is documented to hold. Global/SGN selects the transform direction,
 * since it records the sign convention the collector used.
 *
 * Platform positions are converted from earth-centred fixed coordinates into
 * the local planar frame the file defines for its own image area
 * (SceneCoordinates/ReferenceSurface/Planar), with the image area reference
 * point as origin. That is the frame rs_grid_t addresses, so a grid centred on
 * {0,0,0} is centred on the scene.
 *
 * Returns RS_ERR_FORMAT if the header, the XML or any declared block extent is
 * inconsistent with the file's actual size, RS_ERR_UNSUPPORTED for valid CPHD
 * this build cannot yet handle (bistatic collects, sample formats other than
 * CF8, TOA-domain signal arrays), RS_ERR_MISSING_META naming the first required
 * XML element that is absent, and RS_ERR_RANGE if the assembled geometry is not
 * physically sane. Every field that sizes an allocation or a file offset is
 * validated before it is used. */
resonarsat_status_t rs_read_cphd(const char *path, const rs_cphd_read_opts_t *opts,
                                 rs_cphd_t *cphd);

/* Read a SICD product: a focused complex image in a NITF 2.1 container.
 *
 * On success 'img' owns a sample buffer the caller must release with
 * rs_slc_free(); on any failure 'img' is left with no allocation.
 *
 * This is the format the published method starts from -- its flowchart runs SLC
 * image, two-dimensional transform, band-pass filter, inverse transform, pixel
 * tracking -- so a product read here and passed to rs_subaperture_split()
 * follows the described processing directly, where the CPHD route refines it.
 *
 * NITF carries no key names: every field is a fixed width at an implied offset,
 * so a miscount anywhere shifts everything after it. The reader validates by
 * arithmetic rather than by trust, requiring the field walk to land exactly on
 * the declared header length and the pixel dimensions to account exactly for the
 * declared segment size. The SICD XML sits after the pixel data, so the geometry
 * is read by seeking past the image before any sample is loaded.
 *
 * SICD's Row axis is range and its Col axis is azimuth, so the file is ordered
 * [range][azimuth] while rs_slc_t is [azimuth][range]. The reader transposes.
 * Loading it verbatim would place cross-track data on the azimuth axis, which
 * still looks like a scene and is therefore the kind of error that survives
 * inspection.
 *
 * Returns RS_ERR_FORMAT if the container is inconsistent with itself,
 * RS_ERR_UNSUPPORTED for valid NITF this build does not interpret -- multiple
 * image segments, compressed or blocked pixel data, anything but 32-bit IEEE
 * float I/Q -- and RS_ERR_MISSING_META if the SICD XML or its transmit frequency
 * is absent. */
resonarsat_status_t rs_read_sicd(const char *path, rs_slc_t *img);

/* Read a SICD product's metadata, leaving its pixels on disk.
 *
 * Every field rs_read_sicd() sets is set here too -- geometry, timing, carrier,
 * spacings, dimensions -- except that '*img' owns no sample buffer and 'data'
 * stays NULL. Do not pass the result to anything that reads pixels.
 * rs_slc_free() is still the correct way to release it.
 *
 * The distinction is worth having because the two costs differ by orders of
 * magnitude. A spotlight SICD of 40000 x 40000 samples is 12.8 GB as complex
 * float, and the geometry needed to decide whether a collect can support a
 * measurement is a few hundred bytes of XML sitting behind it. 'resonarsat
 * validate' uses this so that screening a collect does not cost more than
 * processing one.
 *
 * Fails exactly where rs_read_sicd() fails, minus the pixel pass. */
resonarsat_status_t rs_read_sicd_meta(const char *path, rs_slc_t *img);

#endif /* RESONARSAT_READERS_H */
