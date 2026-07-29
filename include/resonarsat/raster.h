/* Raster output.
 *
 * Dependency-free by design: PNG and PGM are both written directly, and raw
 * float cubes carry a text sidecar. GeoTIFF export needs GDAL and is a later
 * addition; keeping it out means the whole project builds with no external
 * libraries at all, which matters more at this stage than georeferencing does.
 *
 * CONTAINER IS CHOSEN BY THE OUTPUT PATH. A path ending in '.png' writes PNG,
 * anything else writes PGM. Callers pass a filename and do not select a format,
 * which is what lets the command line accept either without a flag that could
 * disagree with the extension the user typed. */

#ifndef RESONARSAT_RASTER_H
#define RESONARSAT_RASTER_H

#include <complex.h>
#include <stddef.h>

#include "resonarsat/resonarsat.h"
#include "resonarsat/slc.h"

/* How a real-valued map is coloured.
 *
 * RS_PALETTE_GRAY is the honest default for anything whose absolute level
 * matters, because a reader can judge relative brightness without a key.
 *
 * RS_PALETTE_VIRIDIS exists for maps that are read as structure rather than as
 * level -- depth sections above all, where the question is where a feature sits
 * rather than how bright it is. Greyscale hides small differences in the middle
 * of the range, which is exactly where a weak depth feature would appear.
 * Viridis is perceptually uniform and keeps its ordering in greyscale print and
 * for the common forms of colour blindness, so it does not manufacture contrast
 * that the data does not contain.
 *
 * RS_PALETTE_ENERGY is the conventional energy ramp: blue for low, green for
 * intermediate, red for high. It is what most readers of a tomogram expect, and
 * being able to say "the red band sits at nine metres" without consulting a key
 * is worth a great deal when the figure is being discussed rather than measured.
 *
 * The implementation is Turbo rather than the older jet. Both run blue to red,
 * but jet's luminance rises and falls several times across the ramp, so it
 * creates bright bands at cyan and yellow that read as edges wherever they land
 * -- structure the data does not contain. Turbo was constructed to keep
 * luminance monotonic while preserving the blue-low, red-high ordering, which
 * removes the false banding without giving up the intuition. It is still a
 * rainbow: it is not safe for red-green colour blindness, and it is the wrong
 * choice for a figure that must survive greyscale reproduction. Prefer viridis
 * where the reader will measure from the image rather than read it.
 *
 * RS_PALETTE_JET is the MATLAB default that the SAR Doppler tomography
 * literature is drawn in, including the tomograms of Biondi & Malanga (2022),
 * whose colour bars run dark blue at the low end through cyan, green and yellow
 * to dark red at the high end over a NORMALISED amplitude axis. It exists here
 * for one purpose: putting our output beside a published figure without the
 * colours themselves being a difference. For that comparison it is the correct
 * choice and RS_PALETTE_ENERGY is the wrong one, because a reader matching two
 * images cannot tell a genuine disagreement from a change of ramp.
 *
 * It is otherwise the worst of the three. Jet's luminance is not monotonic, so
 * it invents bright bands at cyan and yellow that read as boundaries wherever
 * they fall, and a normalised colour bar makes every figure's red mean
 * something different. Use it to compare, not to conclude. */
typedef enum {
    RS_PALETTE_GRAY = 0,
    RS_PALETTE_VIRIDIS = 1,
    RS_PALETTE_ENERGY = 2,
    RS_PALETTE_JET = 3
} rs_palette_t;

/* Resolve a palette name to its enum, for a command-line option.
 *
 * Accepts "gray"/"grey", "viridis", "energy" and "jet". An unrecognised name
 * returns 'fallback' after a warning on stderr rather than failing the run: a
 * misspelt colour must not throw away an expensive computation. */
rs_palette_t rs_palette_from_name(const char *name, rs_palette_t fallback);

/* Write an 8-bit greyscale amplitude quicklook.
 *
 * Amplitude is displayed on a decibel scale clipped to 'dyn_range_db' below the
 * 99th percentile, which is the convention that makes SAR imagery legible: a
 * linear stretch is dominated by a handful of bright scatterers and shows
 * nothing else. Percentile rather than maximum for the same reason.
 *
 * Format follows the extension of 'path'; see the note at the top of this file.
 * Returns RS_ERR_IO if the file cannot be written. */
resonarsat_status_t rs_raster_write_quicklook(const rs_slc_t *img, const char *path,
                                              double dyn_range_db);

/* Write a real-valued map with a linear stretch between the given limits.
 *
 * Values outside the limits clamp rather than wrap, so a single outlier cannot
 * alias into the middle of the range and read as signal. Pass lo == hi to
 * autoscale to the data's own extremes.
 *
 * Used for the dominant-frequency and quality maps and for tomographic depth
 * sections, where a decibel stretch would be meaningless. 'palette' is ignored
 * for PGM, which has no way to carry colour: asking for viridis while writing a
 * .pgm silently yields the greyscale the container supports rather than
 * failing, because the alternative is refusing to write an image over a
 * cosmetic request. */
resonarsat_status_t rs_raster_write_map(const double *map, size_t n_row, size_t n_col,
                                        const char *path, double lo, double hi,
                                        rs_palette_t palette);

/* Write a raw float32 cube with a text sidecar describing its shape.
 *
 * The sidecar carries the dimensions and axis descriptions, so the binary is
 * interpretable without reference to the code that wrote it. Written for the
 * tomographic profiles, which are three-dimensional and have no natural image
 * representation. 'axis_desc' is copied verbatim into the sidecar. */
resonarsat_status_t rs_raster_write_cube(const double *data, size_t n_plane, size_t n_row,
                                         size_t n_col, const char *path,
                                         const char *axis_desc);

#endif /* RESONARSAT_RASTER_H */
