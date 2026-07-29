/* Minimal PNG encoder, private to src/io.
 *
 * Not in include/resonarsat/ on purpose. This is an implementation detail of
 * raster.c, which chooses a container from the output path; nothing else in the
 * project should be aware that PNG exists, in the same way that nothing outside
 * src/core/fft.c names the vendored FFT.
 *
 * WHY HAND-ROLLED RATHER THAN libpng. The project builds with no external
 * libraries, and that property is worth more than the few hundred bytes this
 * costs. PNG makes it cheap to keep: the format's compressed stream is a zlib
 * wrapper around DEFLATE, and DEFLATE defines a STORED block type carrying raw
 * bytes verbatim. Emitting only stored blocks produces a completely valid PNG
 * that every decoder reads, and needs two checksums rather than an entropy
 * coder.
 *
 * WHAT THAT COSTS. The output is not compressed. A stored stream is the raw
 * data plus one filter byte per row and five bytes per 65535-byte block, so
 * roughly 0.01% overhead rather than the 40-60% a real encoder would save on
 * this kind of imagery. For quicklooks and tomogram sections -- the images this
 * writes -- that is not a consideration. If it ever becomes one, the fix is a
 * fixed-Huffman block type, not a dependency. */

#ifndef RS_IO_PNG_H
#define RS_IO_PNG_H

#include <stddef.h>

#include "resonarsat/resonarsat.h"

/* Write an 8-bit PNG.
 *
 * 'pixels' is row-major with 'channels' bytes per pixel, tightly packed and
 * with no row padding: 1 for greyscale, 3 for RGB. Any other value is refused
 * rather than guessed at.
 *
 * The image is written non-interlaced with every row using filter type 0
 * (None). Filtering exists to help the compressor, and there is no compressor
 * here, so a predictor would cost time and change nothing.
 *
 * Returns RS_ERR_ARG on a NULL argument, zero dimension or unsupported channel
 * count, RS_ERR_ALLOC if the stream buffer cannot be sized, and RS_ERR_IO if
 * the file cannot be opened or written. */
resonarsat_status_t rs_png_write(const char *path, const unsigned char *pixels,
                                 size_t n_col, size_t n_row, int channels);

#endif /* RS_IO_PNG_H */
