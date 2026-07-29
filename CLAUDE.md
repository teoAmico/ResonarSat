# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build, test, run

No external dependencies. C11 compiler and CMake are enough; the FFT is vendored.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

```sh
ctest --test-dir build -R test_nullmotion --output-on-failure   # one test
./build/test_nullmotion                                        # or run it directly
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=ASAN               # ASAN+UBSAN, a first-class config
python3 scripts/check_docs.py src tools                        # the docs_coverage gate, standalone
resonarsat <command>                                           # usage text, with no arguments
```

`ctest` passes on a fresh clone with `data/` empty — every test builds its own fixture. Only real-data
runs need a download.

Two build facts that bite:

- **OpenMP is optional and silently absent on macOS.** Apple clang ships without it and Homebrew keeps
  `libomp` keg-only. CMake hunts for the Homebrew prefix and then prints either "OpenMP enabled" or a
  loud warning — read that line. Without it backprojection runs on one core and the binary looks
  identical to a threaded one.
- **`-ffast-math` is deliberately absent** (`CMakeLists.txt`). It permits reassociation and flushes
  denormals, which perturbs sub-pixel correlation peaks and interferometric phase — the exact
  quantities this project measures. Do not add it.

Every function in `src/` and `tools/` must carry a block comment; `docs_coverage` is a ctest that
fails otherwise. New functions need one.

## Architecture

A C11 library with a thin CLI on top. `src/main.c` is ~2600 lines and contains **no science** — it
parses options, sequences stages, prints results, and enforces the guards that need to see more than
one stage at once. Every reported quantity is computed in the library. Put physics in `src/core/`.

```
src/readers/   CPHD, SICD, UAVSAR  ->  in-memory products
src/core/      every processing stage
src/util/      geometry, geodesy, orbits, error state
src/io/        raster output (PNG/PGM written directly, no graphics library)
tests/         one executable per subsystem
```

**One type per stage**, each produced by one stage and consumed by the next, each carrying the
metadata needed to interpret it — so no stage reaches back past its input:

```
rs_cphd_t         pulses, platform positions, times, range-compressed signal.
     |            The only form still containing individual pulses.
     |  rs_focus_backproject() over an rs_grid_t
rs_slc_t          focused complex image + the geometry it was focused with
     |  rs_subaperture_from_cphd() / rs_subaperture_split()
rs_subap_stack_t  N images of the same ground at different instants. The film.
     |  rs_microm_track()
rs_microm_t       per-window displacement/velocity series + quality per window
     |  rs_spectrum_compute_opts()
rs_spectrum_t     per-window power spectra, dominant frequency, prominence
     |  rs_tomo_focus()
depth cube
```

`docs/ARCHITECTURE.md` carries the full version. Subcommands: `feasibility`, `info`, `focus`,
`mmotion`, `tomo`, `sweep`.

Nothing outside `src/core/fft.c` may include or name PocketFFT — see `include/resonarsat/fft.h`.

## Running the Giza measurement with `--patent-exact`

### 1. Get the collect (39 GB)

`data/README.md` has the `curl` line. `-C -` is not optional; the peer reset the original transfer at
12 GB. **Verify the byte count is exactly 39,180,270,944** — a truncated CPHD looks entirely ordinary
on disk, and the reader only rejects it when run.

### 2. Point the grid correctly

This is where the project has repeatedly lost days. `--at 29.979175,31.134186` resolves to
`--offset -152,-552` and lands on Khufu. Use coordinates, not a hand-typed offset. Capella declares
`SGN = +1` while shipping the opposite convention; `rs_read_cphd()` compensates for Capella products
only, keyed on `CollectorName`, and says so on stderr.

**A misplaced grid is invisible in the output** — it still produces a complete, well-focused image of
the wrong ground. Focus and look before measuring. `data/README.md` records the confirmed positions of
all three pyramids and how they were established (by structure size, not position — the near-symmetric
line fits two assignments equally well at 85 m rms).

Never use Khafre as the registration check: at `Y = -88` it sits almost on the mirror axis, barely
moves, and confirmed every wrong hypothesis in turn.

### 3. Three settings that each fail by producing a plausible empty image

- **Match the aperture to the cell.** The collect focuses to 0.051 m; sampling that onto a 2 m grid
  aliases into fully developed speckle that reads as "the target is not there". Shorten with
  `--max-pulses` (cell 2.0 m → 8348, cell 1.0 m → 16696). Also 10× faster.
- **`--rbins` is a range window**, N/2 bins either side of the scene reference point at 0.2498 m each.
  Khufu is 552 m out in Y, so it needs `--rbins 4096` (±512 m); `--rbins 2048` cannot contain it at any
  offset. `--rbins 8192` at full PRF needs 22 GB.
- **`--pulse-stride` bounds grid width** and lowers the effective PRF, and therefore the observable
  vibration band. **No measurement run may use it.**

### 4. The command

```sh
resonarsat tomo --cphd data/giza.cphd --at 29.979175,31.134186 \
    --size 256 --grid-cell 2.0 --n 512 --rbins 4096 \
    --velocity 6600 --frequency 22000 \
    --cell 1.3 --depth 60 \
    --patent-exact --out khufu.f32 --geocode khufu.csv
```

`--velocity` and `--frequency` are **required and have no defaults** (except model D, which measures
the frequency). They are *assumed* constants that scale the depth axis; nothing in the pipeline
measures either, and `--frequency` is the paper's "investigation frequency" (12.5 kHz in the paper,
22 kHz in the patent), not a measured vibration frequency.

`--patent-exact` exists so a caller need not know four separate flags to get the published method. It
forces model A, no taper, no mean removal, no regularisation, the true pseudoinverse (Eq. 24),
`--y shifts` (Eq. 21), `--subap paper`, rectangular sub-aperture filters, no coherence mask, the
master-slave pair swept rigidly, **and `lambda = v/f`** — the patent convention. That last one is part
of the flag on purpose: the paper uses `v/2f`, so leaving it out would produce a depth axis twice the
patent's while announcing itself as exact. Incompatible overrides are rejected rather than silently
ignored. Expect it to be noisier; that is the point.

For the micro-motion half alone, the equivalent chain in `mmotion` is
`--subap paper --reference pair --estimator correlation --coherence 0 --no-detrend`.

### 5. Annotate the output

`--out khufu.f32` writes a raw float cube plus `.hdr` (dimensions) and `.meta` (geometry and the
assumed constants) sidecars. `--section-out` writes a section as **bare pixels**, one per window by one
per depth cell — unreadable at native size and carrying no scale. A depth product without an axis
invites being read as a picture of the ground, so render the annotated version:

```sh
python3 tools/plot_tomo_section.py khufu.f32 khufu_tomogram_annotated.png        # middle row
python3 tools/plot_tomo_section.py khufu.f32 khufu_row40.png 40                  # a chosen row
```

It draws the depth and window axes, the colour scale, and the assumed `(v, f)` read straight from the
`.meta` sidecar — which is what keeps a tomogram from circulating without the constants that set its
depth axis. Always publish this rather than the bare `--section-out` pixels.

To annotate the **SAR image** instead — the check that the grid is on the pyramid, per step 2 — dump it
with `focus --raw` and overlay the labelled offset grid:

```sh
python3 tools/annotate_grid.py raw.f32 256 2.0 -152 -552 khufu_grid.png 100 \
    "Khufu:-152:-552,Khafre:-87:-88,Menkaure:77:350"
```

Arguments are `RAW.f32 N CELL OX OY OUT.png [STEP] [TARGETS]`, with target offsets in metres in the
same frame as `--offset`. Gridlines are labelled in the metres `--offset` and `--at` speak, so a
feature can be pointed at and converted back.

Both scripts need `numpy` and `Pillow` — the only place in the project with Python dependencies. The C
build and `ctest` never need them.

**Read positions off the grid; do not infer them from orientation.** Image row is the grid's X and
image column is its Y, and on the Giza collect the product's image-area axes are left-handed with
respect to local up (`(uIAX x uIAY) . up = -0.9999958`), so the layout is mirrored relative to a
north-up map. Compass intuition about which end of a line is north will be wrong.

## Working rules specific to this project

**Null tests are the credibility check, not an afterthought.** `--shuffle-looks SEED`, `--null-trials`,
`--null-static` in `mmotion`; `sweep` for parameter sensitivity. Every micro-motion measurement made
with this implementation so far has failed its null test — see the Status and Honest caveat sections of
`README.md`. Do not present a result as a demonstrated sensitivity.

**Do not read a measurement out of `--reference pair` or `adjacent`.** Neither recovers a frequency on
the synthetic single-target fixture; both return the lowest spectral bin whatever is injected. `pair` is
exposed because it is what the sources describe and should be testable, not because it works.

**Frequencies are trustworthy, amplitudes are not.** Amplitude maps are labelled qualitative in the
metadata and should stay that way.

**Every real-data run gets a directory and a `RUN.md`.** Use `tools/new-run.sh <scene> <suffix>
"<question>"`, which seeds the manifest — including the question the run is meant to answer, recorded
*before* the result. A run that produced a null result keeps its directory: "deleting the runs that did
not work is how a body of evidence quietly becomes a highlight reel" (`images/README.md`).

**`docs/` is a working area and mostly gitignored.** Only an allowlist is tracked (`.gitignore`).
Development notes, review transcripts and plans stay local. The source papers and the patent are
deliberately **not** tracked — `*.pdf` is ignored, and they are cited by DOI and publication number
instead. Keep local copies in `docs/`; do not add them to git.

**`data/` never enters git** beyond its README.
