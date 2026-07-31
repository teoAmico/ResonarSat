# Architecture

How the software is put together: its layers, the data structures that flow
between them, and the rules that keep the stages honest about what they know.

For what the commands do and how to run them, see the command reference at the
end and `resonarsat <command>` with no arguments.

---

## Shape of the system

A C11 library with a thin command-line program on top. No external
dependencies — the only third-party code is a vendored FFT, and image output is
written directly rather than through a graphics library. That is a deliberate
constraint, not an accident of scope: it means the whole thing builds anywhere
with a compiler and CMake, and there is nothing to install before a result can
be checked.

```
  src/readers/     file formats  ->  in-memory products
  src/core/        every processing stage
  src/util/        geometry, geodesy, orbits, error state
  src/io/          raster output
  src/main.c       argument parsing and stage sequencing only
  tests/           one executable per subsystem
```

`src/main.c` is large but contains no science. It parses options, sequences
stages, prints results and enforces the guards that need to see more than one
stage at once. Every quantity it reports is computed in the library.

---

## The data path

One type per stage. Each is produced by one stage and consumed by the next, and
each carries the metadata needed to interpret it — so no stage has to reach back
past its input to find out what it is holding.

```
   file on disk
        |
        |  rs_read_cphd() / rs_read_sicd() / rs_read_uavsar()
        v
   rs_cphd_t          pulses, platform positions, times, range-compressed
        |             signal. The only form that still contains individual
        |             pulses; everything downstream is derived.
        |
        |  rs_focus_backproject(), over an rs_grid_t
        v
   rs_slc_t           a focused complex image, plus the geometry it was
        |             focused with
        |
        |  rs_subaperture_from_cphd() or rs_subaperture_split()
        v
   rs_subap_stack_t   N images of the same ground at different instants,
        |  \          with the time between them. This is the film.
        |   \
        |    \  rs_ccd_locate()
        |     v
        |    rs_ccd_t        one map: where in the scene something changed
        |                    between consecutive sub-apertures. A branch, not
        |                    a link in the chain — it answers "where", not
        |                    "how much", and nothing downstream consumes it.
        |
        |  rs_microm_track()
        v
   rs_microm_t        per-window displacement and velocity series,
        |             plus a quality value per window
        |
        |  rs_spectrum_compute_opts()
        v
   rs_spectrum_t      per-window power spectra, dominant frequency,
        |             prominence
        |
        |  rs_tomo_focus()
        v
   rs_tomo_t          per-window depth profiles, plus the parameters and
                      derived constants that produced them
        |
        |  rs_raster_write_map() / rs_raster_write_cube()
        v
   PNG, PGM, float32 cube + sidecar
```

`rs_grid_t` sits beside this rather than in it: a rectangular patch of ground in
a scene-local frame, with an origin, cell size and extent. It is what makes
`--offset` possible, and it is passed to focusing rather than being implied by
the product.

### Why the types carry their own metadata

`rs_slc_t` holds wavelength, slant range, platform speed, azimuth spacing and
sampling rate alongside the samples. `rs_tomo_t` holds the full parameter set it
was focused with. The alternative — passing bare arrays and remembering what
they meant — is how a depth axis ends up scaled by a constant nobody can later
identify.

The same reasoning drives a rule in `slc.h`: a reader must never populate the
derived timing fields itself. It fills what the file states, then calls the
shared derivation. That prevents a product's transmit pulse rate being written
onto the azimuth sampling axis, which is a different quantity and produces
plausible, wrong answers everywhere downstream.

---

## Layers

### Readers — `src/readers/`

Turn vendor files into `rs_cphd_t` or `rs_slc_t`. Three formats: CPHD 1.x phase
history, SICD/NITF focused imagery, and UAVSAR SLC with its annotation file.

Readers parse untrusted binary. Every field that sizes an allocation or a file
offset is validated before use, block extents are checked against the actual
file length, and the declared dimensions must account for exactly the declared
blocks — which is what distinguishes a truncated download from a valid file
before an hour of processing is spent on it.

CPHD supports both sample formats in the wild: CF8 (two big-endian float32) and
CI4 (two big-endian int16). Sample width is threaded through the size checks and
the decode rather than assumed, because assuming it fits the reader to whichever
vendor it was written against.

### Core — `src/core/`

| module | responsibility |
|---|---|
| `fft.c` | the only file that names the vendored FFT; everything else uses `rs_fft_*` |
| `focus.c` | time-domain backprojection onto an `rs_grid_t` |
| `subaperture.c` | dividing the aperture, in time or in Doppler |
| `coreg.c` | sub-pixel shift between two complex patches |
| `phaselink.c` | split-band phase linking across a whole stack |
| `microm.c` | per-window tracking; owns the estimator choice |
| `ccd.c` | scale-invariant change detection over a sub-aperture stack |
| `validate.c` | pre-flight checks of a collect against an intended measurement |
| `spectrum.c` | windowed periodograms, dominant frequency, prominence |
| `tomo.c` | the four depth models and their solvers, and the alignment null |
| `simulate.c` | synthetic phase history over a real collect's geometry |

Backprojection is chosen over frequency-domain focusing because it makes no
far-field or planar-wavefront approximation — which matters over the
multi-second dwells this method needs — and because it accepts an arbitrary
subset of pulses without any change of formulation. That property is what lets
the sub-aperture stage form its looks by calling it repeatedly.

The dependency graph is shallow and acyclic. Core modules depend on `fft` and
`geom`; nothing in core depends on `readers`, `raster` or `main`.

`ccd.c` is the one stage that is not a link in the chain. It consumes the
sub-aperture stack and produces a map, and nothing consumes the map — it answers
"where in this scene did something change", not "by how much and at what
frequency". It exists because every measurement in the chain sweeps hundreds of
windows with no prior on where a peak should be, and the number of tries is
itself a source of peaks. It is also the one stage that does not depend on the
tracker succeeding, which matters on scenes where the tracker does not.

Both stages that make a claim now carry a null test built for that claim, and
the two are deliberately different instruments. `rs_null_static()` and
`rs_shuffle_looks()` bound the micro-motion stage; `rs_tomo_alignment_null()`
bounds the depth stage by circularly shifting each window's profile, which holds
every per-window artefact fixed and destroys only the agreement between windows.
Using the wrong one of these is not a minor error — a shuffled null passed a
phase measurement at p = 0.03 that a static null then placed within 1% of a
motionless scene. See `runs/giza/2026-07-30-uniform-phase-khufu/`.

### Utilities — `src/util/`

`geom.c` holds the constants chain — acoustic wavelength, orbital aperture,
depth resolution — computed from explicit parameters every time, with no
published value ever inherited as a constant. `geocode.c` converts between the
scene frame, earth-centred coordinates and geodetic latitude and longitude,
including range-Doppler geolocation for slant-plane grids. `orbit.c`
interpolates state vectors. `status.c` holds the error string.

### I/O — `src/io/`

`raster.c` chooses a container from the output path's extension and quantises
once, so callers never select a format. `png.c` writes PNG using stored DEFLATE
blocks, which needs only two checksums rather than a compression library — the
files are slightly larger and the project keeps its zero-dependency property.

---

## Cross-cutting design rules

**Every entry point returns a status.** `resonarsat_status_t`, with a
human-readable message retrievable separately. Output structures are zeroed
before any validation can fail, so "check the status, then free the struct" is
safe on every path — and the output pointer is checked before that zeroing, so a
null pointer returns an error rather than dereferencing.

**Assumptions have no defaults.** The constants that scale the depth axis must
be supplied explicitly — velocity always, and the investigation frequency for
every model that assumes one rather than measuring it. A depth product that came
from library defaults would look identical to one that came from a considered
choice.

**Guards fail loudly, and name the numbers.** A depth cell finer than the
geometry supports is refused, because a finer grid interpolates rather than
resolving. A depth range beyond the unambiguous limit is refused, because
features fold. A sub-aperture route that cannot form the requested decomposition
refuses rather than substituting a different one.

**Products carry their provenance.** Every depth cube is written with a sidecar
recording the model, solver, assumed constants, derived resolutions, and the
full measurement chain — sub-aperture route, estimator, look count, overlap,
window, coherence threshold, grid offset, band floor, detrend. Outputs are
written one directory per run with a manifest of the verbatim commands.

**Null tests are part of the pipeline, not an afterthought.** Two are built in:
shuffling the sub-look time order, and simulating a motionless scene over the
real collect's geometry and running the identical chain. The second exists
because the first cannot control for heavily overlapped sub-looks, whose
measurements vary smoothly before anything moves.

---

## Testing

One executable per subsystem, run by `ctest`, plus a documentation-coverage
check registered as a test of its own: `scripts/check_docs.py` fails if any
function in `src/` or `tools/` lacks a block comment. It is a `ctest` failure
rather than a build failure, and it is registered only when CMake finds a
Python 3 interpreter — so a green build is not evidence the check ran.

Tests build their own fixtures; none reads anything from `data/`, so the suite
passes on a fresh clone. The synthetic generator produces phase history for
targets whose positions and vibration frequencies are known exactly, which makes
it the only place ground truth exists.

Readers are additionally tested against a corpus of deliberately malformed
files — truncated blocks, impossible dimensions, offsets past the end — where
the required behaviour is a status code and a message, never a crash.

ASan and UBSan are a first-class build configuration rather than an
afterthought, because the readers parse untrusted external formats.

---

## Command reference

Each command is one stage of the data path above. `resonarsat <command>` with no
arguments prints full usage.

| command | consumes | produces |
|---|---|---|
| `feasibility` | acquisition parameters | the observable vibration band and its cost in resolution |
| `info` | any supported file | geometry and timing; validates the file |
| `validate` | `rs_cphd_t` + an intended measurement | whether the collect can support it, check by check, before any processing |
| `focus` | `rs_cphd_t` | a focused image, written as PNG or PGM |
| `mmotion` | `rs_cphd_t` | vibration spectra per window, null-test results, and with `--ccd-out` a change-detection map |
| `tomo` | `rs_cphd_t` | depth profiles, a section image, a float32 cube, and with `--geocode` a CSV of window centres in latitude and longitude |
| `sweep` | `rs_cphd_t` | how recovered depth responds to the assumed constants |

**Options that change what is measured**, rather than how well:

- `--offset X,Y` — moves the processing grid. It defaults to the scene centre,
  and a grid smaller than the collect otherwise covers only the middle of it.
- `--subap pulse|uniform|paper` — how the aperture is divided. `paper` follows
  the published decomposition, which holds out part of the Doppler band;
  `pulse` divides in time; `uniform` into equal bands.
- `--estimator correlation|phase|splitband` — what is tracked. Correlation
  measures patch displacement and has an unambiguous range beyond which it
  folds; phase reads a single scatterer's phase and wraps beyond about 16 mm of
  line-of-sight motion between frames.
- `--model A|B|C|D` — how depth is derived. `A` is the published method, `B` a
  standing-wave heuristic that fails differently and so serves as a cross-check,
  `C` conventional multi-baseline tomography as an uncontested reference, and
  `D` a resonance condition reading depth from the measured frequency rather
  than from a baseline. `--solver dft|lstsq` selects how `A` inverts.
- `--convention paper|patent` — which acoustic-wavelength convention applies,
  `v/(2f)` or `v/f`. The two published sources disagree here, and the choice
  scales every depth by a factor of two, so it is never assumed silently.
- `--y shifts|los` — what the tomographic observable is built from. The paper's
  `Y = dx + i·dy` is defined on coregistrator shifts; `los` builds it from
  phase-derived line-of-sight displacement instead, which is a departure from
  the paper and is recorded as one.
- `--velocity`, `--frequency` — the assumed constants that scale the depth axis.
  No defaults. `--frequency` is required for models A, B and C; model D measures
  it from the spectrum instead, leaving velocity as its only free constant.
- `--patent-exact` — the historical name for the unconditioned patent-chain
  interpretation. It refuses hidden departures from the selected chain: raw
  `Y`, no depth taper, unregularised
  pseudoinverse, `lambda = v/f`, Model A only, shift-based complex `Y`,
  spectral `--subap paper`, rectangular sub-aperture filters, no coherence mask,
  no spectrum/detrend stage, and `--reference pair` — the master-slave pair held
  `B_shift` apart and swept rigidly, which is what WO2024008365A1 [0004] and its
  Figures 0.2 and 0.3 describe. It uses the conventional
  `A[i,j]=exp(j*Kz_i*z_j)` repair because rendered Eq. 22's `2*pi*t` and matrix
  dimensions are inconsistent. Geocoding remains the optional export block.
- `--eq22-literal-t VALUE` — experimental fixed-scalar reading of the rendered
  Eq. 22, available only for Model A with `--solver lstsq`. The papers define
  `t` earlier as acquisition time but do not index it in the matrix, state its
  sampling rule, or include it in the resolution. The option therefore tests
  sensitivity to the printed factor; it is not a uniquely sourced time model.

**Options that check whether a result is real:**

- `--null-trials N` — shuffles the sub-look time order N times.
- `--null-static N` — simulates N motionless scenes over the real geometry and
  runs the identical chain.
- `--coherence F` — masks windows whose sub-looks do not correlate.
- `--fmin HZ` — excludes the lowest frequency bins, where drift lives.

**Order.** `feasibility` before downloading; `info` to confirm the file arrived
intact; `focus` to confirm the grid is over the intended ground; `mmotion` with
a null test; `tomo` only if that null was cleared; `sweep` to test whether the
depth axis responds to the ground or to the operator.
