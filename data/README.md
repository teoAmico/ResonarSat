# data/

Empty on purpose. SAR collects are tens of gigabytes and never enter git; this
file exists so the directory does, and so that a fresh clone records where the
data came from rather than leaving it to be rediscovered.

The test suite needs nothing here — every test builds its own fixture — so
`ctest` passes on a clone with this directory empty. Only the real-data runs
need downloads.

## The Giza collect

The scene the current work is built on. Capella Space open data, CC BY 4.0.

```
CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd
  39,180,270,944 bytes exactly    CPHD 1.1.0, spotlight, HH
  centre 29.97500 N, 31.13068 E   -- 576 m from the Great Pyramid
  dwell 32.87 s   335,149 pulses x 25,073 range bins   PRF 10,196 Hz
  scene about 5 km square: Khufu, Khafre, Menkaure and the Sphinx all inside
```

```sh
curl -L -C - --retry 10 --retry-all-errors -o data/giza.cphd \
  "https://capella-open-data.s3.us-west-2.amazonaws.com/data/2024/10/4/CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012/CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd"
```

`-C -` is not optional. The peer reset this transfer at 12 GB during the
original download, and without continuation a 39 GB fetch restarts from zero on
every reset. **Check the byte count above before using the file.** A truncated
CPHD looks entirely ordinary on disk; the reader rejects it in under a second
with a block-extent error, but only if it is run.

## Where the pyramids actually are in this collect

Established 2026-07-29 by focusing and looking, after three earlier attempts
put grids on empty desert. Read this before choosing an offset.

| target | `--offset X,Y` | confidence |
|---|---|---|
| Khufu | `-300,400` | confirmed by size |
| Khafre | `-100,0` | confirmed by size |
| Menkaure | `100,-400` | confirmed by size |

**Settled by measuring the structures, not by position.** The three sit on a
near-symmetric line, so their positions alone fit two different assignments
equally well (85 m rms either way, one under a negated X and one under a negated
Y). Size breaks the tie, because Menkaure's base is 105 m against 215 and 230 for
the others -- half. Focused at 1 m with the aperture matched to the cell, the
largest connected bright structure at each is:

| position | bright structure | base it matches |
|---|---|---|
| `100,-400` | **132 x 104 m** | Menkaure, 105 m |
| `-100,0` | 267 x 224 m | Khafre, 215 m |
| `-300,400` | 248 x 217 m | Khufu, 230 m |

All three show the chevron of two radar-facing edges; the one at `100,-400`
needs zooming to see, being half the size and sitting among the mastaba field.
Spans exceed the true bases because layover extends a pyramid toward the radar
by `h/tan(theta)` -- 174 m for Khufu, 81 m for Menkaure.

A caution paid for twice: a bright-span statistic over a fixed box gets this
backwards. Measured that way, `100,-400` looks like the LARGEST of the three,
because the box captures the surrounding structures rather than the pyramid.
Isolate the connected structure, or zoom in and look.

### SOLVED: Capella mislabels SGN, and the reader now compensates

The projection disagreed with the imagery all day because the images were
**mirrored in range** about the scene reference point -- coherently, so they
focused sharply and looked entirely correct while placing the Great Pyramid
about a kilometre from its own coordinates.

The cause is a vendor metadata error, and it took three sources to establish:

- **Capella's own notebook** (`CPHD_by_Example.ipynb`) range-compresses with
  `np.fft.ifft` and never reads `SGN` -- so the inverse is what their data
  needs.
- **SARPy**, NGA's reference implementation, documents the standard as
  `Phase(fx) = SGN * fx * dTOA` in cycles, so `SGN = +1` means the FORWARD
  transform. This reader was standard-conformant all along.
- **The imagery** confirms which is right in practice.

So Capella declares `SGN = +1` while shipping `exp(-j2pi fx dTau)`, which the
standard would label `SGN = -1`. `rs_read_cphd()` now inverts the transform for
Capella products only, keyed on `CollectorName`, and says so on stderr when it
does. A conformant CPHD from any other source is unaffected.

Any one source alone would have produced a wrong fix: the notebook alone
suggests inverting globally, which breaks conformant data; the standard alone
says do nothing, which leaves every Capella image mirrored.

**`--at` now works.** `--at 29.979175,31.134186` resolves to `-152,-552` and
lands on Khufu, with the mastaba field's tomb rows visible around it.

Residual offsets of roughly 55-135 m remain between the projected and measured
positions. Layover accounts for the right order and direction -- it displaces a
pyramid's bright edges toward the radar by `h/tan(theta)`, 174 m for Khufu --
but that has not been measured.

**Khafre was the trap.** At `Y = -88` it sits almost on the mirror axis, so it
barely moved and appeared to confirm every wrong hypothesis in turn, including
two attempts to blame the image-plane axes -- which were correct throughout.
The lesson is in `runs/README.md` terms: a target on a symmetry axis cannot
discriminate, and should never be used as the check.

| target | lat, lon | projected `--offset`, as declared |
|---|---|---|
| Great Pyramid (Khufu) | `29.979175,31.134186` | `-152,-552` |
| Khafre | `29.976111,31.130833` | `-87,-88` |
| Menkaure | `29.972500,31.128333` | `+77,+350` |
| Great Sphinx | `29.975278,31.137778` | `+401,-556` |

These are what `--at` produces from the file's declared axes, and they are the
values the imagery disputes -- see above. They come from
this collect's IARP (29.975002 N, 31.130680 E) and are specific to it: another
collect of the same ground has its own axes and its own numbers, which is why
the coordinates are the durable column.

**A misplaced grid is invisible in the output.** It still produces a complete,
well-focused image, just of the wrong ground — there is nothing in a figure to
say it is the wrong place. Register before measuring, and treat any offset that
did not come from the product's own axes or from imagery as a claim to check.

## Three settings that decide whether you see anything at all

All three cost this project a day, and each fails by producing a plausible
empty image rather than an error.

**1. Match the aperture to the cell.** This collect focuses to 0.051 m. Sampling
that onto a 2 m grid aliases into fully developed speckle with no structure —
which reads as "the target is not there". Shorten the aperture with
`--max-pulses` so the resolution matches the cell:

```
aperture for cell d:  L = lambda*R/(2d),  pulses = L/V * PRF
  cell 2.0 m -> --max-pulses  8348      cell 1.0 m -> --max-pulses 16696
```

It is also 10x faster: 26 s against 4 minutes. `focus` prints the cell the
collect supports; if it warns, act on it.

**2. `--rbins` is a range WINDOW**, N/2 bins either side of the scene reference
point, i.e. `+/-N/2 x 0.2498 m` of slant range. A target 550 m out in Y needs
about 470 m of slant window, so `--rbins 4096` (+/-512 m); `--rbins 2048` cannot
contain it at any offset. At full PRF, `--rbins 8192` needs 22 GB — check the
memory before asking for it.

**3. `--pulse-stride` bounds how wide a grid may be.** Azimuth sampling is
unambiguous only over `lambda*R*PRF/(2V)`: 16.7 km at full PRF, 4.2 km at
stride 4, **1.0 km at stride 16**. Beyond it the image fills with aliased energy
indistinguishable from speckle. The reader now prints the limit.

## Other sources

- **Capella open data**, `s3://capella-open-data/data/`, CC BY 4.0, no AWS
  account needed. 413 open CPHD collects. Spotlight dwells run 7–60 s. Note that
  filename timestamps are the *collection* window and not always the dwell, and
  that stripmap collects give under a second on any given target however long
  they run — read `TxTime1`/`TxTime2` from the XML, which a 16 KB range request
  is enough for.
- **Umbra open data**, `s3://umbra-open-data-catalog/sar-data/tasks/`, CC BY 4.0.
  Shorter dwells, 2–8 s, and CF8 rather than Capella's CI4.

## Synthetic fixtures

```sh
./build/sim_cphd     # writes sim.cphd with a known vibrating target
```

Ground truth is written alongside it, which makes this the only data here where
the answer is known in advance.
