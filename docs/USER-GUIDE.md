# User guide

A walkthrough from a fresh clone to a depth profile you can defend.

---

## Build

Needs a C11 compiler and CMake 3.20. Nothing else.

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build          # optional; needs no data
```

The binary is `build/resonarsat`. Run it with no arguments for the command list,
or any command with no arguments for its full usage.

**Install OpenMP before your first real run.** Without it every parallel region
compiles away and focusing uses a single core; with it, all of them. On this
class of machine that is roughly a sevenfold difference, and the build gives no
other sign which one you have.

```sh
brew install libomp          # macOS; Linux compilers usually ship with it
cmake -S . -B build          # reconfigure: it prints which outcome you got
```

Configure prints either `OpenMP enabled` or a warning that the build will be
single-threaded. Check which.

---

## Get data

You need a radar product that still contains individual pulses. That means CPHD.
A focused image (SICD, SLC, GeoTIFF) has already integrated the whole
observation, and the sub-aperture stage cannot recover what it needs from one.

Two free sources, both CC BY 4.0 and neither needing an account:

```sh
aws s3 ls --no-sign-request s3://capella-open-data/data/
aws s3 ls --no-sign-request s3://umbra-open-data-catalog/sar-data/tasks/
```

**Choose a spotlight collect.** In spotlight mode the antenna stares at one
place, so every pulse contributes to your target. In stripmap the antenna is
fixed and the beam sweeps past, so a given point is illuminated for well under a
second no matter how long the collection lasts. Stripmap files look attractive
when sorted by duration and are unusable for this.

**Check the real dwell before downloading.** Filenames often carry a collection
window that is longer than the actual pulse span. Read the first 16 KB and look
at the XML:

```sh
curl -s -r 0-16383 "$URL" | strings | grep -E 'TxTime[12]|ModeType'
```

Longer is better. The dwell sets how many frames you can cut and therefore the
highest vibration frequency you can see.

Files run from about 1 GB to 60 GB. For large ones, resume is essential:

```sh
curl -L -C - --retry 10 --retry-all-errors -o data/scene.cphd "$URL"
```

Then confirm the byte count matches the listing. A truncated CPHD looks entirely
ordinary on disk.

---

## Step 1 — is the measurement possible?

Do this before downloading anything large.

```sh
resonarsat feasibility --dwell 20 --range 650000 --velocity 7500 \
    --wavelength 0.031
```

You get a table: for each number of sub-apertures, the highest vibration
frequency observable and the azimuth resolution it costs. Pick the row where the
frequency you care about is comfortably below the limit and the resolution is
still usable.

The trade is unavoidable. Reaching higher frequencies means cutting the
observation into more, shorter pieces, and shorter pieces make blurrier images.
Overlapping the pieces buys back a great deal of this and is what `--overlap`
does later.

---

## Step 2 — check the file

```sh
resonarsat info --cphd data/scene.cphd
```

Confirms the file is intact and prints pulses, range bins, carrier frequency,
pulse rate, dwell and near range. If it reports that a declared block extends
past the end of the file, the download is incomplete — resume it rather than
proceeding.

Note the dwell and pulse count. You need them to choose sub-aperture settings.

---

## Step 2½ — ask whether the measurement is possible at all

```sh
resonarsat validate --cphd data/scene.cphd --frequency 2.0 --amplitude 1.0 \
    --alpha 0.05 --overlap 0.88 --cell 1.0 --size 512
```

`feasibility` works from parameters you type; this works from the collect, and
from the measurement you actually intend. It checks the frequency against the
dwell's resolution, against the sub-aperture Nyquist, against the observation
ratio, and against the displacement-averaging nulls; it checks the aperture
fraction against the range anyone has validated, the grid width against the
spectral route's requirement, the PRF stability against the uniform-timing
assumption, and the amplitude against the tracker's quantisation floor.

**Run it before every measurement.** It costs a metadata read; the thing it
prevents costs twenty minutes and produces a complete, plausible, wrong answer.
The configuration used in `runs/giza/2026-07-30-validated-spot-khufu` fails it
on one line:

```
FAIL  observation ratio   t_sap 1.6434 s at 2.000 Hz gives eta 3.287. The
      sub-look RESOLVES the paired echoes, the tracked feature fragments, and
      correlation tracking fails at any amplitude. At this aperture fraction
      eta stays under 0.20 only below 0.122 Hz.
```

That run was made before this command existed, took twenty minutes, and its
conclusion had to be withdrawn.

**Every check is necessary and none is sufficient.** The last line always reads
`UNKNOWN ground truth`, because whether anything in the scene moves is not a
property of the file. A clean verdict means the configuration is capable, not
that the measurement is real.

---

## Step 3 — find your target

The processing grid is much smaller than the area the satellite covered, and by
default it sits at the scene centre. If your target is elsewhere you will get a
clean image of the wrong ground with nothing to indicate it.

Start with a coarse image of the whole scene to see what is there:

```sh
resonarsat focus --cphd data/scene.cphd --size 1024 --cell 4.0 \
    --out overview.png --dyn-range 45
```

`--size` is the grid in cells and `--cell` its spacing in metres, so this covers
4 km. Water appears dark, because a flat surface reflects away from the
satellite. Buildings and machinery appear bright, often with cross-shaped
sidelobes.

To point the grid somewhere specific, `--offset X,Y` moves it in metres along
and across the satellite's track. If you know your target's latitude and
longitude, `resonarsat info` prints the scene corners in the same coordinates,
which lets you work out the offset.

Then focus finely on the target to check you have it:

```sh
resonarsat focus --cphd data/scene.cphd --size 512 --cell 0.5 \
    --offset 1700,-1700 --out target.png
```

**If a fine image looks like noise, check the stretch before concluding
anything.** The default display clips 40 dB below the 99th percentile, so a few
very bright objects can drive their surroundings to black. Try `--dyn-range 60`,
or `--raw out.f32` to write the amplitudes themselves.

---

## Step 3½ — narrow down *where* to measure

Optional, and worth it on a scene where you do not already know which object you
care about.

```sh
resonarsat mmotion --cphd data/scene.cphd --offset 1700,-1700 \
    --size 512 --cell 1.0 --n 159 --overlap 0.88 \
    --ccd-out scene
```

This runs a change detector over consecutive frames and writes a map. It answers
"where in this scene did something change", not "at what frequency" — a
different question from everything else here, and one that **does not need the
tracker to succeed**, because it never tracks anything.

Its useful property is that it ignores brightness. A bright but stationary
target scores the same as a dim one: only a change in the *structure* of the
signal counts, not its level. That is the false alarm an amplitude-based check
cannot suppress.

Read it with two numbers in mind:

- **1.0 is the no-change value.** The map's minimum tells you more than its
  maximum. On a synthetic scene with one vibrating target, the background sits at
  1.0 and only the target rises. On real distributed clutter it will not: at
  Giza the *minimum over the entire map* was 7.35, because sub-look speckle
  decorrelates everywhere and decorrelation is a change.
- **A map has no threshold.** The method it comes from implements none. Compare
  against a motionless scene through the same chain before reading structure into
  a picture, and check where the extremes actually sit — at Giza the top 0.1% of
  pixels formed a five-column stripe spanning the image, which is a processing
  artefact and not a target.

`--ccd-win` sets the sliding window (default 5) and `--ccd-loading` the noise
floor. The map is written as a PNG and a float32 cube with a sidecar.

---

## Step 4 — measure the motion

```sh
resonarsat mmotion --cphd data/scene.cphd --offset 1700,-1700 \
    --size 256 --cell 2.0 --n 512 --overlap 0.99 \
    --subap paper --estimator phase --win 8 \
    --null-trials 24 --out run/scene
```

Reading the output:

- **observable band** — the highest frequency these settings can represent. If
  what you are looking for is above it, raise `--n`.
- **tracked N windows; M pass the coherence mask** — how many points in the
  scene held together well enough to measure. If M is zero, see below.
- **strongest peak** — frequency, prominence, and a quality value.
- **sub-aperture response** — how much the sub-aperture length attenuates a
  signal at that frequency. Low values mean the amplitude is understated; the
  frequency is unaffected.
- **null floor** — what the same processing produces with the time order
  destroyed.

### Choosing `--n` and `--overlap`

`--n` is the number of frames and `--overlap` how much consecutive frames share.
Together they set both the frequency reach and the sharpness of each frame.

High overlap raises the sampling rate without shortening each frame. It does
**not** raise the frequency reach: each frame averages the motion over its own
length, so the band stops at `1/(2 × frame length)` however finely the frames are
spaced. Overlap buys a finer look at the same band, not a wider one — quoting the
step's Nyquist instead overstates the reach by `1/(1 − overlap)`, which is a
factor of 100 at 99% overlap. To reach higher, shorten the frame. And **overlap
is the second thing to choose, not the first.**

**Set the frame length first.** The quantity that decides whether anything tracks
at all is the *aperture fraction*: the length of one frame divided by the whole
acquisition. Published work validated against ground truth operates at **3 to 8
percent**. Below that the target stops being a distinct feature in each frame and
the tracking collapses — quietly, with a fully populated result.

That is measured, not theoretical. On a 32.9 s Giza collect, frames of 0.30 s
(0.91%) gave a maximum tracking quality of 0.108 across 961 windows, with nothing
above the 0.4 mask. Frames of 1.64 s (5.0%) roughly doubled the quality and still
nothing reached 0.4 — see `runs/giza/2026-07-30-validated-spot-khufu/`.

**Then buy the sampling rate with overlap.** Frame length and sampling rate pull
against each other, and overlap is what separates them. On a long acquisition
this forces high overlap: holding 5% on a 32.9 s collect means 1.64 s frames, and
sampling near 5 Hz then needs an overlap around 0.88. On the 6–15 s acquisitions
the validated work used, 23–49% was enough for the same aperture fraction.

**The sampling rule.** The interval between frame centres must be well under the
period of the vibration you are looking for — a quarter of it is the value the
published work uses. Sampling more slowly than the target vibrates produces a
complete, plausible, entirely wrong result.

If you want a specific frame length in seconds, work back from the pulse rate:
frames of `T` seconds need `T x PRF` pulses each, and a step of `s` seconds
between them means an overlap of `1 - s/T`.

### Choosing `--estimator`

- `correlation` measures how far a patch appears to slide between frames. Robust
  for large, slow motion. It has an unambiguous range: motion beyond roughly
  half the window folds and reports a clean multiple of the true frequency.
- `phase` reads one bright scatterer's phase directly. Far more sensitive, and
  the right choice for small fast motion. It wraps if the target moves more than
  about 16 mm along the line of sight *between consecutive frames*, which heavy
  overlap makes unlikely.
- `splitband` uses all frame pairs at once and needs good coherence throughout.

### If no windows pass the coherence mask

In order of likelihood:

1. **Your tracking window is smaller than one resolution cell.** Each frame is
   blurrier than the full image by roughly the frame count. If `--win` times
   `--cell` is below that, nothing can correlate. Raise `--win`, coarsen
   `--cell`, or lower `--n`.
2. **The grid is over featureless ground.** Check the focused image first.
3. **The threshold is wrong for your scene.** `--coherence 0` disables it, which
   is the right way to inspect what is there before deciding.

---

## Step 5 — test whether the result is real

A pattern in noise looks like a pattern. Two checks are built in.

**Shuffle the frame order** with `--null-trials N`. The intent is that everything
is preserved except the sequence in time, which no real vibration survives. The
output reports how many shuffles matched or beat your measurement.

> **The shuffle does not preserve everything except time, and with
> `--estimator phase` it is not a valid test.** Reordering the frames puts
> non-consecutive looks next to each other, which is exactly when a phase series
> steps furthest. Measured over Giza: the median largest step between frames is
> 0.052 rad in the real order and 1.878 rad shuffled — the shuffle inflates the
> per-step noise 36-fold. That leaves a smooth series being compared against a
> deliberately roughened one, so a slow drift beats its own shuffles easily and
> the test reports a strong detection where there is none. It passed one at
> p = 0.03 that `--null-static` then placed within 1% of a motionless scene
> (`runs/giza/2026-07-30-uniform-phase-khufu/`). Use the shuffle for
> `correlation`, where the observable is a patch position and reordering does not
> change the per-frame noise. For `phase`, treat a passed shuffle as necessary
> and nowhere near sufficient, and go to the static test.

**Simulate a motionless scene** with `--null-static N`. This builds synthetic
data of completely still objects using your satellite's actual path and timing,
and runs the identical processing. Whatever it produces is what this pipeline
yields from a world where nothing moves.

```sh
resonarsat mmotion --cphd data/scene.cphd --offset 1700,-1700 \
    --size 256 --cell 2.0 --n 512 --overlap 0.99 \
    --subap paper --estimator phase --coherence 0 --null-static 8
```

The static test is slower — each trial refocuses — but it is the one to trust
when the two disagree, because it reproduces the processing's own artefacts as
well as its noise. Use it whenever frames overlap heavily.

**If your measurement does not beat the static floor, stop here.** The depth
stage has nothing to invert, and a tomogram built on it will look no different
from one built on a real signal.

Two further checks worth running on any peak you believe:

- Re-run with `--fmin` set above the lowest few bins. A peak that vanishes was
  drift, not vibration — but watch for the peak that *moves* instead. A drifting
  series has a red spectrum, so raising the floor relocates the peak to just
  above the new cut rather than removing it. Run a ladder (`--fmin` 0.3, 1.0,
  2.0) and look at where the peak lands each time: a mode stays put, drift
  follows the cut. Watch the prominence too, which decays as the cut moves away
  from DC when the energy is really at DC.
- Check the reported observation ratio. If each frame spans many cycles of the
  frequency you found, the amplitude is heavily understated.

### `--no-optimize` — checking the arithmetic rather than the result

The null tests ask whether your *result* is real. `--no-optimize` asks a narrower
question: whether the fast code paths gave you the same answer the obvious ones
would. Accepted by `focus`, `mmotion`, `tomo` and `sweep`.

```sh
resonarsat mmotion --cphd data/scene.cphd --offset 1700,-1700 \
    --size 256 --cell 2.0 --n 512 --win 32 --no-optimize
```

It changes two things, worth very different amounts:

- **The correlation peak search becomes exhaustive.** Instead of upsampling the
  neighbourhood of the strongest integer sample, the whole cross-correlation
  surface is zero-padded to the full 1/upsample lattice and searched for a global
  maximum. **This can change a reported shift.** It is the reason the flag exists.
- **Backprojection and the window loop become serial, and change nothing at all.**
  The focused samples are *bitwise identical* — the parallel loop is over grid
  cells and each cell accumulates privately over pulses in chronological order,
  sharing no accumulator, so there is no threading drift to remove.
  `tests/test_focus.c` asserts this to the bit. Treat this half as a
  reproducibility check, not a correction.

**Cost is affordable — single digits, not orders of magnitude.** Measured: the
exhaustive correlator is 1.7–3.2x the optimised one per call (2.5x at the default
32×32 window and 10×20 upsampling), and losing the threads costs about 4x more on
eight cores. You can run this on a full-scale scene; you do not need to crop.

**What to expect from the comparison.** On the synthetic two-target fixture, 10 of
108 tracked samples differ — all of them in windows of tracking quality 0.39–0.61,
with every window above that agreeing on every look. That is the right shape. The
fast path's integer peak is *already* a global maximum over the sampled surface,
so a global search can only beat it where the surface has competing lobes of
comparable height, which is what a window with no target looks like. **Where
neither answer is a measurement, the difference between them is not a
correction** — do not read a `--no-optimize` shift as the more accurate one.

**One trap before you count differences.** The correlation surface is periodic in
the patch size and the two modes fold at slightly different points, so a genuine
agreement can appear as a difference of exactly one patch width — on a 24-pixel
patch, +12.5 and −11.5 are the same place. Fold onto the period first, or your
worst "disagreement" will be an artefact of the wrap.

Products record which mode made them: the `.meta` sidecar carries an
`arithmetic_mode` line and the `.hdr` axis description an `[UNOPTIMIZED]` tag, so
two cubes that differ for this reason cannot be mistaken for two cubes that differ
because the ground did.

---

## Step 6 — depth

```sh
resonarsat tomo --cphd data/scene.cphd --offset 1700,-1700 \
    --size 256 --grid-cell 2.0 --n 512 --overlap 0.99 \
    --subap paper --estimator phase --y los \
    --velocity 3000 --frequency 200 --model A \
    --cell 3.2 --depth 70.4 --null-align 200 \
    --section 24 --section-out tomogram.png --out cube.f32
```

`--null-align` is in that line deliberately: it is the depth stage's null test,
it costs nothing next to the focusing, and a tomogram exported without it has
nothing to say about whether its windows agree. See below.

`--velocity` and `--frequency` have no defaults and must be supplied. They are
assumptions about the ground, not measurements, and together they set the entire
depth scale. Choose them from what you know about the site's material — seismic
velocity in rock runs roughly 1500 to 6000 m/s — and record what you chose. The
one exception is `--model D`, which measures the frequency from the spectrum
rather than taking it from you, leaving velocity as its only free constant.

`--cell` is the depth cell and `--depth` the range. Both are checked: a cell
finer than the geometry supports is refused, because a finer grid interpolates
rather than resolving, and a range beyond the unambiguous limit is refused
because features there fold back and appear shallow. The error messages give the
supported values.

Outputs are a float32 cube with a header, a `.meta` sidecar recording every
parameter, and a section image.

### Choosing `--model`

Four are implemented, and running more than one on the same data is the point:
they fail differently, so agreement between them means something.

- **`A`** — the published method. Steering-matrix inversion over the sub-aperture
  index of a single pass, Biondi & Malanga Eqs. 21-24. This is the contested one.
- **`B`** — a standing-wave heuristic: map each spectral peak to a wavelength and
  accumulate energy at its harmonics. Simpler, and wrong in different ways from
  `A`, which is what makes it worth running beside it.
- **`C`** — classic multi-baseline tomography over genuine perpendicular
  baselines and the radar wavelength. Established physics, no disputed step.
  **Run this whenever the collect allows it.** It is the reference the other
  three are judged against, and a feature that appears only under `A` and never
  under `C` is telling you something.
- **`D`** — resonance mapping, `z = v / (2f)`, reading depth from the measured
  vibration frequency rather than from a baseline. Depth is inverse in
  frequency here, not linear in a bin index, which is why it reaches hundreds of
  metres where `A` reaches tens. It is a materially different reading of the
  source, implemented so the disagreement can be settled by measurement.

`--solver dft|lstsq` selects how `A` performs its inversion; the other models
ignore it.

### The factor of two

The paper and the patent do not agree on how an acoustic wavelength is formed:
the paper uses `lambda = v/(2f)` and the patent `lambda = v/f`. Every depth you
compute is scaled by two depending on which one applies.

```sh
resonarsat tomo ... --convention paper      # v/(2f), the default
resonarsat tomo ... --convention patent     # v/f
```

The software will not pick for you silently, and the choice is written into the
sidecar. `--patent-exact` implies `--convention patent`. Any depth quoted
without stating the convention is ambiguous by a factor of two.

### Patent-chain tomography

Use the historically named `--patent-exact` when you want the unconditioned
patent-chain interpretation, not the conditioned/default processing path:

```sh
resonarsat tomo --cphd data/scene.cphd --offset 1700,-1700 \
    --size 256 --grid-cell 2.0 --n 512 \
    --velocity 6600 --frequency 22000 \
    --cell 1.3 --depth 60 \
    --patent-exact --out patent-eq21-24.f32
```

This selects the following interpreted chain:

```text
tracked complex shifts -> raw Y -> h(z) = A(Kz,z)^dagger Y
```

The flag deliberately refuses hidden departures from that path. It uses Model A
only, complex coregistrator shifts for `Y`, the patent wavelength convention
`lambda = v/f`, spectral `--subap paper`, rectangular sub-aperture filters,
`--reference pair`, no coherence mask, no spectrum/periodogram stage, no
detrend, no mean subtraction, no depth taper, and the unregularised
pseudoinverse.

By default it does **not** implement rendered Eq. 22 literally. The patent PDF prints
`exp(j 2*pi*Kz*t*z)` while defining `Kz` with `4*pi` already. Earlier sections
of the related papers call `t` the acquisition-time variable, but Eq. 22 does
not say whether it means a scalar, `t_i`, a sample interval, or normalized
time, and its resolution formula contains no `t`. The source also declares `A`
as `k x F` while drawing it as `F x k`, and declares `h` as a row where Eq. 23
requires a column. The implementation makes the dimensionally consistent
repair `exp(j*Kz*z)`.

For controlled sensitivity testing only, a fixed-scalar reading of the printed
exponent can be enabled with an explicit value:

```sh
resonarsat tomo ... --model A --solver lstsq --eq22-literal-t 1.0
```

This experimental option multiplies every steering wavenumber by `2*pi*t`,
and consequently scales the reported resolution and unambiguous depth by its
reciprocal. It is rejected for every model except A and for the DFT shortcut.
The value, effective scale, and experimental status are written to metadata.
The software never infers the scalar from acquisition time, look spacing, or
the depth grid; the sources do not provide a reproducible mapping.
The steering matrix remains the repaired, dimensionally consistent `k x F`
orientation. The flag name is retained
for CLI compatibility; metadata labels it `PATENT-CHAIN INTERPRETATION`, not
“equations as written.”

`--reference pair` is part of the flag because the front end is part of the
method. Blocks 3, 4 and 7 of Figure 0.5 are two band-pass filters feeding one
pixel-tracking stage, and [0004] holds that pair "rigidly at a distance
`B_shift`" while sweeping it across the band, so each sample of `Y` is the
offset of one look's slave from *its own* master. Tracking everything against a
single fixed look computes a different observable, and no exactness in
Eqs. 22-24 downstream repairs that. The sidecar records which was used as
`tracking_reference` and withholds the patent-chain label unless it was the
pair.

**That pairing does not recover a frequency on this project's synthetic
fixture**, and the run prints a warning saying so. It returns the lowest
spectral bin whatever is injected. `--patent-exact` therefore preserves the
source's stated chain choices but is *not* a literal or validated
implementation of inconsistent Eq. 22. It is the right flag for testing this
documented patent interpretation, and the wrong one for producing a
measurement. Use the default path for that, and state which you used.

That also means these combinations are rejected: `--model B`, `--model C`,
`--model D`, `--solver dft`, `--y los`, `--estimator phase`,
`--estimator splitband`, `--subap pulse`, `--subap uniform`,
`--reference first`, `--reference adjacent`, and nonzero `--coherence`.

`--patent-exact` is not the whole Figure 0.5 export chain. The patent's block 11
is geocoding, and this code implements it, but keeps it as an export choice. Add
`--geocode FILE.csv` when you want the full Figure 0.5 product form:

```sh
resonarsat tomo --cphd data/scene.cphd --offset 1700,-1700 \
    --size 256 --grid-cell 2.0 --n 512 \
    --velocity 6600 --frequency 22000 \
    --cell 1.3 --depth 60 \
    --patent-exact --out patent-full.f32 --geocode patent-full.csv
```

The geocode export requires the input product to carry a valid scene plane. If
that metadata is absent, exact tomography can still run, but block 11 cannot be
written.

### The depth stage's own null test

The nulls in Step 5 test the *motion* measurement. This one tests the depth
claim, and it is a different instrument because the depth claim is a different
claim.

```sh
resonarsat tomo ... --null-align 200
```

A tomogram is a stack of per-window depth profiles. Every window's profile has
peaks in it whatever the ground does — that is what a steering matrix and a
window function produce. The only thing that distinguishes a structure at depth
from each window's own artefacts is whether the windows **agree** about which
depth.

So `--null-align` circularly shifts each window's profile by an independent
random amount and restacks. Every profile is preserved exactly — same peaks,
same widths, same sidelobes — and the only thing destroyed is the agreement. If
the stacked contrast survives, the windows independently picked the same depth.
If it collapses to the null, the tomogram was showing you the average shape of
per-window artefacts.

**Do not substitute a sub-look shuffle here.** Shuffling rebuilds every profile
as well as their alignment, so it tests the whole chain at once and cannot say
which half carried the result — and it is measurably too weak: an independent
reproduction of this method reports a wrong-shallow-depth case with high contrast
that the shuffle does not catch.

Clearing this null is **necessary and not sufficient**. The depth axis is still
scaled by the assumed velocity and frequency, which nothing measures. Run the
next check too.

### The one check that costs nothing

Run it twice, changing only the assumed velocity, with the depth grid scaled to
match:

```sh
resonarsat tomo ... --velocity 3000 --cell 3.2 --depth 70.4 --out a.f32
resonarsat tomo ... --velocity 6000 --cell 6.4 --depth 140.8 --out b.f32
cmp a.f32 b.f32
```

If the cubes are identical, the depth axis is your assumption rescaled and the
data did not contribute to it. Do this before quoting any depth in metres.

`sweep` automates the same idea across many values and several grid positions:

```sh
resonarsat sweep --cphd data/scene.cphd --velocity 3000 --frequency 200 \
    --range 0.5,2.0 --steps 5 --offsets 0,0:400,-300:-350,250
```

A slope near 1 means depth tracks the assumption exactly.

---

## Reading the images

Depth sections use a colour scale normalised to each image's own maximum, so red
marks the largest value *in that picture*, not a fixed level. Red therefore
appears in every image, including one containing only noise.

This means a tomogram cannot be judged by eye. Read it beside the null-test
result and the velocity check, not after them.

`--palette` selects the ramp: `jet` matches the convention used in the
literature, `viridis` is the choice when someone will measure from the image,
`energy` is a blue-to-red ramp with even brightness, `gray` for print.

---

## Keeping runs straight

```sh
tools/new-run.sh myscene paper-phase "does anything appear below the west wall?"
```

Creates `runs/myscene/<date>-paper-phase/` with a manifest to fill in. Write
outputs there and paste the commands into it. Sidecars record the processing
chain automatically, but only the manifest records why you ran it.

Keep the runs that found nothing. They are most of the evidence.

---

## Troubleshooting

**"a declared block extends past the end"** — the file is truncated. Resume the
download; do not process it.

**"sample format is neither CF8 nor CI4"** — an unsupported variant. Those two
cover Umbra and Capella respectively; anything else is refused rather than
guessed at, because guessing a sample width produces plausible noise.

**"depth cell X is finer than the resolution Y"** — use the value the message
gives, or increase the aperture used by processing more of the collect.

**"depth extent exceeds the unambiguous range"** — reduce `--depth` or raise
`--n`. Beyond that limit, deep features fold and appear shallow.

**"--estimator phase with --y shifts"** — those measure different things and
would produce an empty result. Use `--y los` with the phase estimator.

**"--patent-exact requires ..."** — patent-chain mode is intentionally strict.
Remove the incompatible option, or drop `--patent-exact` if you want the
conditioned processing path. `--no-optimize` is *not* incompatible with it: it
changes how the correlation peak is found, not which model, observable, solver or
wavelength convention is used, so the sidecar still certifies the chain.

**"exhaustive search over ... needs a ... surface, above the ... ceiling"** —
`--no-optimize` builds a surface of `win × upsample` on each axis, and your window
or `--upsample` makes it too large. Lower either. The check fires before anything
is allocated and before tracking starts, deliberately: raised per window instead,
it would return a complete result with every window zero and nothing saying why.

**Need full patent Figure 0.5 output** — add `--geocode FILE.csv` to the
`--patent-exact` command. Without it, the run still applies the interpreted
Eq. 21–24 transform but does not write block 11.

**A uniformly zero tomogram** — no window survived the coherence mask. Re-run
the motion stage with `--coherence 0` to see what is actually there.

**Focusing is slow** — cost scales with grid cells times pulses. Halving `--size`
quarters the work. Use `--subap paper`, which focuses once and filters, rather
than `--subap pulse`, which focuses every frame separately.
