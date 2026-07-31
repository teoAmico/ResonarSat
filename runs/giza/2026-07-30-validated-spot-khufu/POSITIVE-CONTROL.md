# Positive control, and the bisection that followed

Synthetic, so this is not a run directory in the sense `runs/README.md` means --
it is kept beside the run whose central claim it withdrew.

## Why

`RUN.md` claimed its Khufu null was attributable to the target rather than the
configuration. Nothing had ever put that configuration in front of a target
known to vibrate.

## The fixture

`sim_cphd` collect, 20 s dwell, 8000 pulses, one vibrating point at the scene
centre. Target 0.5 Hz, 4 mm vertical, giving 4.19 px peak-to-peak of azimuth
displacement on a 0.4 m grid -- 69 times the 0.061 px quantisation floor at
`--upsample 40`. Observation ratio 0.50, inside the 0.39-0.69 range the published
validation operates at. A cluttered variant adds 400 Rayleigh scatterers.

Grid 320 x 320 at 0.4 m, so the target falls in windows 180, 181, 199 and 200.
**Those four windows are the measurement.** Anything reported elsewhere in the
scene is not a recovery of the target, for reasons the next section makes
quantitative.

## The bisection

Seven configurations against the same target, one thing changed at a time.
`t_sap` is held at 1.0 s wherever the comparison requires it -- the two routes
share `denom = N - (N-1)*Omega`, so equal `(N, Omega)` gives equal `t_sap` and
`dt` and the route can be varied alone.

| configuration | quality at target | p-p excursion | frequency at the four target windows |
|---|---|---|---|
| C1 uniform N=159 ov=.88 up=40 clean | 0.08-0.12 | 10-32 px | 1.20 1.57 1.05 0.79 — miss |
| C4 uniform N=159 ov=.88 up=**10** clean | 0.08-0.12 | 10-32 px | identical to C1 — miss |
| C3 uniform N=**33** ov=**.40** up=40 clean | 0.07-0.11 | 9-32 px | 0.46 0.10 0.77 0.77 — marginal |
| C2 **pulse** N=159 ov=.88 up=40 clean | **0.24-0.28** | 22-32 px | 0.31 0.31 **0.52** 0.31 — marginal |
| C0 **pulse** N=128 ov=**.00** up=40 clean | 0.07-0.09 | 28-32 px | 0.30 0.30 0.20 2.97 — miss |
| K1 uniform N=159 ov=.88 up=40 **clutter** | **0.28-0.32** | 26-32 px | 0.11 0.11 0.05 0.05 — miss |
| K2 pulse N=159 ov=.88 up=40 **clutter** | 0.15-0.16 | 29-32 px | 0.31 0.31 0.79 0.31 — miss |

### What is ruled out

- **Sub-pixel upsampling.** C4 against C1 is not merely similar, it is
  numerically identical to three decimals. 1/10 px and 1/40 px produce the same
  answer, so the refinement is not what fails.
- **Overlap.** C3 at 0.40 with `t_sap` held behaves like C1 at 0.88.
- **The route, as a cause.** The pulse route triples tracking quality on the
  clean scene (0.24-0.28 against 0.08-0.12) and the uniform route triples it on
  the cluttered one (0.28-0.32 against 0.15-0.16), so the routes differ -- but
  neither recovers the target in either scene.
- **Clutter as the missing ingredient.** `microm.c` warns that isolated point
  targets on an empty scene score low "even when tracking perfectly", so a bare
  point was suspected of being the wrong fixture. Adding 400 scatterers raises
  quality by a factor of three and moves the reported frequency to the LOWEST
  bins, 0.05-0.11 Hz. Better coherence, worse answer.

### What is left

**The correlation tracking stage itself.** In all seven configurations the
peak-to-peak excursion at the target sits between 22 and 32 pixels in a
32-pixel window -- the correlation peak wandering the full patch -- against an
injected 4.19 px. That signature is invariant to route, overlap, upsampling and
scene content, and it is the same signature `RUN.md` recorded over Khufu
(median 31.98 px in a 32 px window).

`tests/test_coreg.c` passes, so the correlator recovers known sub-pixel offsets
between ideal patches. The failure is therefore in what reaches it -- patch
extraction, the reference look, or the sub-look images themselves -- and not in
the primitive. That is the next thing to bisect and it is not done here.

## The finding that matters most: the existing positive control is vacuous

`tests/test_nullmotion.c` scores a detection by searching **every** window for
one whose dominant frequency lands within two bins of the injected value, and
comparing its prominence against a static-scene floor.

Applied to this ladder at 361 windows, **all seven configurations "recover" the
target** -- including C1, which misses at every one of the four windows that
actually contain it. The recovering windows are 217, 107, 162 and 255, scattered
across the scene, none of them near the target.

The arithmetic is not subtle. A five-bin window out of eighty is 6.25% per
window, so among 361 windows:

```
expected chance matches      22.6
P(at least one match)        1.000000
```

A criterion satisfied with certainty by chance is not a criterion. The existing
test runs 9 windows rather than 361, where the same figure is P = 0.44 -- weak
rather than vacuous, and partly rescued by the comparison against a static
floor -- but the criterion does not become sound at any window count, it only
becomes less obviously unsound.

**Consequence, and the fix, which did not go the way this section predicted.**

Restricting the search to windows containing the target -- the fix proposed here
first -- turns out not to bite on this test's geometry: its grid is 64 cells and
its windows are 32 wide at stride 16, so all nine of them contain the centre
target. The criterion's weakness is not where the search looks but what it
searches FOR. It asks whether SOME window reports approximately the right
frequency, which is a question about a set of candidates.

The criterion that cannot be satisfied by a lucky window is what the tool would
PRINT: the dominant frequency of the most prominent window, one value per run,
selected by rs_spectrum_best_window(). That is now asserted alongside the old
measure, and **it passes**:

```
   injected  matched prom     floor  clears? TOOL reports
      0.3 Hz         24.4      16.4      yes    0.302 Hz ok
      0.5 Hz         32.0      16.4      yes    0.504 Hz ok
      0.7 Hz         30.6      16.4      yes    0.706 Hz ok
      0.9 Hz         17.9      16.4      yes    0.907 Hz ok
```

So the sentence this section originally carried -- that nothing in the
repository demonstrates the chain recovers a correct frequency -- was wrong, and
is retracted. It does, at 128 looks, zero overlap, the pulse route, a 20 mm
target on a 64-cell grid. Four frequencies, four correct answers, each within
half a bin.

What remains true is that the demonstration does not extend to the operating
point this run used. The gap between the two is now the whole question, and the
most conspicuous difference is how much of the correlation window the target's
response fills: at 128 looks the sub-look resolution is 8.26 m against a 0.5 m
cell, so the target spans about sixteen cells of a thirty-two cell window; at
159 looks with 0.88 overlap it is 1.27 m against 0.4 m, about three cells of
thirty-two. That is the next thing to vary, and it is not one of the four
candidates this bisection set out to test.


---

# The second bisection: what actually breaks it

The first bisection ruled out upsampling, overlap, route and clutter, and left
"the correlation tracking stage". Varying the remaining geometry isolates two
separate limits, and the second one has a closed form.

## Limit 1: amplitude

The working configuration (128 looks, zero overlap, pulse route) recovers a
0.5 Hz target at 20 mm and misses the identical target at 4 mm:

```
  C0, 20 mm:  win 199 -> 0.504 Hz, prominence 15.9   RECOVERED
  C0,  4 mm:  0.30 0.30 0.20 2.97 Hz                 missed
```

So there is a sensitivity threshold between 4 and 20 mm for that geometry --
between 4.19 and 21 px of azimuth displacement against a tracking noise that
saturates the search window. Ordinary, and worth knowing.

## Limit 2: the observation ratio, and it is not about signal level

The Khufu operating point misses at BOTH amplitudes:

```
  Khufu config, 20 mm:  1.569 1.569 1.569 0.785 Hz   MISSED
  Khufu config,  4 mm:  1.203 1.569 1.046 0.785 Hz   MISSED
```

Five times the signal changes nothing, so this is not a sensitivity limit. What
changes between the two configurations is what the target LOOKS like in a
sub-look.

A vibrating target does not image as a point. It images as a paired-echo train
with spacing `f * lambda*R/(2*v_p)`, and the sub-look resolves that train when
the spacing exceeds its own resolution `lambda*R/(2*v_p*t_sap)`. The ratio is:

```
   [ f * lambda*R/(2 v_p) ]  /  [ lambda*R/(2 v_p t_sap) ]  =  f * t_sap  =  ETA
```

**The ghost-resolution ratio is exactly the observation ratio.** Every case
measured falls in line:

| case | `t_sap` | f | eta | ghosts | result |
|---|---|---|---|---|---|
| passing control 0.3 Hz | 0.156 s | 0.3 | 0.047 | unresolved | ok |
| passing control 0.9 Hz | 0.156 s | 0.9 | 0.141 | unresolved | ok |
| C0, 20 mm | 0.156 s | 0.5 | 0.078 | unresolved | **recovered** |
| C0, 4 mm | 0.156 s | 0.5 | 0.078 | unresolved | missed — amplitude |
| Khufu config, 4 mm | 1.002 s | 0.5 | **0.501** | **resolved** | missed |
| Khufu config, 20 mm | 1.002 s | 0.5 | **0.501** | **resolved** | missed |

Once `eta` approaches 0.5 the tracked feature fragments into a train of
comparably bright spots, the correlation surface acquires competing peaks of
similar height, and the argmax hops between them. The shift series shows exactly
that -- not wander but discrete jumps, returning repeatedly to the same handful
of offsets:

```
  window 180:  0.00  0.05  0.15  13.12  15.25 -12.20 -14.38  16.52 -0.28  3.05  5.20 ...
```

and `quality`, which is the mean correlation PEAK HEIGHT (`microm.c:510`), sits
at 0.08-0.30 while the argmax spans the entire window. A modest peak that is not
the highest peak is precisely the signature of competing lobes.

This is not a defect in this implementation. It is the upper bound on `t_sap`
that Vattulainen et al. state qualitatively as their third competing effect --
"longer sub-aperture times can make the feature less distinct... since the
changing signal phase is integrated over a greater period" -- given a
quantitative form. Their validated runs sit at `eta` 0.39-0.69, right at the
boundary, and they report performance degrading across that range.

## What this means for the Khufu run

`t_sap` there was 1.643 s, so `eta` reaches 0.5 at **0.30 Hz**. Above that the
paired echoes resolve and correlation tracking degrades. Combined with the
displacement-averaging nulls at 0.609 Hz and above, the configuration's usable
band was roughly **0 to 0.3 Hz** -- not the 2.53 Hz its own output reported as
the observable band.

The run was misconfigured, in a way none of the four original candidates named
and which the printed band figure actively concealed.

## The most dangerous single result here

Narrowing the correlation window raises quality and prominence monotonically
while leaving the answer wrong:

| win | fill fraction | quality | p-p excursion | frequency | prominence |
|---|---|---|---|---|---|
| 32 | 0.10 | 0.08-0.12 | 10-32 px | 1.20 1.57 1.05 0.79 | 7-11 |
| 16 | 0.20 | 0.16-0.21 | 15.8-16.0 px | 0.42 0.26 0.42 0.79 | 9-26 |
| 8 | 0.40 | 0.30 | 7.7-7.9 px | **0.209 in all four** | **47-57** |

At `--win 8` all four windows containing the target agree on 0.209 Hz with
prominence near 50, against an injected 0.5 Hz. **Confident, spatially
coherent, and wrong** -- four adjacent windows concurring, high prominence, and
the excursion still saturating the window. Every heuristic this project uses to
decide a peak is real would pass it.

Note also that the excursion equals the window width at every size tried, which
is what a correlation argmax does when no peak dominates: the reported shift is
bounded only by the search extent.

---

# 2026-07-31: the observation-ratio threshold, attempted and withdrawn

The section above localised the failure to the observation ratio `eta = f*t_sap`
and left the threshold as a bracket between two data points. This is the attempt
to measure it, and what it found instead.

## Design

Hold the modulation index `B` fixed at 1.648 -- 5 mm, the amplitude
`tests/test_pairedecho.c` verifies -- so the ghost train's *structure* is
constant, and vary only frequency. Then `eta = f*t_sap` is the sole variable.
N=159, overlap 0.88, `t_sap` 1.0020 s, `dt` 0.1202 s, `df` 0.0523 Hz, cell 0.4 m,
size 320, win 32, upsample 10. Analysis restricted to the four windows containing
the target. "RECOVERED" means the reported dominant frequency is within two bins
of the injected one.

## The ladder

| f (Hz) | eta | uniform reports | pulse reports |
|---|---|---|---|
| 0.2 | 0.200 | 1.569 | 0.314 |
| 0.3 | 0.301 | 1.569 | 0.314 *(within 2 bins)* |
| 0.4 | 0.401 | 1.569 | 0.314 *(within 2 bins)* |
| 0.5 | 0.501 | 1.046 | 0.314 |
| 0.7 | 0.701 | 1.569 | 0.314 |
| 1.0 | 1.002 | 1.569 | 0.367 |
| 1.4 | 1.403 | 0.105 | 0.314 |

**The reported frequency barely depends on the injected one.** The uniform route
answers 1.569 Hz at five of seven; the pulse route answers 0.314 Hz at six of
seven. The two apparent recoveries are that fixed value landing within tolerance
of the target -- 0.314 is within two bins of both 0.3 and 0.4 -- and are
coincidences, not detections.

## The static control

Run the identical chain on a target with no vibration at all:

| configuration | reports | prominence |
|---|---|---|
| pulse N=128 overlap 0 | 1.260 Hz | 10.3 |
| pulse N=159 overlap 0.88 | 0.157 Hz | 13.0 |
| **uniform N=159 overlap 0.88** | **1.569 Hz** | **27.9** |

The uniform route returns **the same 1.569 Hz, at the highest prominence in the
whole experiment, from a scene containing no motion.** It is a processing
artefact, and it is what the ladder above was reading at every frequency.

## Why the eta threshold is not measurable this way

The known-working shape -- pulse, N=128, zero overlap, `t_sap` 0.156 s -- was run
on the same collects and grid:

| f (Hz) | eta | p-p displacement | reports | |
|---|---|---|---|---|
| 0.2 | 0.031 | 2.1 px | 0.302 | miss |
| 0.3 | 0.047 | 3.1 px | 0.302 | coincidence |
| 0.5 | 0.078 | 5.2 px | 0.302 | miss |
| 0.7 | 0.109 | 7.3 px | **0.706** | recovered |
| 1.0 | 0.156 | 10.5 px | **1.008** | recovered |
| 1.4 | 0.219 | 14.7 px | **1.411** | recovered |

This shape does track the injected frequency, but only once the displacement
exceeds roughly 7 px peak-to-peak -- below that its own artefact at ~0.302 Hz
wins. So within this configuration the discriminator is **displacement**, not
`eta`.

And the two configurations do not overlap in `eta` at all: 0.031-0.219 against
0.200-1.403. They meet only at `eta ~ 0.2`, where 128/0 recovers (14.7 px) and
159/0.88 misses (2.1 px) -- a comparison confounded by a factor of seven in
displacement.

The confound is structural, not an oversight in the design. Displacement goes as
`f*A`, the modulation index `B` as `A`, and `eta` as `f*t_sap`. Fixing any two
determines the third. Matching `eta` and displacement across two values of
`t_sap` forces `B` to differ by a factor of six, which changes the ghost train
being tracked. **These three cannot be separated by varying the target.**

## What stands and what does not

Stands: the paired-echo train exists and sits where theory puts it
(`tests/test_pairedecho.c`); `eta = f*t_sap` is exactly the ratio of ghost
spacing to sub-look resolution; `eta` grows with dwell.

Does not: that `eta ~ 0.5` is where the tracker breaks, and that `eta` is what
defeated the Khufu run. Both were inferred from the confounded comparison above.
`RS_ETA_BAD` was demoted from FAIL to WARN in `src/core/validate.c` accordingly,
and the README's "nothing above about 0.12 Hz is reachable" was withdrawn.

## The finding that displaced it

A fixed spurious frequency, configuration-dependent and present with no motion in
the scene, at prominence exceeding anything the real signal produces. Separating
`eta` from displacement needs a target-independent handle -- varying `t_sap` at
fixed `f` and `A`, which moves `eta` alone -- but that is not worth running until
the artefact is understood, because it is what both routes are reporting.

---

# 2026-07-31, continued: what the controlling variable actually is

The withdrawal above left the artefact as the thing to chase. Chasing it found
the constraint that governs the whole chain, and it is not `eta`.

## The artefact is a saturating argmax

The static-scene shift series is a sawtooth pinned to the edges of the search
extent -- `+-15.6 px` against a `--win 32` half-width of 16 -- clean for the
first ten looks and wrapping thereafter. The reported "frequency" is the wrap
rate. Nothing about the scene enters it, which is why it is fixed.

## Two bounds, and they can cross

- **Ceiling.** `rs_microm_recommend_looks()` already requires the peak shift to
  stay inside three quarters of a **sub-look** resolution cell. Driving a target
  up through it at 0.5 Hz: 1x and 2x recover, 4x reports **1.512 Hz -- exactly
  three times the true 0.504** -- and 16x collapses to the lowest bin. A wrapping
  sawtooth generating odd harmonics is the expected signature. Not monotone, so
  the ceiling is a boundary of reliability rather than a cliff.
- **Floor.** Below roughly 7 px peak-to-peak the artefact wins.

Sub-look resolution is `lambda*R/(2*v_p*t_sap)`, so a **short** sub-aperture
raises the ceiling. Overlap shortens nothing -- it widens each look's band and
therefore *lowers* it. For the three configurations run:

| configuration | sub-look res | ceiling | floor | window |
|---|---|---|---|---|
| pulse N=128 ov=0 | 8.13 m | 30.5 px | 7 px | 7-30.5 px |
| pulse N=159 ov=0.88 | 1.27 m | 4.8 px | 7 px | **none** |
| uniform N=159 ov=0.88 | 1.27 m | 4.8 px | 7 px | **none** |

This predicts **every row of both ladders above**: the 128/0 recoveries sat at
7.3, 10.5 and 14.7 px and its misses at 2.1, 3.1 and 5.2 px.

## Verified in both directions

**Fixed displacement, sweeping frequency.** Scaling amplitude as `1/f` holds the
excursion constant while `eta` varies 14-fold:

| f (Hz) | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 | 1.4 |
|---|---|---|---|---|---|---|---|
| eta | 0.016 | 0.031 | 0.047 | 0.078 | 0.109 | 0.156 | 0.219 |
| reported | 0.101 | 0.202 | 0.302 | 0.504 | 0.706 | 1.008 | 1.411 |

**7 of 7**, including four frequencies that missed at fixed amplitude. Across
this range the controlling variable is displacement, not `eta`.

**The crossed window really is empty.** Four targets placed *under* the
0.88 configuration's 4.8 px ceiling -- 1, 2, 3 and 4 px p2p -- all miss and all
report 1.569 Hz. Prominence falls monotonically (26.2, 22.1, 18.1, 12.6) as the
real signal grows without ever displacing the artefact. Its floor is above its
ceiling, so nothing can be measured with it at any amplitude or frequency.

## Two consequences worth stating separately

**The aperture fraction does not transfer.** Giza at alpha 5% over a 32.9 s
dwell gives `t_sap` 1.64 s, a 1.03 m sub-look and a 1.5 px ceiling -- unusable
at every frequency. Opening a window needs alpha near 0.8%, far below the
4.5-7.6% the published campaigns validate. Those campaigns collect much shorter
dwells, so **`t_sap` in seconds is the transferable quantity, not alpha.** The
validator now reports both and lets them disagree.

**The sensitivity figure was optimistic by a factor of 57.** The old check
reported `2.449/upsample` px -- the sub-pixel *interpolation* limit -- as the
detection floor. The measured floor is 7 px p2p. At 0.1 Hz on Giza that moves
the smallest visible vertical amplitude from 1.2 mm to **68.7 mm**, and even a
short sub-aperture only gives an admissible band of 69-95 mm. Seven centimetres.
Both numbers are now printed, with the interpolation limit labelled as such.

---

# 2026-07-31, third pass: the floor depends on the scene, and prominence does not discriminate

The 7 px floor above -- which sets the whole sensitivity figure -- was measured
against an isolated point on an empty background. `sim_cphd --help` says not to
do that: correlation tracking "is biased by its own point response when the
window is otherwise empty," and `--clutter` is "whenever the question is about
the tracker rather than about focusing." The question here was about the
tracker. Repeating it properly.

## The floor, both ways

Pulse route, 128 looks, zero overlap, 0.5 Hz, sweeping amplitude at 0.4 m cells.

| px p2p | isolated point | vibrating clutter patch |
|---|---|---|
| 1 | 1.260 | 2.671 |
| 2 | 1.260 | 2.369 |
| 3 | 0.302 | 0.050 |
| 4 | 0.302 | **0.504** |
| 5 | 0.302 | **0.504** |
| 6 | 0.302 | **0.504** |
| 7 | **0.504** | **0.504** |
| 10 | **0.504** | **0.504** |

Injected 0.500 Hz. **Isolated 7 px, textured 4 px** -- a factor of 1.75. The
sensitivity figure is scene-dependent by that factor: 68.7 mm at 0.1 Hz on Giza
becomes about 39 mm over textured ground. Still centimetres.

## A marginally open window is not an open window

Lowering the constant to 4 px would give the uniform 159/0.88 configuration --
ceiling 4.8 px -- a nominal 4.0-4.8 px window. Targets placed at 4.0, 4.4 and
4.8 px inside it, with clutter, **all miss**, reporting the lowest one or two
spectral bins. A ceiling has to clear the floor by a real margin; a ratio of 1.2
does not, and 4.4 (the working configuration's) does. Where between those the
boundary lies is not measured, and is not being guessed at. The validator keeps
the conservative 7 px, which refuses this configuration for the right reason.

## Prominence does not separate right answers from wrong ones

Worth stating on its own, because it is a heuristic this project leans on:

| case | reported | injected | prominence |
|---|---|---|---|
| isolated, 1 px | 1.260 | 0.500 | 8.8 |
| isolated, 2 px | 1.260 | 0.500 | 9.6 |
| textured, 4 px | **0.504** | 0.500 | **4.6** |
| textured, 5 px | **0.504** | 0.500 | 6.3 |
| uniform 0.88, 4 px | 0.105 | 0.500 | **27.7** |

The correct answers carry prominence 4.6 and 6.3. The wrong ones carry 8.8, 9.6
and 27.7. **Prominence is anti-correlated with correctness across this set.** It
measures how cleanly a peak stands out, and a saturating argmax produces a very
clean peak at the wrong frequency -- cleaner than a real signal near the floor.
It cannot be used to decide whether a detection is real.

---

# 2026-07-31, fourth pass: what the artefact actually is

## A guard that was built and did not work

The sawtooth pinned to `+-win_az/2` looked like the mechanism, so a saturation
metric was added to `rs_microm_t`: per window, the fraction of looks whose shift
sits within a pixel of the search extent. **It does not discriminate and was
removed.** Over the static uniform run, the 56 windows reporting the 1.569 Hz
artefact have a *lower* mean saturation (0.023) than the 305 that do not
(0.052), and both have median zero. The pinned sawtooth is real in the windows
where it appears, but it is not what produces the artefact across the grid.
Recorded because the metric is an obvious thing to reach for twice.

## What it is: the band was computed from the wrong quantity

The step between sub-apertures sets how finely the series is sampled. It does
not set what the series can carry. Each sub-aperture **averages** the motion
over its own duration, a lowpass whose first null is at `1/t_sap` -- the same
null `RS_VALIDATE_AVERAGING_NULL` already reported without the connection being
drawn. So the band is `1/(2*t_sap)`, and overlap buys resolution in time rather
than bandwidth. Quoting `1/(2*dt)` overstates it by `1/(1-overlap)`.

For the configuration that produces the artefact: `t_sap` 1.002 s, first
averaging null 0.998 Hz, band 0.499 Hz. **The 1.569 Hz it reports on a
motionless scene is past the first null**, where the sub-aperture has
essentially no response. Nothing measured can live there. Its dt-based `f_max`
of 4.158 Hz called that comfortably in band, which is why it went unnoticed
through an entire ladder.

The working configuration is untouched: `t_sap` 0.156 s gives a 3.20 Hz band,
and all seven frequencies it recovered fall inside it. The check is not
rejecting everything.

## And it puts eta = 0.5 back, derived

`eta = f*t_sap`, so `f < 1/(2*t_sap)` **is** `eta < 0.5`. The bracket withdrawn
earlier was the right number reached by the wrong argument: it is not about
resolving paired echoes, it is the sub-aperture's own averaging response, and it
follows from arithmetic rather than from the confounded comparison. The
paired-echo mechanism is still real and still confirmed by
`tests/test_pairedecho.c`; it simply is not what sets this bound.

## Retroactive

`runs/giza/2026-07-30-uniform-phase-khufu` configuration A ran at `--overlap
0.99`, where `denom = N - (N-1)*overlap` collapses to 2.27 and each sub-aperture
spans 14.5 s of a 32.9 s dwell. Its band is 0.0345 Hz. **Every frequency it
reported is above it** -- 0.054 by 1.6x, 0.324 by 9.4x, 1.026 by 30x -- and that
is the mechanism behind its own P5 observation that the peak follows `--fmin`.
Amended in that run's RUN.md. Configuration B is inside its band and unaffected.

---

# 2026-07-31, fifth pass: how the floor scales is not known

`resonarsat validate` was run on the real collect for the first time. It works
end to end -- 13 checks, `VERDICT: FAIL` for a 5 mm target at 0.1 Hz, because
that is 0.1x the floor it computes. Everything below came from trying to answer
the question that verdict raises: **is there any operating point on this collect
that passes?**

## The optimisation, and why it was wrong

Both amplitude bounds move with configuration:

```
floor   = FLOOR_PX * cell / k(f)          k(f) = R*2*pi*f*cos(theta)/v
ceiling = 0.75 * res_sap / k(f)
```

Read that way the floor scales with the CELL and the ceiling does not, so a
finer grid buys sensitivity for free. Sweeping alpha and cell over the Giza
geometry with `f` at the eta = 0.20 limit, the best point came out at **0.3 mm
vertical at 3.04 Hz** (alpha 0.002, cell 0.125 m) with a 44x window -- which
would have made the collect capable of the published claim after all.

It rests entirely on the floor being a fixed pixel count. That was measured at
one cell size and never tested.

## Tested, and refuted

Holding the physical scene at 128 m and the physical window at 12.8 m, so only
the cell changes:

| target p2p | cell 0.2 m (win 64) | cell 0.4 m (win 32) |
|---|---|---|
| 1.2 m | 0.958 miss | 0.756 miss |

At the finer cell 1.2 m is **6 px**, comfortably above the 4 px textured floor,
and it misses. Two clean misses where the pixel reading predicts recovery.
**A finer cell does not buy sensitivity, and the 0.3 mm figure is withdrawn.**

An earlier sweep in the same design confounded cell with window size in pixels
-- halving the cell at fixed physical window drops `--win` from 32 to 16, and
narrow pixel windows were already recorded above as producing confident wrong
answers. The 0.8 m column of that sweep is uninterpretable for that reason.

## The obvious replacement is also refuted

If the floor were a fixed fraction of the sub-look resolution cell it would
transfer across collects, which is what a validator wants. Varying `n_looks` at
fixed cell is suggestive -- 1.2 m p2p recovers at 64 looks (res 4.06 m) and
misses at 128 and 256 (8.13 m, 16.26 m).

But it fails against data already in hand. Both bounds would be proportional to
`res_sap`, making the window always open at a fixed 7.5x ratio, while the
uniform route at 159 looks and 0.88 overlap was measured **closed** -- targets
at 4.0, 4.4 and 4.8 px all missing inside the window it predicts. And the
64-look recovery is weak: `df` is 0.050 Hz there and it reported 0.450 against
0.500, a full bin off, where every clean recovery in this file reads 0.504.

## Where it stands

- pixel floor: **refuted** (two clean misses at the finer cell)
- fraction-of-`res_sap` floor: **refuted** (predicts an open window measured closed)
- fixed in metres: survives, but only by discounting the marginal 64-look result

`src/core/validate.c` implements the pixel reading. Its comment now says the
scaling is refuted and that figures derived at other cells are not established;
the constant itself is unchanged, because it is still what was measured at the
cell it was measured at.

**Cell size, sub-look resolution and window size in pixels are confounded in
every run made so far**, including all of today's. Separating them needs a
designed experiment that varies each independently -- which means accepting
physical windows smaller than a sub-look resolution cell in some arms, and
correlation windows of 128 px or more in others. Until that is done, the
question this pass opened -- whether any operating point on this collect
passes -- has no answer.

## Addendum: the 0.88 evidence, redone in band

The runs above that refuted the fraction model were all injected at 0.500 Hz
against a configuration whose band edge is 0.499 Hz -- eta 0.501, exactly at the
limit. That was a fair objection to them, so they were redone at 0.200 Hz,
eta 0.20, well inside the band.

The standard offset-tracking CRB (Bamler & Eineder 2005) is in units of the
resolution cell and falls as `1/sqrt(N)` with `N` the independent samples in the
window. That predicts the 0.88 configuration should be **41x better** than the
working one -- `res/sqrt(N)` of 0.126 m against 5.164 m -- because its finer
sub-look packs 102 independent samples into the same 12.8 m window against 2.5.

Excursion fixed at 2 px, frequency swept, all in band:

| injected | true bin | reported | prominence |
|---|---|---|---|
| 0.15 Hz | 2.87 | 0.157 | 116.3 |
| 0.20 Hz | 3.82 | 0.105 | 145.9 |
| 0.30 Hz | 5.74 | 0.209 | 66.1 |
| 0.40 Hz | 7.65 | 0.157 | 27.1 |

Reported values sit in bins 2-4 while the truth runs 3 to 8. **It does not
track.** The apparent recoveries at 1, 2 and 3 px in the amplitude sweep were a
low-bin artefact landing near a low injected frequency, at prominences up to
330 -- the anti-correlation between prominence and correctness again.

So the conclusion survives the objection, and the CRB reading is eliminated
cleanly: it predicts recovery at 1 px and the configuration tracks nothing at
any excursion tried.

**This also moves overlap into the frame.** The failing and working
configurations share cell, window and pixel window size; they differ in
sub-look resolution, independent-sample count, route and overlap. Theory says
the failing one should be better on the first two. Both routes fail at 0.88
overlap -- uniform reports a fixed 1.569 Hz, pulse a fixed 0.314 Hz -- so the
route is not the discriminator either. High overlap is implicated on its own,
and by a mechanism none of the three floor models describes.

The designed experiment therefore needs four factors, not three: cell,
sub-look resolution, pixel window size, and overlap.

---

# The criterion that would have caught all of it

Every withdrawal in this file came from the same mistake: judging recovery by
`|reported - injected| < 2 bins` at ONE frequency. A chain emitting a fixed
spurious frequency passes that wherever the fixed value happens to fall near
the injection, which at these bin spacings is a wide range. The four cases:

| configuration | reports | scored "recovered" at |
|---|---|---|
| uniform 159 / 0.88 | 1.569 Hz for 0.2-1.4 Hz injected | wherever tolerance reached |
| pulse 159 / 0.88 | 0.314 Hz for six of seven | 0.3 and 0.4 Hz |
| pulse 128 / 0 | 0.302 Hz below the floor | 0.3 Hz |
| uniform 159 / 0.88 in band | bins 2-4 for bins 3-8 | 0.15 Hz |

None is a detection. All four passed the per-point test.

`rs_track_fit()` in `tests/rs_test.h` replaces it. Sweep the injected frequency,
fit the reported against it, and check two numbers:

- **slope** -- 1 for a chain that follows the target, 0 for a fixed artefact.
  No single point can distinguish these, which is the whole problem.
- **rms** of `reported - injected`, bounded at half a bin. Four times tighter
  than the per-point tolerance, because a chain that genuinely tracks has no
  reason to be looser.

Both are asserted in `tests/test_nullmotion.c`, alongside a **negative control**
that feeds the criterion a constant 1.569 Hz and requires it to fail -- so the
criterion cannot quietly stop discriminating.

At the working operating point: slope 1.008, rms 0.0052 Hz against a 0.0252 Hz
bound, so it passes with five times the required margin. The fixed 1.569 Hz
series gives slope 0.000 and rms 0.9945 Hz.

**Every result in this file that was later withdrawn fails the new criterion,
and the one operating point that survives passes it comfortably.** Future
sweeps should report slope and rms rather than a count of per-point matches.

---

# Overlap isolated, and a fifth withdrawal

Overlap was the one factor worth testing before the full four-factor sweep: the
failing 0.88 configuration beats the working one on sub-look resolution and on
independent-sample count, so theory favours it, and it fails anyway. If overlap
alone accounts for that, most of the sweep is unnecessary.

It isolates cleanly. Holding `denom = n - (n-1)*overlap` fixed pins `t_sap`, and
with it the sub-look resolution, the observable band and the ambiguity ceiling.
Only the number and spacing of looks then changes:

| overlap | n_looks | denom | t_sap | dt |
|---|---|---|---|---|
| 0.00 | 128 | 128 | 0.1562 s | 0.1562 s |
| 0.50 | 255 | 128 | 0.1562 s | 0.0781 s |
| 0.75 | 509 | 128 | 0.1562 s | 0.0391 s |

Excursion held at 10 px p2p, frequency swept 0.3 to 1.1 Hz, scored with
`rs_track_fit()`:

| overlap | slope | rms (Hz) | verdict |
|---|---|---|---|
| 0.00 | -0.050 | 0.4696 | does not track (0.050 reported at 1.1 Hz) |
| 0.50 | **+1.004** | **0.0030** | tracks, 5 of 5 |
| 0.75 | +0.955 | 0.0359 | does not track (one bin low at 0.5 and 0.9 Hz) |

Read alone this says 0.5 overlap is markedly better than zero -- which would
contradict `rs_microm_reference_t`'s stated reasoning and the project's default.

## It does not survive replication

The zero-overlap verdict rests entirely on one anomaly at 1.1 Hz. Repeating that
point over three clutter realisations, with 0.5 overlap on the identical scenes:

| seed | overlap 0.00 | overlap 0.50 |
|---|---|---|
| 7 | 0.050 miss | 1.104 ok |
| 11 | 1.100 ok | 1.104 ok |
| 23 | 1.150 ok | **6.275 miss** |

Both fail on some realisations. Zero overlap's failure is seed-specific, and
0.5 overlap -- which scored slope 1.004 and rms 0.0030 on seed 7 -- reports
6.275 Hz against an injected 1.1 Hz on seed 23. **The overlap finding is
withdrawn; nothing here establishes an ordering between them.**

## What that means for the criterion

`rs_track_fit()` was added earlier today precisely to stop single-point matches
producing false conclusions, and it did: every earlier withdrawal fails it. But
it sweeps FREQUENCY and not SPECKLE, and a five-point sweep on one realisation
still passed a configuration that fails on another.

The criterion's documentation now says so. A passing slope and rms establish
that a chain tracked *on that scene*; a claim about a *configuration* needs the
sweep repeated over independent `--seed` realisations and the verdicts pooled.

The four-factor experiment therefore has a fifth dimension, and it is not
optional: **realisation**. Any arm evaluated on one seed can say nothing.
