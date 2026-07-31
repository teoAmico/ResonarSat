# Run: 2026-07-30 giza / uniform-phase-khufu

**Question this run is meant to answer:** Does the independently validated chain
-- uniform filter bank, heavy overlap, pixel-phase estimator -- recover a
vibration at Khufu where the patent chain measured a complete null?

- git commit: `8decc82`
- started:    2026-07-30T14:32:38Z
- host:       Darwin arm64 (8 cores, 25.8 GB RAM), OpenMP enabled

## Why this run exists

The two runs before this one (`2026-07-29-patent-exact-true-khufu`,
`2026-07-30-no-optimize-khufu`) both measured a complete null at the Great
Pyramid, and both used the patent's own chain: `--subap paper --reference pair
--estimator correlation`. That null was *predicted analytically* and has nothing
to do with the scene. In paper mode `N_D * dt = t_sap` identically, so every
resolvable bin sits at an integer observation ratio, where a
displacement-averaging observable has no response.

So those runs bound the patent's method, not this software's sensitivity, and
they say nothing about whether micro-motion is measurable at Giza by the route
that has published accelerometer validation. That route differs in all three
places that matter:

| | patent chain (previous runs) | this run |
|---|---|---|
| sub-apertures | `paper` — hold out B_DL, sweep the remainder | `uniform` — plain filter bank with overlap |
| reference | `pair` — master vs its own swept slave | not applicable |
| observable | correlation peak position (displacement-averaging) | **pixel phase** (no averaging) |

`RS_MICROM_EST_PHASE` is the observable of Clemente et al. (EuRAD 2025), the only
one in `microm.h` with published accelerometer validation on this class of free
X-band data, and the only one that is not blind at integer observation ratios --
that work recovers 36 Hz at an observation ratio of exactly 18, where the
averaging model predicts precisely zero.

## Collect

```
CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd  ->  data/giza.cphd
  Capella open data, CC BY 4.0
  39,180,270,944 bytes (verified)
  dwell 32.869 s   335,141 pulses x 25,073 range bins   PRF 10,196.35 Hz
  carrier 9.3000 GHz   wavelength 0.0322 m
```

Same collect, same grid as both previous runs: `--at 29.979175,31.134186` ->
`--offset -152,-552`, Khufu. Not a hand-typed offset.

## The band this chain can reach, and why there are two configurations

For a uniform bank of `N` looks at fractional overlap `F` over a dwell `T`, the
band layout in `rs_subaperture_split()` gives `dt = T(1-F)/(N(1-F)+F)`, so

```
f_max = 1/(2 dt) = N/(2T) + F/(2T(1-F))
```

The first term needs look count; the second needs overlap pushed very close to 1,
which buys sampling rate only by making every look nearly the whole aperture.

**`--subap uniform` is capped by the grid, not by the collect.** It splits a
focused image, and `rs_subaperture_split()` requires `n_az >= 2*n_looks`
(`src/core/subaperture.c:192`). At `--size 256` that is `--n 128` maximum, which
puts `f_max` at 3.45 Hz even at `--overlap 0.99`. Every validated result this
estimator comes from is one to two orders of magnitude above that: 36 Hz from a
ship's engine, 87 Hz from an idling van.

Reaching those needs `N ~ 2 T f_max`, thousands of looks, and therefore a grid
thousands of cells wide in azimuth -- tens of millions of cells to backproject.
The `pulse` route is the same `RS_SUBAP_UNIFORM` mode built from pulse windows
instead of Doppler band-passes, carries no `n_az` constraint, and costs
`n_cells * n_pulse/(1-F)` rather than a grid blow-up. So:

- **A** is the configuration as asked for, literally: `--subap uniform`.
- **B** is the same mode at the look count the validated band needs, via the
  route that can afford it.

B also happens to be the better-conditioned image: its sub-looks are ~0.16 s
long, so their azimuth resolution is around 10 m and a 2.0 m cell oversamples
them. A's sub-looks are 14.5 s long and resolve 0.115 m, which a 2.0 m cell
aliases -- the same warning both previous runs carried, kept identical here for
comparability.

## Commands

Paths shortened to the run directory as `$D`; every option is otherwise verbatim.

```sh
D=runs/giza/2026-07-30-uniform-phase-khufu

# A -- uniform filter bank, as asked. f_max 3.45 Hz.
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 256 --cell 2.0 \
    --subap uniform --overlap 0.99 --estimator phase \
    --n 128 --win 32 --coherence 0.4 \
    --shifts $D/A_uniform_n128_phase.csv --out $D/A_uniform_n128_phase

# A-null -- identical, sub-look time order destroyed.
#   (A + --shuffle-looks 20260730, outputs renamed Anull_*)

# A2 -- A with a 32-shuffle null floor instead of one shuffle.
#   (A + --null-trials 32, no --shifts/--out)

# The fmin ladder -- A2 with the search's lower limit moved up.
#   (A + --fmin 0.3 --null-trials 32), then --fmin 1.0, then --fmin 2.0

# B -- uniform mode at the validated band. f_max 31.47 Hz.
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 256 --cell 2.0 \
    --subap pulse --overlap 0.9 --estimator phase \
    --n 2048 --win 32 --coherence 0.4 --null-trials 16 \
    --shifts $D/B_pulse_n2048_phase.csv --out $D/B_pulse_n2048_phase
```

`--null-trials` reuses the stack already built, so its trials cost no refocus --
which is why the shuffle floors here are 32 and 16 samples rather than one.

Two analysis scripts sit in this directory because RUN.md quotes numbers the CLI
does not print, and a figure whose derivation is not reproducible is the thing
`runs/README.md` exists to prevent:

```sh
python3 shift_stats.py A_uniform_n128_phase.csv   # fold-ceiling and unwrap check
python3 freq_map.py    A_uniform_n128_phase.csv 0.3   # per-window frequency map

gunzip -k B_pulse_n2048_phase.csv.gz              # B's series, 225 x 2048 looks
python3 shift_stats.py B_pulse_n2048_phase.csv
```

B's shift series is stored gzipped: 225 windows x 2048 looks is 29.7 MB of text,
five times every other run in `runs/` combined, and it compresses to a third of
that. A's two are small enough to leave alone.

`freq_map.py` reproduces `rs_spectrum_compute_opts()`'s path -- linear detrend,
periodogram of LOS displacement, dominant bin above a floor -- so the fmin ladder
can be run offline from a `--shifts` CSV without refocusing.

**It is not a reimplementation of the CLI's window selection, and does not agree
with it on which window wins.** At `--fmin 0.3` the CLI names window 217 and the
script names window 43, both at 0.324 Hz and prominence ~55. The frequencies and
their spatial pattern are what it is used for here and those match; the winner
differs because `rs_spectrum_best_window()` applies the coherence gate and its
own prominence definition, and the script applies neither. Read it for the map,
not for the ranking.

## THE PREDICTION, recorded before the result

**This chain will report a frequency, and the frequency will not survive its null
test.** Those are two separate claims and the second is the one that matters.

It will report *something* where the patent chain refused, for a reason that is
structural rather than encouraging: the phase estimator has no correlation
surface, so `quant_px` is zero and the quantisation floor that produced
`RS_ERR_RANGE` twice does not apply to it at all. Only the coherence gate runs.
A chain that cannot decline is not more sensitive than one that can.

The failure mode expected is the one `microm.h` documents against this exact
configuration. Sub-look coherence is very nearly the fraction of pulses two looks
share, it levels off near 0.9 rather than approaching 1, and the unwrap
accumulates per-step noise as a random walk of `sigma*sqrt(N)`. Over a 33 s dwell
that is tens of radians at every overlap tried, where a usable unwrap needs it
well under pi. Raising the overlap makes this *worse*, not better, and both
configurations here are heavily overlapped.

- **P1** Both A and B report a peak rather than `RS_ERR_RANGE`.
- **P2** The peak sits in the lowest few bins after linear detrending -- a random
  walk is red -- and its prominence is low, under about 3, so the run prints its
  own flat-spectrum warning.
- **P3** The tell fires: peak-to-peak line-of-sight velocity at or near
  `lambda/(2 dt)`, which is 111 mm/s for A and about 1.0 m/s for B. That is the
  ceiling the phase fold imposes; hitting it means the series wrapped rather than
  that the target moved.
- **P4** The shuffled null matches or beats the measurement. Destroying the time
  order does not lower a random walk's floor by much, because the floor is set by
  per-step noise, which shuffling preserves.

**What would falsify it:** a peak clear of the lowest bins, in a contiguous patch
of at least four windows (they overlap at the tracking stride, so a resolved
target must occupy a 2x2 block), with prominence well above the shuffled floor
*and* peak-to-peak velocity well below `lambda/(2 dt)`. All four conditions
together. Any one of them alone is producible by unwrap noise.

## A SECOND PREDICTION, recorded after A and A-null but before the fmin ladder

P4 was wrong. A reports prominence 60.3 and a 32-shuffle null floor of mean 10.8,
worst 13.9 -- the measurement beats every shuffle, empirical p = 0.0303. Taken at
face value that is the first result in this project to pass a shuffle null test.

It should not be taken at face value, and the reason is visible in where the peak
sits. A's peak is bin 1 of 65 (0.054 Hz = 1/18.5 s, one cycle over the observed
span). Re-running with `--fmin 0.3` moves it to 0.324 Hz -- which is bin 6, the
first bin above the cut -- at prominence 55.8 against its own null floor of 9.9,
the same 5.6x ratio. A peak that relocates to wherever the search is told to
start is not a mode; it is a red spectrum, and a shuffle cannot test one, because
shuffling destroys the very time-ordering that makes a series red. The measured
per-step phase noise says the same thing from the other side: median |step| is
0.052 rad in A and 1.878 rad in A-null, so the shuffle does not preserve the
noise it is supposed to hold constant -- it inflates it 36-fold.

- **P5** The peak follows the floor. At `--fmin 1.0` it lands at 1.026 Hz (bin
  19) and at `--fmin 2.0` at 2.052 Hz (bin 38) -- in each case the first bin
  above the cut -- with prominence again roughly 5-6x its own shuffled floor.

**What would falsify it:** the peak staying put at some frequency as the floor
moves past it, or the prominence-to-null ratio collapsing when the floor is
raised. Either would mean there is a real mode there that the lowest bins were
merely outshouting.

## Result

### A — the chain reports a frequency, and it is drift

```
sub-apertures: 128 looks, dt 0.1448 s
  observable band  f_max 3.45 Hz   AT sub-look resolution 0.12 m
tracked 225 windows (15 x 15); 189 pass the 0.40 coherence mask
spectrum taken of line-of-sight DISPLACEMENT
spectra: 65 bins, 0.0540 Hz resolution
strongest peak in window 30: 0.054 Hz, prominence 60.3, quality 0.807,
                             peak-to-peak velocity 111.4 mm/s
  170 of 225 windows were eligible for selection
  sub-aperture response 0.2585 (-11.8 dB) at an observation ratio of 0.78
```

**0.054 Hz is bin 1 of 65** -- one cycle over the 18.5 s the 128 look centres
actually span. And the window it came from is wrapped, not moving: its
peak-to-peak line-of-sight velocity is 111.4 mm/s against a fold ceiling
`lambda/(2 dt)` of 111.2 mm/s, which is 100.2% of it, and its largest phase step
between consecutive looks is 6.277 rad -- 2 pi to within 0.1%.

Across the whole grid the series splits cleanly in two:

| | value |
|---|---|
| median peak-to-peak LOS velocity | 1.1 mm/s (1% of the fold ceiling) |
| windows at or above the fold ceiling | **29 of 225** |
| windows whose largest phase step exceeds pi | **29 of 225** — the same ones |
| median largest phase step | 0.052 rad |
| median peak-to-peak LOS displacement | 4.35 mm |

So 196 windows hold a smooth, slowly varying phase, and 29 have wrapped outright.
The reported winner is one of the 29. The estimator has no way to tell the two
apart, and prominence actively prefers the wrapped ones.

### The shuffle null is passed, and does not mean what it looks like

```
NULL FLOOR from 32 shuffles of the sub-look time order:
  mean 10.8, sd 1.5, worst 13.9
  detection 60.3 is 5.60x the mean and 4.34x the worst null
  0 of 32 nulls reached it -- empirical p = 0.0303
```

This is the first measurement in this project to beat its shuffle floor. It is
still not a detection, for a reason the shuffled run's own shift series shows:
**shuffling does not hold the noise constant here.** Median largest phase step is
0.052 rad in A and 1.878 rad in A-null -- destroying the time order inflates the
per-step noise 36-fold and whitens a red series. A shuffle test asks "could this
peak arise if the frames were in any other order"; against a drift-dominated
phase series it instead compares a smooth series with a deliberately roughened
one, and the smooth one wins whatever it contains. The test is not valid at this
operating point.

### The peak has no home: it follows wherever the search is told to start

| `--fmin` | peak | window | prominence | null mean / worst | ratio to mean | p-p velocity |
|---|---|---|---|---|---|---|
| none | **0.054 Hz** (bin 1) | 30 | 60.3 | 10.8 / 13.9 | 5.60x | 111.4 mm/s (at the fold) |
| 0.3 Hz | **0.324 Hz** (bin 6, first allowed) | 217 | 55.8 | 9.9 / 12.8 | 5.61x | 1.2 mm/s |
| 1.0 Hz | **1.187 Hz** | 193 | 34.3 | 8.1 / 10.7 | 4.24x | 1.1 mm/s |
| 2.0 Hz | **3.075 Hz** | 88 | 14.1 | 6.3 / 8.6 | 2.25x | 3.0 mm/s |

Four different frequencies, four different windows, and every one of them clears
its own shuffled floor with the same empirical p = 0.0303. A mode does not move
when you change where the search begins.

The per-window map settles it. Recomputing the spectra from
`A_uniform_n128_phase.csv` above 0.3 Hz, **all 225 windows peak between 0.32 and
0.59 Hz**: 161 of them at exactly 0.324 Hz, the first allowed bin, and 214 of 225
within the first three. The entire grid -- pyramid, mastaba field and open desert
alike -- reports the lowest frequency it is permitted to report. That is a red
spectrum everywhere, which is what a slowly drifting phase gives, and it is not a
property of the ground.

The prominence-to-null ratio decaying monotonically as the floor rises (5.60,
5.61, 4.24, 2.25) says the same thing: the significance lives at DC and thins out
as the cut moves away from it.

### How the predictions did

| | predicted | measured |
|---|---|---|
| **P1** | reports a peak, not `RS_ERR_RANGE` | **held** — the quantisation floor does not apply to this estimator, as expected |
| **P2** | peak in the lowest bins **and** prominence under ~3 | **half wrong.** Bin 1, yes. Prominence 60.3, not under 3 |
| **P3** | peak-to-peak velocity at the fold ceiling | **held, exactly** — 111.4 vs 111.2 mm/s, largest step 2 pi to 0.1%. But it describes 29 of 225 windows, not all of them |
| **P4** | the shuffled null matches or beats the measurement | **wrong** — 0 of 32 shuffles reached it, p = 0.0303 |
| **P5** | at `--fmin` 1.0 and 2.0 the peak lands on the first allowed bin | **wrong in form, right in substance** — it lands at 1.187 and 3.075 Hz, not 1.026 and 2.052, but it still moves with the cut and never stays put |

**P2's failure is the useful one.** Prominence measures how concentrated a
spectrum is, and drift is maximally concentrated -- all of it in bin 1. So a high
prominence is not evidence against an artefact here; it is what the artefact
produces. The flat-spectrum warning at prominence < 3 cannot catch this, because
the failure mode is the opposite of flat.

**P5's falsification clause was mis-stated and should not be read as written.** It
named "the prominence-to-null ratio collapsing when the floor is raised" as
evidence of a real mode being outshouted. That is backwards: a ratio that decays
monotonically as the cut rises is exactly what a red spectrum gives, and the
measured 5.60 -> 5.61 -> 4.24 -> 2.25 confirms the drift reading rather than
refuting it. The clause as written would have been satisfied by the very thing it
was meant to rule out.

### What A establishes

The validated chain's front end runs on this collect and produces a well-formed,
fully populated measurement where the patent chain produced nothing at all: 189
of 225 windows coherent, a spectrum in every one, and a peak that survives 32
shuffles. **None of that is a detection.** The observable is a smooth phase drift
in 87% of windows and an outright 2 pi wrap in the other 13%, and the reported
frequency is set by the search's lower limit rather than by the ground.

Two things would have to change before this configuration could carry a
measurement: an unwrap that survives the dwell (`microm.h` predicts it does not,
and the 29 wrapped windows are that prediction landing), and a null test that
holds per-step noise fixed while destroying the time structure -- which the
shuffle demonstrably does not do.

### B — the validated band is reachable, and the phase there is pure wrap

```
sub-apertures: 2048 looks, dt 0.0159 s
  observable band  f_max 31.47 Hz   AT sub-look resolution 10.59 m
tracked 225 windows (15 x 15); 178 pass the 0.40 coherence mask
spectra: 1025 bins, 0.0307 Hz resolution
strongest peak in window 200: 0.830 Hz, prominence 89.8, quality 0.473,
                             peak-to-peak velocity 1981.7 mm/s
  193 of 225 windows were eligible for selection
  sub-aperture response 0.9713 (-0.3 dB) at an observation ratio of 0.13

NULL FLOOR from 16 shuffles of the sub-look time order:
  mean 20.5, sd 2.2, worst 24.7
  detection 89.8 is 4.37x the mean and 3.64x the worst null
  0 of 16 nulls reached it -- empirical p = 0.0588
```

The band arrives as designed: 31.47 Hz, which covers the 36 Hz and 87 Hz class of
validated result far better than A's 3.45 Hz, and the sub-aperture response is
-0.3 dB rather than -11.8 dB because a 0.16 s look averages away almost nothing.
The 10.59 m sub-look resolution is properly oversampled by the 2.0 m cell, so
none of A's aliasing applies.

**And the phase series is saturated wrap in every window that passes the mask.**

| | value | what it should be |
|---|---|---|
| median peak-to-peak LOS velocity | 2004.4 mm/s | below the 1013.3 mm/s fold ceiling |
| as a fraction of the fold ceiling | **1.978x** | under 1 |
| windows at or above the ceiling | **178 of 225** — every one that passes coherence | few |
| median peak-to-peak LOS displacement | **16.10 mm** | anything below lambda/2 |
| lambda/2, the full ambiguity interval | 16.1 mm | — |
| median largest phase step | **6.226 rad** | well under pi |
| windows with a step above pi | 178 of 225 | none |

The displacement is pinned at exactly lambda/2 in every window, and the velocity
at twice the fold ceiling. Those are not measurements that came out large; they
are the arithmetic of a phase uniformly distributed over [-pi, pi]. The per-step
noise itself exceeds pi, so there is no unwrap to attempt -- the series was
ambiguous from the first step, not merely after accumulating.

`microm.h` predicts this and B is where the prediction bites hardest. Fine time
sampling and per-look SNR are in direct conflict on one collect: reaching 31 Hz
means each look uses 1/206th of the aperture, and at 0.16 s the pixel phase on
this scene is noise. Median tracking quality is 0.454 against a 0.531 maximum --
nothing on this grid is a strong enough scatterer to hold phase across looks that
short.

An offline recomputation of the spectra from `B_pulse_n2048_phase.csv` -- the same
detrend-and-periodogram path the tool runs, so no refocus needed -- shows the same
floor-following as A, now across a 31 Hz band:

| search floor | where the top windows land | top prominence |
|---|---|---|
| none | 0.34 – 0.92 Hz | 78.5 |
| 2 Hz | 2.09 – 2.89 Hz | 57.6 |
| 10 Hz | 10.05 – 11.40 Hz | 12.9 |

Nothing sits still. In a 1025-bin spectrum reaching 31.47 Hz, the strongest peak
is always within a bin or two of wherever the search is told to begin, and its
prominence decays as that point moves up. Neighbouring windows agreeing is not
independent confirmation either: the tracking stride is half the window, so
adjacent windows share half their pixels by construction.

### What this run establishes

**The independently validated chain does not detect micro-motion at Giza on this
collect, in either configuration, and the two fail differently.**

- **A** (`--subap uniform`, 128 looks, 0.99 overlap) is phase-stable -- median
  step 0.052 rad -- but can only reach 3.45 Hz, and what it measures there is a
  drift whose reported frequency is set by the search's lower limit. 29 of 225
  windows wrap outright.
- **B** (2048 looks, 0.9 overlap) reaches 31.47 Hz, the band where the validated
  results live, and has no usable phase there at all: every tracked window is
  saturated at the ambiguity interval.

Between them they bracket the trade this collect imposes. There is no setting of
this chain on this data that has both the band and the phase stability a
micro-motion measurement needs.

**A caveat about the project's standing claim.** `README.md` says every
micro-motion measurement made with this implementation has failed its null test.
After this run that is no longer true as written: A clears 32 shuffles at
p = 0.0303 and B clears 16 at p = 0.0588. Neither is a detection, and the
evidence above is that the shuffle test is simply not valid against a
drift-dominated or wrap-saturated phase series -- it inflates the per-step noise
36-fold in A while claiming to hold everything but the ordering constant. The
honest statement is now narrower and worse for the method: *this chain can pass
the shuffle null test without measuring anything*, which makes the shuffle
insufficient as the project's credibility check for phase observables. The
`--null-static` test, which reproduces the processing's own artefacts rather than
just its noise, was not run here and is the one that should adjudicate.

### A THIRD PREDICTION, recorded before the static null

`--null-static 8` was then run on A: eight simulated motionless collects built on
this collect's own geometry, 400 scatterers each, put through the identical
chain. `rs_null_static()` is the test the code's own comment calls the one "an
overlapping decomposition cannot walk over the way it can walk over a shuffle",
and it is the right adjudicator here precisely because it reproduces the
processing's artefacts rather than only its noise.

- **P6** A motionless scene reaches prominence 60.3. `nge > 0`, and the run
  prints "A MOTIONLESS SCENE REACHED THIS MEASUREMENT". The drift is
  manufactured by the chain -- overlapping sub-apertures, unwrap, linear detrend
  -- and a world where nothing moves has all of those.

The way this could come out otherwise, stated in advance so it is not
rationalised afterwards: the simulated scene is 400 discrete scatterers, not
real clutter, so its per-look SNR may be far better than the pyramid's. Cleaner
phase means less drift, and a static floor well *below* 60.3 would be a statement
about the simulator's scene statistics, not a clean bill of health for the
measurement. If that happens the honest reading is that neither null available in
this tool is valid for this observable, not that the peak survived.

### Result of the static null — the two nulls disagree by 4.3x, and the honest one wins

```
strongest peak in window 30: 0.054 Hz, prominence 60.3

STATIC-SCENE NULL FLOOR from 8 simulated motionless collects
with this collect's own geometry, through the identical chain:
  mean 50.5, sd 4.6, worst 59.9
  detection 60.3 is 1.19x the mean and 1.01x the worst
  0 of 8 reached it -- empirical p = 0.1111
  No motionless realisation reached it.
```

Every trial, without exception:

```
  static trial 1/8: prominence 50.9 at 0.054 Hz     trial 5/8: 49.6 at 0.054 Hz
  static trial 2/8: prominence 52.8 at 0.054 Hz     trial 6/8: 51.1 at 0.054 Hz
  static trial 3/8: prominence 50.9 at 0.054 Hz     trial 7/8: 44.1 at 0.054 Hz
  static trial 4/8: prominence 59.9 at 0.054 Hz     trial 8/8: 44.7 at 0.054 Hz
```

**Eight simulated worlds in which nothing moves all produce a peak at 0.054 Hz --
the same bin as the measurement -- at 84% of its prominence on average and 99% of
it at worst.** The measurement stands 2.1 sd above the motionless mean, p = 0.11.

Put the two nulls side by side:

| null | floor (mean / worst) | measurement relative to it |
|---|---|---|
| 32 shuffles of the time order | 10.8 / 13.9 | **5.60x / 4.34x**, p = 0.03 |
| 8 simulated motionless collects | **50.5 / 59.9** | **1.19x / 1.01x**, p = 0.11 |

They disagree by a factor of 4.3, and `docs/USER-GUIDE.md` says which to believe:
the static test "is the one to trust when the two disagree, because it reproduces
the processing's own artefacts as well as its noise. Use it whenever frames
overlap heavily." At `--overlap 0.99` they overlap about as heavily as they can.

So the 0.054 Hz peak is manufactured by the chain. Not by the pyramid, not by the
desert, not by the collect -- by overlapping sub-apertures, phase unwrapping and
a linear detrend applied to a series 18.5 s long, all of which a motionless world
has too. This is the same conclusion the fmin ladder and the all-225-window
frequency map reached, now from a third direction and with a simulated ground
truth behind it.

**P6 was wrong in letter and right in substance.** It predicted `nge > 0` and the
"A MOTIONLESS SCENE REACHED THIS MEASUREMENT" warning. No trial reached 60.3, so
the tool printed "No motionless realisation reached it" instead -- which, read
alone, sounds like the measurement survived. It did not: the margin is 1.01x on
the worst of eight trials, at an identical frequency. A threshold test counting
`>=` was the wrong instrument to have predicted with, and quoting its verdict
without the 50.5 beside it would be the most misleading single line this run
could produce.

The alternative outcome flagged in advance -- a static floor far *below* the
measurement, indicating only that the simulator's 400 discrete scatterers are
cleaner than real clutter -- did not occur. The floor came out just under the
measurement, which is the outcome that carries information.

### Final statement for this run

**The validated chain does not detect micro-motion at Giza on this collect.** The
peak it reports is at the same frequency, and within 1% of the same prominence,
as the peak the identical processing produces from a simulated world containing
no motion at all.

The shuffle test passed it at p = 0.03. That is now a documented false positive
with a known mechanism, and the more useful result of this run than the null
itself: **`--shuffle-looks` and `--null-trials` are not valid nulls for a phase
observable under heavy overlap.** Shuffling inflates the median per-step phase
noise from 0.052 to 1.878 rad -- 36-fold -- so it destroys far more than the time
ordering it is supposed to isolate, and any red series beats the roughened
version of itself. `--null-static` is not fooled, because a motionless scene goes
through the same overlap, the same unwrap and the same detrend.

## B's static null

Run after the above was committed, because B had been left with no static floor:
its phase was saturated wrap in every window, and the reading recorded above was
that there was nothing for a floor to adjudicate. That reasoning is not good
enough. "Saturated" is a claim about *this* scene, and the way to test it is to
put a motionless scene through the identical chain and see whether it saturates
too.

```sh
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 256 --cell 2.0 \
    --subap pulse --overlap 0.9 --estimator phase \
    --n 2048 --win 32 --coherence 0.4 --null-static 8
```

### A FOURTH PREDICTION, recorded before the result

B's measurement is prominence 89.8 at 0.830 Hz, from a series pinned at lambda/2
in every window. The static scene is 400 point scatterers with Rayleigh
reflectivity over this collect's geometry -- discrete bright targets, no
distributed clutter, so its per-look SNR is probably *better* than the real
desert's.

- **P7** The static floor lands close to 89.8 -- within a factor of about 1.5 --
  and its trials report frequencies at the bottom of the band, as every red
  series here has. Whether any single trial reaches 89.8 is close to a coin
  flip and is not the interesting part.
- **P8** The interesting part is what it implies about saturation. At 0.16 s per
  look, phase noise per step is set by how much energy one look collects, and
  400 bright scatterers collect more than desert does. So if the static floor
  comes back *near* B's measurement, the wrap is a property of the chain at this
  look count and no target on this collect would do better. If it comes back
  well *below*, the wrap is a property of Khufu's backscatter, and a brighter
  structure might not saturate -- which would make B's failure specific rather
  than structural.

**What would falsify both:** a static floor far below 89.8 *and* trials reporting
frequencies spread across the 31 Hz band rather than bunched at the bottom. That
would mean the motionless chain produces a white spectrum, that B's red one came
from the scene, and that the saturation reading is wrong.

`--null-static` reports prominences only, not the per-trial shift series, so it
cannot measure the simulated scene's peak-to-peak directly. P8 is therefore
inferred from the floor rather than measured, and a direct test would need a
static collect written to disk and put through `mmotion --shifts` separately.

### Result: ABANDONED after one trial of eight

Stopped deliberately, not because it failed. P7 and P8 stand recorded and
unanswered, and this section exists so that a prediction with no result cannot
quietly disappear.

**Why it was stopped.** Two reasons, and the runtime is the weaker one.

1. **The test could not answer the question it was run for.** As noted when it
   was launched, `--null-static` reports prominences only, never the per-trial
   shift series, so it cannot measure whether the *simulated* scene saturates.
   P8 -- is the wrap a property of the chain or of Khufu's backscatter? -- would
   have been inferred from a floor rather than measured. Three hours for an
   inference that would still need hedging.
2. **B needs no null to be disqualified.** A null floor asks whether a
   measurement exceeds what a motionless world gives. B has no measurement to
   compare: median peak-to-peak displacement pinned at exactly lambda/2, median
   largest phase step 6.226 rad, in all 178 windows passing coherence. That is
   the arithmetic of a phase uniform on [-pi, pi], and it is already conclusive.

Reading Vattulainen et al. (2026) in the middle of the run also made B look like
the wrong question. Their nine ground-truth tests operate at **26-81
sub-apertures and 23-49% overlap with correlation tracking**. B's 2048 looks at
0.9 overlap with a phase estimator was a construction of this project's to reach
31 Hz, not a configuration anyone has validated. The successor run is
`runs/giza/2026-07-30-validated-spot-khufu`.

**The one trial that did complete**, recorded because it exists:

```
  static trial 1/8: prominence 0.0 at 0.031 Hz
```

**Do not read this as evidence against saturation.** Prominence 0.0 is what
`rs_null_static()` stores when `rs_spectrum_best_window()` returns
`RS_ERR_RANGE`, i.e. when *no window was eligible at all*. The simulated scene is
400 point scatterers spread over the grid extent, and `mmotion`'s own `--coherence`
help says isolated point targets on an otherwise empty scene score below 0.4 even
when tracking perfectly. The overwhelmingly likely reading is that the 0.4
coherence mask emptied the simulated scene, not that a motionless world produced
a clean spectrum. One trial, with a confound, is not a floor.

Anyone resuming this should either drop the mask to `--coherence 0` for the
static trials or give the simulator enough scatterers to look like clutter. Both
change what the floor means, which is itself worth thinking about before
quoting a number from it.

### Not done in this run

- A tomogram. Neither configuration produced a micro-motion result worth
  focusing into depth, and building one would have meant inverting a drift.

---

# Amendment, 2026-07-31: configuration A was outside its own band

`f_max` above is computed from the **step** between sub-apertures. That is the
wrong quantity. Each sub-aperture averages the target's motion over its own
duration `t_sap`, which is a lowpass whose first null is at `1/t_sap`, so the
series carries nothing above `1/(2*t_sap)` however finely overlap samples it.
Overlap buys resolution in time, not bandwidth.

| config | t_sap | f_max as quoted | band it can carry | overstated |
|---|---|---|---|---|
| A, uniform N=128 ov=0.99 | 14.480 s | 3.45 Hz | **0.0345 Hz** | 100x |
| B, pulse N=2048 ov=0.9 | 0.160 s | 31.47 Hz | 3.129 Hz | 10x |

**Every frequency A reported is above the band A could carry** -- 0.054 Hz by a
factor of 1.6, 0.324 Hz by 9.4, and 1.026 Hz by 30. The high overlap is what did
it: at fixed `n_looks`, `denom = N - (N-1)*overlap` collapses to 2.27 at 0.99,
so each sub-aperture spans 14.5 seconds of a 32.9 second dwell.

This does not overturn the run's conclusion, which was already a null. It
supplies the mechanism, and it explains **P5** directly: the peak follows
`--fmin` because there is no signal anywhere in the search region to compete
with the noise, the entire region being out of band.

B is inside its band at 0.324 Hz and is unaffected by this.

`resonarsat validate` now computes the band from `t_sap`, and `mmotion` warns at
run time when the peak it just printed came from beyond it.

---

# Amendment 2, 2026-07-31: configuration A's coherence mask was vacuous

`rs_splitband_shift()` estimates coherence as the mean magnitude over every pair
of sub-looks -- the standard estimator (ESA TM-19 Part C, Eq. 1.14) applied to a
stack repeat-pass interferometry never has: **overlapping** sub-apertures. Two
looks whose bands overlap share spectral content by construction, so their
pairwise coherence is high whatever the scene does. For a white scene it is
roughly the fraction of band they share, `1 - d*(1-overlap)` for looks `d` apart.

Averaged over all pairs that is a floor the estimator cannot report below:

| config | pairs sharing band | floor on an incoherent scene | mask used |
|---|---|---|---|
| A, N=128 ov=0.99 | 95% | **0.574** | 0.4 |
| B, N=2048 ov=0.9 | 1% | 0.004 | 0.4 |

**A's 0.4 mask sits below its own floor.** It passed every window it was given,
regardless of the scene, and the `--coherence 0.4` in the command above did no
filtering at all. B's identical mask is meaningful.

This is the same phenomenon as TM-19's coherence bias (Part C, Eq. 1.15,
`sqrt(pi/4N)` at true zero) but far larger and architecture-specific: the sample
count contributes about 0.03 at a 32x32 window, where band overlap contributes
0.57.

It compounds with the band error in Amendment 1. A ran with a sub-aperture too
long to carry any of the frequencies it reported, behind a mask that could not
reject anything. Neither was visible in the output.

`resonarsat validate` now computes this floor and calls a mask at or below it
vacuous.
