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
