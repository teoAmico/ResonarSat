# Open follow-ups

Work identified and deliberately not folded into the change that found it.
Each entry says what the defect is, what evidence exists, and what would have
to be true for the fix to be correct.

Tracked, unlike most of `docs/`. The allow-list exists to keep development notes
and plans local, and this began as one. It is not one now: it is the record of
what has been tried and disproven, including entries that withdraw earlier
entries, and losing it means the next person repeats the work rather than reading
it.

Resolved entries are not kept here: once a finding is settled it belongs next to
the code it constrains, and the reasoning that produced it is in the commit that
made the change. Four were retired that way -- the correlation-transcendental identity,
the quantisation floor's derivation, the distributed-texture fixture questions,
and the sub-look coherence measurement behind the phase estimator's
non-accumulating design. They now live in `rs_microm_estimator_t`,
`rs_spectrum_best_window()` and `rs_microm_track()`.

Questions for the method's author are in
[`CORE-QUESTIONS.md`](CORE-QUESTIONS.md), with the detailed working list in
[`OPEN-QUESTIONS.md`](OPEN-QUESTIONS.md).

---

## 1. The quantisation floor is a per-window test with no multiplicity correction

`rs_spectrum_best_window()` applies the 2.449-step floor to each window
independently and returns the most prominent survivor. Nothing accounts for how
many windows were tested, so the chance that *some* window crosses the floor on
quantisation noise alone grows with the grid.

**Measured on the real Giza collect**, same scene and same chain at two grid
sizes:

| | 225 windows | 961 windows |
|---|---|---|
| excursion exactly zero | 225 (100%) | 958 (99.7%) |
| cleared the floor | 0 | 2 |
| reported | `RS_ERR_RANGE` | 0.183 Hz, prominence 29.9 |

The observable is the same in both. Only the number of opportunities changed,
and that flipped an honest refusal into a confident-looking frequency. The
winner sits at an observation ratio of exactly 3.00 -- one of the frequencies
where this observable's response is exactly zero -- and its series is a
five-level staircase: an artefact by every available test.

**Implemented so far:** the count is reported. `rs_spectrum_best_window()` takes
`out_n_candidates` and `mmotion` prints `N of M windows were eligible` beside
every selection, warning below four on a bound taken from the window geometry --
windows overlap at a stride of half their width, so a resolvable target occupies
a 2x2 block at minimum.

**Not implemented, and what a fix has to do.** Not merely raise the threshold;
that trades one arbitrary constant for another. The floor is three sigma for
*one* window. Testing `n` needs either a correction over the EFFECTIVE number of
independent windows -- fewer than `n`, because overlapping windows are
correlated -- or a different statistic: the *fraction* of windows clearing the
floor and its spatial contiguity is far more informative than the best single
one, since a real localised mode should occupy a patch rather than two scattered
windows.

---

## 2. No distributed-texture fixture exists on which the chain demonstrably works

The current one cannot work, for reasons now understood and recorded in
`rs_microm_estimator_t`: a lone vibrating point inside static clutter is diluted
roughly 9:1 in its own correlation window, and the paper front end's sweep puts
every resolvable bin at an integer observation ratio.

`sim_cphd --clutter-vib` was added so the whole patch can vibrate coherently,
which is what a structure's surface does. A configuration that actually recovers
an injected frequency on distributed texture -- coherent motion, a window
spanning several sub-look resolution cells, sub-looks that share most of their
spectrum, or the phase estimator -- has **not** been demonstrated.

**ANSWERED 2026-07-31, and the answer is negative for both estimators.**

Swept with `sim_cphd --clutter-vib` -- 400 coherently vibrating Rayleigh
scatterers, so the patch moves as a whole -- at the one operating point ever
demonstrated to track anything: pulse route, 128 looks, zero overlap, 0.5 m
cell, 32 px window, upsample 40. Displacement held constant at 10 px
peak-to-peak by scaling amplitude as 1/f, so amplitude is not confounded with
frequency. Five frequencies from 0.3 to 1.1 Hz, three clutter seeds, each with a
static control on the identical clutter.

```
                 0.3    0.5    0.7    0.9    1.1   static   slope    rms
correlation s7  2.671  2.117  0.353  0.806  1.109  2.974   -2.217  1.2935
correlation s11 1.159  0.504  1.109  2.772  2.369  0.554   +2.344  1.0973
correlation s23 0.806  1.462  0.302  1.411  2.319  0.806   +1.487  0.7857
phase       s7  2.369  1.008  1.058  2.319  1.462  1.008   -0.252  1.1672
phase       s11 1.764  1.764  1.764  1.764  1.764  1.663   +0.000  1.1010
phase       s23 0.857  2.974  1.411  1.815  0.302  1.109   -1.135  1.2970
```

Neither tracks. Correlation's slope swings from -2.2 to +2.3 across seeds, which
is noise rather than a trend; phase on seed 11 returns a fixed 1.764 Hz for every
injection with its static control at 1.663. Every static control lands inside the
spread of the moving reports. Doubling the displacement to 20 px changes nothing
-- phase then returns 1.109 Hz for all five injections AND for the static
control.

**The operating point was admissible, which is what makes this decisive.** The
tool reported sub-look resolution 8.26 m and an observable band to 3.23 Hz, so
the ambiguity ceiling is 24.8 px against a 4 px textured floor: the injected
10 px sits mid-window, and 20 px is still inside it. Every injection is well
inside the band, and the observation ratio at 0.5 Hz is 0.33. This is not a
misconfigured run.

**Note the near-misses.** Correlation seed 11 reports 0.504 against an injected
0.5, and seed 23 reports 0.302 against 0.3. Read one at a time those are
recoveries; read across the sweep they are the coincidences that make a
per-point criterion useless, and they are why `rs_track_fit()` exists.

**What this closes.** The synthetic evidence for the tracker now rests on
exactly one configuration -- correlation, isolated point target, empty scene,
128 looks, zero overlap -- which is the easiest case that exists and the
furthest from a real structure. Distributed texture was the last standing
explanation for the negatives that was not "the chain does not work". It no
longer stands.

**Not written as a regression test, deliberately.** The correlation sweep costs
2 minutes 48 seconds for three seeds, which would quadruple `ctest`. The result
is recorded here and in `IMPLEMENTATION-VERIFICATION.md` instead. If it is ever
worth locking in, one seed at three frequencies is about a minute and would
catch a change of sign.

---

## 3. The Capella SGN override is keyed on a vendor string

`rs_read_cphd()` inverts the FX-to-delay transform for Capella products because
their signal does not follow their declared `Global/SGN`; read to the standard,
their imagery comes back mirrored in range. The trigger is a substring match on
`CollectorName`.

That is blunt in both directions. If Capella corrects the metadata the override
keeps firing and silently starts mirroring; if another vendor ships the same
defect it is not caught.

A robust check would test the data rather than the label -- compressing a pulse
both ways and choosing the direction whose energy lands inside the declared
`TOA1`/`TOA2` support would decide it per product, without a name. Not
implemented.

---

## 4. Long dwells may need to be deliberately truncated, and never have been

Every run this project has made feeds the whole collect in: 32.869 s at Giza,
and the Istanbul candidate would be 25 or 60 s. The published work appears not
to do that.

**The direct evidence.** `MODIFIED-BACKPROJECTION.md` records that the Trento
corner-reflector validation reconstructs displacement from **the first three
seconds of a roughly 20 s dwell**, and that the authors note the spectrogram's
clarity degrades after about 11 s. It is already flagged there as "a design
choice this project has not tried."

**The arithmetic that independently implies the same number.** The tracker's
amplitude window reaches its 4.4x margin -- the only ratio ever demonstrated to
work -- at `t_sap` of roughly 0.16-0.21 s across every X-band sensor surveyed in
`DATASETS.md`, because the window depends on `t_sap` in seconds and not on the
aperture fraction. The published campaigns state `alpha` of 4.5-7.6%. Since
`alpha = t_sap / T`, those two together give

```
T = t_sap / alpha = 0.16/0.076 .. 0.21/0.045  =  2.1 .. 4.7 s
```

So the literature's *effective* dwell is a few seconds, arrived at from its
stated aperture fraction without reference to Trento, and Trento's three seconds
sits inside that range.

**Why this is not merely a curiosity.** If the usable dwell is ~3 s, then a
32.9 s collect is not a better version of a 3 s one -- it is eleven windows, and
processing it as a single span buys frequency resolution while risking whatever
degrades after 11 s. Candidates for what that is: orbital curvature the
backprojection grid does not model over the full arc, atmospheric phase drift,
target decorrelation, or the phase-unwrap accumulation `microm.h` already
predicts and which run A exhibited.

**What would settle it.** Cheap, and no new machinery: cut a collect into
consecutive ~3 s spans, run the identical configuration on each, and compare.
Three outcomes distinguish themselves:

- the reported frequency is stable across spans -- the long dwell was fine, and
  the truncation buys nothing but costs `df`
- it is stable early and wanders late -- something accumulates, and the point at
  which it starts is measurable
- it is unstable everywhere -- the configuration was never measuring anything,
  which is a different finding and worth having

The third is the one the null tests already suggest, and this is a cheaper way
to reach it than another full-dwell run. Scoring must use `rs_track_fit()` over
an injected sweep, not a single reported value, and must be pooled over seeds.

**Where it bites first.** The Istanbul 25 s collect is eight such spans. Deciding
before the run whether to process it whole or in spans is a configuration choice
that belongs in the `RUN.md` question, not a post-hoc reinterpretation.

---

## 5. `rs_slc_t.r0` means different things in different readers

`slc.h:65` documents `r0` as "slant range of first range sample, m". The UAVSAR
reader fills it that way. **The SICD reader does not**: `readers/sicd.c:410`
sets it from the product's `SCPCOA/SlantRange`, which is the range to the scene
centre point.

On a wide swath those differ by half a swath -- tens of metres on a spotlight
product, kilometres on a stripmap one.

**Found by** writing `validate --sicd`. The first version added
`0.5 * n_rg * rg_spacing_m` to reach the scene centre, which is right for the
documented meaning and double-counts for what the SICD reader actually stores.
The command now takes `r0` as-is with a comment saying why, so the defect is
visible rather than compensated for in one caller.

**Why it matters beyond that one call site.** Every consumer of an `rs_slc_t`
has to guess which convention it is holding, and the guess is invisible when
wrong: a slant range off by half a swath still produces a complete image and a
complete spectrum. `rs_geo_slant_to_ground()` and the sub-look resolution both
take `r0`, so the error propagates into geolocation and into the ambiguity
ceiling.

**What a fix has to do.** Not simply re-point the SICD reader at the first
sample, because `SCPCOA/SlantRange` is the better-conditioned quantity and is
what the product actually guarantees. Either:

- keep `r0` as the first-sample range everywhere and derive it in the SICD
  reader from the SCP range, the SCP pixel row and the range spacing, which the
  reader already parses; or
- redefine the field as the scene-centre range, fix UAVSAR to match, and rename
  it so no existing caller keeps its old assumption silently.

The second is cleaner and the rename is what makes it safe. Neither is done.
Until one is, `r0` should be read as "a slant range into this product, reader's
choice which".

---

## 6. The sub-look images are correct; the tracker does not read them

`POSITIVE-CONTROL.md` localised the failure to "what reaches the correlator --
patch extraction, the reference look, or the sub-look images themselves" and
said "that is the next thing to bisect and it is not done here." It is now done
one level further, and the sub-look images are eliminated.

**The images carry the motion.** Taking the plain argmax of `|look[i]|` over the
whole image, for an isolated target injected at 0.5 Hz and 20 mm (predicted
azimuth excursion +/- 8.4 cells at a 0.5 m cell):

| configuration | fitted amplitude | predicted | variance at f | phase |
|---|---|---|---|---|
| 128 looks, overlap 0.00 | 8.24 cells | 8.4 | **0.933** | 0.01 rad |
| 159 looks, overlap 0.88 | 10.65 cells | 10.5 | **0.832** | 0.01 rad |
| 64 looks, overlap 0.50 | 7.58 cells | 8.4 | **0.979** | 0.01 rad |

Right amplitude, right phase, most of the variance at the injected frequency --
and the failing configurations are as good as the working one. The peak
magnitude varies by only 1.16x across the series, so the target is not fading.
An integer argmax with no sub-pixel refinement at all recovers the motion.

**The tracker does not.** Handed the same stack, for the window containing the
target:

```
tracker shift series, first 10:  +0.0 -12.0 -15.2 +14.6 +10.3 +8.3 +7.0 +7.9 +8.7 +14.2
argmax trajectory,   first 10:   +7.0  +7.0  +2.0  -3.0  -3.0 -8.0 -8.0 -8.0 -3.0 -3.0
```

Correlation between them **-0.190**, with **4.1%** of the tracker's variance at
the injected frequency against 93% for the argmax. The tracker's excursions
reach +/-15 px where the target moves +/-8.4, which is the saturating-argmax
signature already recorded.

**Some window does track.** Scanning all windows, the best reaches correlation
**+0.832** and 74% of variance at f -- but it is not the window containing the
target, and which window wins moves with the target's position. So the
information survives into the tracking stage and is then attributed to the wrong
place.

**What this settles.** The fault is downstream of the sub-look images, in the
tracker's consumption of them. It is not missing physics, not a defect in
focusing or sub-aperture formation, and not a sensitivity limit -- the signal is
demonstrably present and recoverable by the crudest possible estimator. That
eliminates two of the three explanations `README.md` leaves open.

**What is NOT established.** Why. Two hypotheses were tried and neither holds:

- *Window-edge crossing.* The idea that the tracker reads the target entering and
  leaving the analysis window rather than its displacement. **Refuted:** with a
  64 px window and 32 px stride the target sits 24 px clear of both edges and a
  window still tracks at +0.861.
- *Excursion as a fraction of window size.* Suggested by the containing window
  failing at 32 px and a window tracking at 64 px. **Not measured** -- the sweep
  written for it mis-identified the containing window and returned nothing
  usable. It remains the most promising next cut.

The next step is to compare, for a single window and a single look pair, the
correlation surface the tracker computes against the shift the argmax implies.
That is one window and two images, and it is where the +0.832 and the -0.190
have to diverge.

### Item 6, resolved: the sub-looks are too decorrelated to correlate

Done as described above -- one window, the same patches `rs_microm_track()`
extracts, handed to the same primitive it calls.

**The primitive is sound.** Look 0 against itself returns shift +0.00 at peak
1.000. Where the correlation peak is high the shift is right: look 16 peaks at
0.571 and returns -10.25 px against an argmax truth of -10.0. Where the peak
collapses to 0.02-0.2, the returned shift is unrelated to the truth and reaches
the search extent.

**The peak collapses because consecutive sub-looks share no pulses.** Measured
mean correlation peak against look 0, over the window containing the target:

| overlap | mean peak | rms(coreg - truth) | looks agreeing within 2 px |
|---|---|---|---|
| 0.00 | 0.090 | 17.76 px | 10 of 128 |
| 0.50 | 0.107 | 16.07 px | 19 of 128 |
| 0.75 | 0.177 | 13.49 px | 28 of 128 |
| 0.90 | 0.310 | 10.29 px | 40 of 128 |
| 0.95 | 0.340 | 8.30 px | 27 of 128 |

`rs_microm_estimator_t` already carries the table this reproduces -- coherence is
very nearly the fraction of pulses two sub-looks have in common, reaching 0.07 at
zero sharing. The measured 0.090 at zero overlap is that number. **The
configuration this project calls its working operating point -- 128 looks, zero
overlap -- correlates sub-looks that share no pulses at all.** Its own header
predicts a coherence of 0.07 for that, and the tracker is being asked to find a
displacement between independent speckle realisations.

**And raising overlap does not rescue it.** At 0.95 the peak only reaches 0.340
and 27 of 128 looks agree. Meanwhile high overlap sharpens the sub-look and
collapses the wrap ceiling -- 0.88 overlap gives a 4.8 px ceiling against a 7 px
floor, measured, so the window closes.

**The two constraints oppose each other on the same knob, and that is new.**
Coherence needs overlap; the ambiguity ceiling forbids it. Both are documented
in this codebase, separately, and neither is stated as the other's opposite.
Together they explain why no configuration has worked: there may be no overlap
that satisfies both, and if so the correlation estimator cannot work on this
geometry at any setting rather than merely failing at the ones tried.

**What would settle that.** Compute both curves over overlap on the same
collect -- coherence from pulse sharing, ceiling from sub-look resolution -- and
find whether they cross. That is arithmetic on quantities `validate` already
derives, needs no processing run, and either produces an admissible overlap band
or proves there is none. It is the single most valuable thing left.

### The curves cross, and the binding constraint is the reference look

Computed on the Giza geometry as `validate` derives it -- T 32.869 s, lambda
0.0322 m, R 762.8 km, V 7263 m/s, 0.4 m cell, 7 px floor. Ceiling from sub-look
resolution, coherence from pulse sharing via the table in
`rs_microm_estimator_t`.

**There is an admissible band.** Seventeen (N, overlap) combinations clear both
constraints -- ambiguity ceiling above the floor AND adjacent-look coherence
above 0.6:

| N | overlap | t_sap | ceiling | gamma | band |
|---|---|---|---|---|---|
| 256 | 0.75 | 0.508 s | 12.5 px | 0.61 | 0.98 Hz |
| 512 | 0.90 | 0.631 s | 10.1 px | 0.78 | 0.79 Hz |
| 1024 | 0.90 | 0.318 s | 19.9 px | 0.78 | 1.57 Hz |
| 2048 | 0.75 | 0.064 s | 98.9 px | 0.61 | 7.80 Hz |

So the two constraints do NOT exclude each other, and the earlier suspicion that
they might is withdrawn. High look counts at moderate-to-high overlap satisfy
both, with margins up to 14x the floor and bands out to 7.8 Hz.

**Every one of them is unreachable with a fixed reference.** The last column of
the full sweep is the fraction of the series that shares any pulses at all with
look 0, and across all seventeen it runs **0% to 3%**. Looks i and j overlap only
while `|i-j| < 1/(1-overlap)`, so with a fixed reference the coherent span is a
handful of looks out of hundreds and the rest is independent speckle. That is
what the measured peak of 0.090 at zero overlap was showing, and raising overlap
does not fix it because the span grows as `1/(1-overlap)` while the series grows
as N.

**So the defect is the reference scheme, not the geometry.** `--reference first`
is the default and it discards the coherence the admissible band exists to
provide.

**And the estimator built for exactly this has never been tested.**
`rs_microm_estimator_t` describes `RS_MICROM_EST_SPLITBAND` as using "all N^2
interferograms rather than the N-1 formed against one reference", coming "within
0.5 dB of the Cramer-Rao bound", and requiring "interferometric coherence between
looks; it has nothing to track if that is absent". That requirement is precisely
what the admissible band supplies (gamma 0.61-0.85) and what zero overlap denies
(gamma 0.07). Correlation and phase have now both been swept and both fail.
Split-band phase linking is the untested third, and it is the one whose stated
design matches the constraint that turns out to bind.

**Next.** Sweep `RS_MICROM_EST_SPLITBAND` with `rs_track_fit()` at N=512,
overlap 0.90 -- ceiling 10.1 px, gamma 0.78, band 0.79 Hz, all comfortably
admissible -- and again at N=1024, overlap 0.90 for a 1.57 Hz band. Static
control on the same clutter, pooled over seeds. If it tracks, the negatives so
far are a reference-scheme defect and the method is recoverable. If it does not,
all three estimators have failed at operating points this project's own
arithmetic calls admissible, and that is the end of the line.

---

## 7. RS_MICROM_REF_LAG implemented; it fixes coherence and exposes a frequency floor

A fourth reference mode: each look against the one `ref_lag` places before it,
with no accumulation. `--reference lag --lag N`. Twenty-odd lines beside the
other three, plus `ref_lag` in `rs_microm_params_t`.

**It does what it was designed to do.** On the distributed fixture at 512 looks
and 0.75 overlap, tracking quality goes from **0.049 to 0.806**. The coherence
the fixed reference was throwing away is recovered exactly as predicted.

**And it is not sufficient, for a reason that is now arithmetic.** Three
constraints have to hold at once:

```
(1) ceiling      1.5 * res_sap / cell  >  2A        excursion fits, no wrap
(2) coherence    gamma(1 - L(1-overlap)) > 0.6      patches correlate at the lag
(3) differential 2A * sin(pi f L dt)   >  7 px      the difference clears the floor
```

(3) is the new one and it is what a differencing observable costs. At lag 1 on
the fixture tested, the differential is `2 x 8.4 x sin(pi x 0.5 x 0.0388)` =
**1.0 px** against a 7 px floor -- the signal is real, coherent and far too small
to measure. Raising the lag raises the differential and destroys (2).

**Swept over N, overlap, lag and frequency on the Giza geometry, 15 combinations
satisfy all three -- and every one of them starts above 1.5 Hz:**

| f_min | f_max | N | overlap | lag | gamma | ceiling |
|---|---|---|---|---|---|---|
| 1.56 | 15.59 Hz | 4096 | 0.75 | 1 | 0.61 | 197.7 px |
| 1.56 | 7.80 Hz | 2048 | 0.75 | 1 | 0.61 | 98.9 px |
| 1.96 | 6.24 Hz | 2048 | 0.80 | 1 | 0.67 | 79.2 px |
| 2.60 | 9.36 Hz | 4096 | 0.85 | 1 | 0.72 | 118.7 px |

**So the differencing observable has a hard low-frequency floor near 1.5 Hz on
this collect**, because the differential goes as `f` at low `f` while the lag is
capped by coherence. Every measurement this project has attempted has been at
0.05-1 Hz, which is exactly where this observable has no sensitivity, and the
absolute reference that does have low-frequency sensitivity is the one that
decorrelates. That is the trade, stated as an inequality rather than a suspicion.

**The constructive part.** 1.56-15.6 Hz is not a useless band -- it contains the
1-3 Hz where medium-span bridge decks live, which is the target class the archive
survey in `DATASETS.md` settled on. The Istanbul and Bratislava candidates were
chosen for structures whose modes sit inside the one band this observable can
reach.

**Untested, and this is the next step.** The table above is arithmetic on
constraints, not a measurement. A confirming sweep needs N=2048 at 0.75 overlap,
lag 1, injected 2-4 Hz at roughly 15 mm, scored with `rs_track_fit()` and pooled
over seeds. At about four minutes a run it is roughly twenty minutes for five
frequencies on one seed. Nothing below 1.5 Hz should be expected to work, and a
sweep that includes sub-Hz points will show the floor rather than a failure.

### Split-band, for completeness

`RS_MICROM_EST_SPLITBAND` was swept at N=512/0.90, N=1024/0.90 and N=2048/0.90.
All three return one fixed frequency for every injection and for the static
control: 0.417 Hz, 0.391 Hz and 0.391 Hz respectively, slope 0.000, identical
across three clutter seeds. So all three estimators fail with the FIRST
reference, which is consistent with the reference being the defect rather than
the estimator -- and is why the lag mode was written.

### Item 7, first measurement: the lag mode tracks above its predicted floor

N=2048, overlap 0.75, lag 1, excursion held at 44 cells, one clutter seed.

```
f=2.5 Hz  differential  6.7 px (BELOW the 7 px floor)  ->  0.391 Hz
f=3.0 Hz  differential  8.1 px                         ->  2.604 Hz
f=3.5 Hz  differential  9.4 px                         ->  3.125 Hz
f=4.0 Hz  differential 10.8 px                         ->  3.581 Hz
f=4.5 Hz  differential 12.1 px                         ->  4.492 Hz
STATIC    no motion                                    ->  0.391 Hz
```

Over the four points above the floor: **slope +1.224, rms 0.3439 Hz**.

**This is the first configuration in the project whose output follows its
input.** Every other estimator and operating point tested -- correlation, phase,
split-band, at every setting -- returned slope 0.000 and a single fixed value.

Three things make it more than a fluke:

- **The static control returns 0.391 Hz**, well clear of every supra-floor
  reading. The rising sequence is not something the processing produces
  regardless of the scene.
- **The one sub-floor point returns that same 0.391 Hz.** The only injection
  whose differential falls below the tracking floor is the only one that reports
  the artefact, and it reports exactly the artefact the static scene does. The
  split lands where the arithmetic put it.
- 0.391 Hz is also what split-band returned at N=512, 1024 and 2048, so it is a
  property of the chain rather than of this mode.

**It does not pass `rs_track_fit()`.** The rms is fourteen times the half-bin
bound, there is a systematic 0.4 Hz low bias on three of four points, and the
slope is 1.224 rather than 1.0. **And it is one seed**, which by this project's
own standard settles nothing about a configuration.

So: tracks, imprecisely, above the frequency floor its own arithmetic predicts,
on one realisation. Much more than "fails", much less than "works".

**Next, in order.** Repeat over seeds 11 and 23 and pool. Then chase the low
bias -- a 0.4 Hz offset that vanishes at 4.5 Hz is not obviously a scale error,
and the differencing response |2 sin(pi f L dt)| is not flat across the sweep, so
compensating for it before peak-picking is the first thing to try.

### Item 7, pooled over three seeds: the trend survives, the precision does not

Same configuration, seeds 7, 11 and 23.

| f | seed 7 | seed 11 | seed 23 | spread |
|---|---|---|---|---|
| 2.5 Hz (sub-floor) | 0.391 | 2.083 | 2.083 | 1.692 |
| 3.0 Hz | 2.604 | 2.604 | 2.604 | **0.000** |
| 3.5 Hz | 3.125 | 3.516 | 3.125 | 0.391 |
| 4.0 Hz | 3.581 | 3.581 | 3.581 | **0.000** |
| 4.5 Hz | 4.492 | 4.492 | 4.102 | 0.390 |
| static | 0.391 | 0.391 | **12.174** | -- |

**Pooled over the supra-floor points: slope +1.120, rms 0.3461 Hz, n = 12.**
Per seed the slopes are +1.836, +1.159 and +1.003 -- all positive, all near 1,
none collapsing or changing sign. No other estimator or operating point in this
project has produced a non-zero slope at all, so the trend is real and
reproducible.

**Four things qualify it, and the third is the one that matters.**

*The rms is 14x the half-bin bound.* It does not pass `rs_track_fit()`. What has
been demonstrated is a response to the injected frequency, not a measurement of
it.

*The static control is not stable.* 0.391, 0.391, 12.174 Hz. The last is near
this configuration's 12.8 Hz Nyquist, which is what a high-pass observable does
with noise. It never lands inside the 2-4.5 Hz band the injections occupy, so it
still separates -- but "the static control returns 0.391" was a seed-7 fact, not
a property.

*There is a systematic offset of about 0.39 Hz, and it is the chain's own
artefact frequency.* The spectrum bins here are 0.0651 Hz, and every reported
value is an exact bin: 32, 40, 48, 54, 55, 63, 69. The injections fall at bins
38.4, 46.1, 53.8, 61.4, 69.1. The reported bin is consistently about six below
the injected one, and 6 x 0.0651 = 0.391 Hz -- the same 0.391 Hz that the static
scenes, the sub-floor injection and every split-band run all return. That is
unlikely to be coincidence and it is a concrete lead.

**A bias survives redundancy**, so this must be understood before a network is
built on top of it. If the offset is the artefact mixing into the peak
selection, a bigger estimator inherits it intact.

*The pooling is weaker than it looks.* `--clutter-vib` moves every scatterer
coherently with the target, so changing the seed changes the speckle but not the
motion. Three seeds sharing an identical signal is a weaker control than three
independent realisations of a lone mover would be, and the exact 0.000 spreads at
3.0 and 4.0 Hz are the visible sign of it. Worth repeating against
`--clutter` without `--clutter-vib` before the result is leaned on.

### Item 7, corrected: the floor claim was computed on the wrong dt

Dumping the shift series with `--shifts` shows the tool's own layout, and it is
not what the arithmetic above assumed:

```
sub-apertures: 2048 looks, dt 0.0075 s
spectra: 1025 bins, 0.0651 Hz resolution
```

`dt` is 0.0075 s, not the 0.00975 s that `t_sap*(1-overlap)` with
`t_sap = T/denom` predicts. The series therefore spans **15.36 s of a 20 s
dwell** -- `1/0.0651` -- leaving 23% of the collect unused, and `df` is 0.0651 Hz
rather than 0.05. Whether the layout is meant to leave that tail is a separate
question and is not answered here.

Recomputed on the tool's dt, the differentials are:

| f | differential | vs 7 px floor |
|---|---|---|
| 2.5 Hz | 5.2 px | below |
| 3.0 Hz | 6.2 px | **below** |
| 3.5 Hz | 7.2 px | above |
| 4.0 Hz | 8.3 px | above |
| 4.5 Hz | 9.3 px | above |

So **two** injections were below the floor, not one -- and the 3.0 Hz point
reported 2.604 Hz on all three seeds, which is a tracking-like response rather
than the artefact.

**The claim that "the split lands where the arithmetic put it" is withdrawn.**
It rested on a dt this code does not use. The pooled slope of +1.120 stands, and
so does the fact that no other configuration produces a non-zero slope; the
supporting story about the floor does not.

### A provenance bug, introduced with the lag mode and now fixed

`--shifts` wrote `reference=first` for a run made with `--reference lag`, and the
`.meta` sidecar did the same: both label sites carried three-way conditionals
that fell through to "first". Every lag run so far is misrecorded in its own
metadata. Fixed at `main.c:1469` and `main.c:2205`.

Worth noting how it was caught -- not by a test, but by reading a dump while
chasing something else. The same three-way fallthrough pattern would silently
mislabel any future mode, and nothing in the suite checks that a run's recorded
provenance matches what produced it.

### Item 7, the bias explained: a 0.391 Hz multiplicative modulation

Dumped the shift series with `--shifts` and took its spectrum independently in
Python. **The tool's peak-picking is correct**: my spectrum of window 42 peaks at
2.6042 Hz, exactly what `mmotion` reported. The spectrum stage is faithfully
reporting the strongest bin of the series it is given, so nothing downstream of
the tracker is at fault.

The series is what carries the defect. Its strongest six bins, for a 3.000 Hz
injection:

```
2.604   0.391   2.995   5.599   3.385   6.250
```

- **2.995 Hz is the injected carrier, and it IS there** -- third strongest.
- **0.391 Hz is the artefact**, second strongest.
- **2.604 = 2.995 - 0.391**, the lower sideband, and it is the strongest bin.
- **3.385 = 2.995 + 0.391**, the upper sideband, fifth strongest.

Both sidebands present, symmetric about the carrier, spaced by the artefact
frequency. That is multiplicative modulation: the tracked series is being
multiplied by something oscillating at 0.391 Hz, and peak selection then picks a
sideband instead of the carrier.

**This explains the bias completely**, including why it is inconsistent. The
reported value is whichever of carrier and lower sideband happens to be stronger,
so most injections come back about 6 bins low and the 4.5 Hz case -- where the
carrier won -- came back exact.

**And it is good news for the method.** The injected frequency is present in the
tracked series. The tracker is recovering the motion; a modulation is stealing
the peak. That is a different and far more tractable problem than "the chain does
not measure anything", which is where this investigation started.

**What 0.391 Hz is remains unidentified.** It is exactly bin 6 of the 0.0651 Hz
axis, a period of 2.558 s fitting 6.0 times into the 15.36 s record, and it is
the same value the static scenes, the sub-floor injections and every split-band
run return. A modulation whose period divides the record length exactly is more
suggestive of the sub-aperture layout than of the scene. The layout here is 12
pulses per look stepping 3 pulses, which is exact and leaves 23% of the dwell
unused.

**Next:** demodulate, or find the source. If the modulation can be identified and
removed, the carrier is already there to be read -- and that is worth more than
the network, because it fixes the answer rather than averaging it.

### Item 7: the multiplicative-modulation explanation is WITHDRAWN

The static scene was dumped and the 0.391 Hz component measured in every window.
If it were a modulation imposed by the processing it would have the same phase
everywhere. It does not:

```
0.391 Hz component across 49 windows, static scene
  amplitude   min 0.0009   median 0.0046   max 0.0141 px
  phase       spread 5.09 rad over a possible 6.28
  phase concentration |R| = 0.174        (1.0 would be identical everywhere)
```

Random phase, and an amplitude of a few thousandths of a pixel. Folding a
window's series at the 341-look period shows no coherent waveform. **There is no
0.391 Hz modulation.** It is where the noise happens to peak in whichever window
wins the prominence contest.

The sideband reading was over-read from a top-six list. Its own evidence was
already against it: multiplicative modulation produces symmetric sidebands, and
the observed pair was not symmetric -- 2.604 Hz was the strongest bin in the
spectrum while 3.385 Hz was only fifth. Picking six bins out of a thousand and
finding two near `f +/- f0` is not the coincidence it looked like.

**What survives, and it is the important part:**

- The spectrum stage and peak-picking are correct. An independent spectrum of
  the dumped series matched the tool's reported 2.6042 Hz exactly.
- **The injected carrier is present in the tracked series** -- 2.995 Hz for a
  3.000 Hz injection, third strongest bin. The tracker recovers the motion.
- The reported peak is a different bin, consistently about 6 bins low across
  three injections and exact on the fourth. That the offsets cluster near 6 is
  still unexplained; it is not a modulation, and it is not spectral leakage or
  the differencing response, both of which bias the other way.

So the problem is **selection, not recovery**: the right frequency is in the
series and is not the largest thing in it. That is a much narrower problem than
where this investigation started, and it points at the peak-picking policy rather
than at the tracker.

**A concrete thing to try**, in preference to another hypothesis: report the
strongest few bins rather than one, and see whether the carrier is reliably among
them across a sweep. If it is, the recovery is real and the failure is entirely
in `rs_spectrum_best_window()` choosing between candidates it has no basis to
rank. That is measurable with the series already dumped, at no processing cost.

### Item 7 resolved: the recovery is real; the SELECTION POLICY discards it

For a 3.000 Hz injection, with the lag reference, the carrier is in the **top ten
bins of all 49 windows** -- rank 1 in 23 of them, median rank 2, worst rank 8.
The tracker recovers the injected frequency essentially everywhere.

The rank-1 frequency of each window, tallied:

```
  2.995 Hz   23 windows      <- the injection
  2.604 Hz   16 windows
  0.391 Hz    8 windows
  6.380 Hz    1 window
 15.169 Hz    1 window
```

**A plurality vote returns 2.995 Hz, within half a bin of the truth.** The tool
returned 2.604 Hz, because `rs_spectrum_best_window()` reports the rank-1 bin of
the single most prominent window and that window was one of the sixteen.

So the ~6-bin bias, the rms of 0.35 Hz and the "does not pass `rs_track_fit()`"
verdict are all artefacts of the selection policy, not of the measurement. The
measurement is right and is being thrown away at the last step.

**This is item 1 of this file, arriving from the other direction.** That item
already argued the case on multiplicity grounds -- "the *fraction* of windows
clearing the floor and its spatial contiguity is far more informative than the
best single one, since a real localised mode should occupy a patch rather than
two scattered windows" -- and proposed exactly this replacement. It was never
implemented. This is the measurement that shows what it costs: the difference
between reporting the injected frequency and reporting a neighbouring bin.

**What to implement.** A consensus statistic beside
`rs_spectrum_best_window()`: for each candidate bin, the number of windows
ranking it first (or within the top few), reported with that count so a caller
can see whether a peak is a consensus or a single window's opinion. Spatial
contiguity of the voting windows is the natural refinement and is what item 1
asks for.

**Before implementing, verify it generalises.** This is one frequency on one
seed. The dumped series for 3.5, 4.0 and 4.5 Hz would settle whether the
plurality lands on the carrier every time or whether 3.0 Hz was lucky -- and the
static scene must be checked too, since a consensus that also produces a
confident plurality on a motionless scene would be worse than what it replaces.

### The consensus statistic has a null behaviour, and it is a good one

The static control was the test that mattered, because a consensus that also
agrees on a motionless scene would be worse than the policy it replaces. It does
not agree:

| | distinct rank-1 winners | plurality | share | runner-up |
|---|---|---|---|---|
| **static, no motion** | **19** of 49 windows | 0.391 Hz | 29% | 12.174 Hz at 27% |
| **injected 3.000 Hz** | **5** | 2.995 Hz | 47% | 2.604 Hz at 33% |

Two things separate them, and the first is the sharper:

- **The number of distinct winners: 19 against 5.** With no signal the windows
  disagree about what the peak is; with signal they concentrate.
- **The plurality's margin.** The static scene's leader beats its runner-up by
  14 to 13 -- a tie, which is what noise looks like. The moving scene's leads
  23 to 16.

So this is not merely a better point estimate. It is a detection statistic with a
measurable null: agreement when there is signal, fragmentation when there is not.
That is the property `rs_spectrum_best_window()` has never had -- a single
window's argmax is equally confident either way, which is why prominence turned
out to be anti-correlated with correctness.

**Report the winner count alongside the frequency.** Nineteen distinct winners
over 49 windows should read as "no consensus" however prominent the leader is.

---

## 8. WITHDRAWN: "the defect is the reference scheme"

Same scene, same configuration, 3.000 Hz injected, the two references compared
by cross-window vote:

| reference | plurality | share | windows voting for the carrier | distinct winners |
|---|---|---|---|---|
| `lag` | 2.995 Hz | 47% | 23 of 49 | 5 |
| **`first`** (default) | 2.995 Hz | **61%** | **30 of 49** | 11 |

**`FIRST` recovers the injected frequency in MORE windows than `LAG` does.** The
reference this file spent item 6 and item 7 indicting is not the binding
constraint, and the mode built to replace it is worse on the measure that
matters.

**What was right.** The decorrelation is real and measured: look 0 against the
rest averages a correlation peak of 0.090 at zero overlap and 0.310 at 0.90,
while adjacent pairs reach 0.913. The sub-look images are correct. The tracker's
series in any ONE window correlates poorly with the truth.

**What was wrong.** Inferring from those that the reference was what stopped the
chain measuring. It does not: with the default reference the carrier is the
top bin in 30 of 49 windows. The chain was recovering the frequency the whole
time, in most of the scene, and `rs_spectrum_best_window()` was reporting one
window's argmax.

**This is item 1 of this file, and it has now been the answer twice.** Item 1
argued on multiplicity grounds that a single window's peak cannot be a
measurement and proposed a consensus statistic. Item 7 arrived at the same place
from the lag experiment. Neither was implemented, and every subsequent
conclusion in items 6 through 7 -- including a new reference mode, now committed
-- was built on top of the unfixed defect.

**The lesson worth keeping.** Three separate estimator "failures" (correlation,
phase, split-band), a reference diagnosis, a coherence analysis, a frequency
floor derivation and a new code path all followed from measurements taken
through a selection policy that cannot express disagreement. The policy should
have been fixed first. It is the cheapest thing in this file and it has been open
the longest.

**`LAG` should not be removed** -- it is documented, tested and harmless, and the
coherence reasoning behind it stands on its own -- but it should not be presented
as a fix for anything until it beats `FIRST` on a consensus measure, which it
currently does not.

---

## 9. Item 1 partly implemented: consensus and contiguity

`rs_spectrum_consensus()` now sits beside `rs_spectrum_best_window()`, returning
the frequency the most windows agree on together with how many agree, how many
distinct answers there are, and **the largest 4-connected block of agreeing
windows**. `mmotion` prints all four and warns on two conditions.

Measured on the same scene with and without motion, 3.000 Hz injected:

| | best_window says | agreement | distinct | largest block |
|---|---|---|---|---|
| moving | 2.995 Hz, prominence 45.5 | 23/49 (47%) | 18 | **15** |
| **static** | **12.565 Hz, prominence 23.1** | 8/49 (16%) | 33 | **3** |

On the motionless scene `rs_spectrum_best_window()` reports 12.565 Hz at
prominence 23.1 with nothing to mark it as noise. The consensus line reports 16%
agreement over 33 distinct answers in a largest block of 3, and both warnings
fire.

**Contiguity is the sharper of the two statistics.** Agreement separates the
cases 47% to 16%; the largest block separates them 15 to 3, and the static case
falls below a bound that comes from the window geometry rather than from tuning
-- overlapping windows put a resolvable target in a 2x2 block at minimum, so a
largest block under four cannot be a spatially resolved mode. That is the bound
`rs_spectrum_best_window()`'s own header already derives for its candidate
count; it applies unchanged here and is not a constant anyone chose.

**What is implemented and what is not.** Item 1 asked for either a multiplicity
correction over the effective number of independent windows, or ranking on the
qualifying fraction and its spatial contiguity. The second is now available as a
statistic a caller can read. Nothing RANKS on it yet: `rs_spectrum_best_window()`
still selects by prominence and is still what the tool reports as its answer, with
the consensus printed beside it. Replacing the selection rather than annotating
it is the remaining half, and it should wait for more evidence than three
detections on one seed.

**The threshold evidence remains thin** and is stated where the thresholds live.
`rs_spectrum_consensus()` applies none. `mmotion` warns below one third
agreement, from measurements putting correct recoveries at 47-61% and everything
wrong at 14-29%, and below a contiguous block of four, which is geometric.

23/23 pass in Release and under ASAN.

---

## 10. Re-scoring the documented results through the consensus gate

Step 3 of the plan in item 9. Every case below is one `IMPLEMENTATION-VERIFICATION.md`
already records, re-run through the agreement gate.

**The documented false positives mostly become declines:**

| recorded | re-scored | agreement |
|---|---|---|
| 0.90 Hz -> 1.811 Hz, "false second harmonic" | **NO FREQUENCY** | 22% |
| 1.10 Hz -> 2.213 Hz, "false second harmonic" | reports 0.102 Hz | **33%** |
| pair, 0.5 Hz -> 0.100 Hz | **NO FREQUENCY** | 30% |
| pair, 1.0 Hz -> 3.300 Hz | **NO FREQUENCY** | 32% |

**And the documented recoveries survive:** 0.3, 0.5, 0.7 and 0.9 Hz at the
working operating point all report, at 75%, 80%, 67% and 57% agreement, with the
consensus frequency correct in every case.

So the gate keeps what worked and refuses three of four things that did not.

**Two honest limits.**

*The 1.1 Hz case sits exactly on the threshold* -- 3 of 9 windows is one third,
the gate does not fire, and it reports 0.102 Hz against a 1.1 Hz injection. A
case landing precisely on a tuned constant is where that constant is guaranteed
to be arbitrary.

*The margin is thinner than item 9 claimed.* Correct cases run 57-80% and false
positives 11-33%: a gap, but the pair cases clear the threshold by two points
rather than comfortably. A threshold of 0.25 would pass both.

**What this means for the verification document.** Several entries reading "the
chain reports a wrong frequency" are more accurately "the chain reports a
frequency that a consensus read refuses". That is a materially different claim
and the document should distinguish them. It does NOT mean the chain now works:
the frequencies are still not recovered, only the false confidence is removed.

**Not yet re-scored:** the phase estimator sweep, the split-band sweeps, and the
distributed-texture results of item 2. Those were all run before the consensus
existed and all reported single-window answers.

---

## 11. The consensus gate is BLIND to common-mode artefacts

Re-scoring the phase estimator found the case that bounds everything in items 9
and 10. At `test_tracking`'s operating point, isolated point target:

```
inj 0.2 Hz -> 0.407 Hz   9 of 9 windows agree (100%)
inj 0.3 Hz -> 0.407 Hz   100%
inj 0.4 Hz -> 0.407 Hz   100%
inj 0.5 Hz -> 0.407 Hz   100%
inj 0.6 Hz -> 0.407 Hz   100%
inj 0.7 Hz -> 0.407 Hz   100%
STATIC     -> 0.407 Hz   100%
```

**Unanimous agreement on a pure artefact, including on a scene with no motion in
it.** The gate passes every one of these, at the highest confidence the statistic
can express.

**This is structural, not a threshold problem.** Agreement detects noise that is
INDEPENDENT across windows -- correlation's scattered peaks, speckle-driven
static answers, everything item 10 successfully refused. An artefact produced by
the PROCESSING rather than by the scene appears identically in every window, so
the windows agree about it unanimously. No value of the constant helps: 100% is
the ceiling.

**So the gate catches one class of false positive and passes another.** Item 9's
framing -- "a fragmented vote and a motionless scene look alike" -- is true only
for scene-driven noise. For common-mode noise a motionless scene looks like a
perfect detection.

**The only thing that catches this is a null control**, because a common-mode
artefact is by definition identical whether or not the scene moves, so the sole
way to see it is to run the identical processing over a scene known to be
motionless and compare. That is `--null-static`, and it is what `README.md`
already names as the credibility check that matters. Today's work does not
replace it and must not be described as though it does.

**A run needs both.** Agreement for scattered noise, a null control for coherent
noise. Neither substitutes for the other, and the phase estimator is the case
that proves it: item 10's gate refused three of four correlation false positives
and would pass every phase-estimator result ever produced.

**Where this leaves item 9's claim.** The consensus statistic is worth having and
its null behaviour is real for the failure mode it addresses. But it is a partial
check, and the sentence in `rs_spectrum_consensus()`'s header about a fragmented
vote and a motionless scene looking alike needs the qualification that it holds
for scene-driven noise only.
