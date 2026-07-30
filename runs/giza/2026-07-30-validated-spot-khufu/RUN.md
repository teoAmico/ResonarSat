# Run: 2026-07-30 giza / validated-spot-khufu

**Question this run is meant to answer:** Does the operating point that
Vattulainen et al. (2026) validated against LVDT ground truth -- a few dozen
sub-apertures at moderate overlap, correlation tracking -- recover anything at
Khufu?

- git commit: `51bfcd6`
- started:    2026-07-30T16:52:13Z
- host:       Darwin arm64 (8 cores, 25.8 GB RAM), OpenMP enabled

## Why this run exists

Every micro-motion run this project has made was outside the envelope anyone has
validated, and it took reading the metrological assessment to notice.

Vattulainen et al., *Assessment of Spaceborne SAR Micro-Motion Measurement for
Vibration-Based SHM*, IEEE Access 14 (2026) 6043, measure a corner reflector on a
programmable shaker with a synchronous LVDT, over nine Umbra spotlight
acquisitions. Their Table 3 gives the processing parameters for all nine:

| | their range | this project's runs so far |
|---|---|---|
| sub-apertures `N` | **26 – 81** | 128, 512, 2048 |
| overlap `Omega` | **23 – 49%** | 40% (defaults), 90%, 99% |
| observation ratio `eta` | **0.39 – 0.69**, always below 1 | 0.78 and 31.25 in the phase runs |
| estimator | **correlation (SPOT)** | correlation under the patent chain; phase in the uniform runs |

`README.md` asserted the opposite -- that "the validated work uses thousands of
sub-apertures at 90-99% overlap rather than a few dozen" -- and that line is what
set `--overlap 0.99` in `2026-07-30-uniform-phase-khufu`. The paper's own
reasoning (Section VI-A4) is that three effects compete: Nyquist sampling wants a
short `t_sap`, cross-correlation precision wants a long one for finer azimuth
resolution, and motion-induced blurring puts a ceiling on it. "The combined
result of all three effects is a moderate overlap value and an optimal `f_s` near
the Nyquist rate."

So this run is the first at an operating point somebody has checked against
ground truth.

## Collect

```
CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd  ->  data/giza.cphd
  Capella open data, CC BY 4.0
  39,180,270,944 bytes (verified)
  dwell 32.869 s   335,141 pulses x 25,073 range bins   PRF 10,196.35 Hz
  carrier 9.3000 GHz   wavelength 0.0322 m
```

Same collect and same grid as every previous Giza run: `--at
29.979175,31.134186` -> `--offset -152,-552`, Khufu.

## Deriving the parameters -- and an unplanned check of the implementation

Targeting `eta ~ 0.6` for a 2 Hz vibration gives `t_sap = 0.3 s`. At
`Omega = 0.4`, `N = 182`. ResonarSat's uniform band layout
(`rs_subaperture_split()`, `denom = N - (N-1)*Omega`) then gives:

```
denom  109.6000     t_sap 0.2999 s     dt 0.1799 s
f_s = 1/dt = 5.557 Hz          f_max = 1/(2 dt) = 2.779 Hz
```

The paper states the same layout in different variables -- aperture fraction
`alpha = t_sap/t_a`, then `N = (alpha^-1 - Omega)/(1 - Omega)` (Eq. 5) and
`dt = (t_a - t_sap)/(N - 1)` (Eq. 6). Substituting `t_sap = 0.2999`:

```
alpha 0.009124   ->   Eq. 5 gives N  = 182.0000    (chosen: 182)
                      Eq. 6 gives dt = 0.1799 s    (layout: 0.1799 s)
```

**Both reproduce exactly.** That was not the point of the run, but it is the
first independent confirmation that this project's uniform sub-aperture layout
agrees with a published formulation of the same decomposition, to four decimal
places, in a paper that never saw this code.

Two consequences of the choice worth stating before the result:

- **`--size 512` is forced.** `rs_subaperture_split()` needs `n_az >= 2*n_looks`,
  so 182 looks need 364 azimuth lines. That costs 4x the backprojection of the
  256-cell runs, and it produces **961 tracking windows** rather than 225.
  `src/main.c` carries a warning about exactly this: on this collect the same
  chain "refused honestly at 225 windows and reported 0.183 Hz, prominence 29.9
  at 961, off two crossings". The eligible-window count matters more than the
  headline number here.
- **The sub-looks are properly sampled for once.** `t_sap = 0.3 s` gives roughly
  5.8 m azimuth resolution, which a 2.0 m cell oversamples. Every previous Giza
  run carried the aliasing warning; this one earns it least.

### Where this still departs from the paper

Recorded so the run is not oversold as a replication:

- **Patch shape.** They use 5 range x 51-131 azimuth pixels, scaled to target
  velocity. `--win` is square, so this is 32 x 32. They report tracking quality
  is "predominantly" set by the azimuth dimension and that too small a patch
  underestimates the movement.
- **Reference frame.** They correlate against the **middle** sub-aperture; the
  reference here is frame 0 (`--reference first`). ResonarSat has no middle-frame
  mode.
- **No sinc interpolation.** Their `B = 20` upsampling of the velocity series has
  no equivalent here. It affects amplitude, not frequency.
- **The target.** Theirs is an isolated corner reflector on a shaker in an open
  field, chosen to remove clutter. Ours is a pyramid.
- **Velocities.** Their nine tests span 1.42 to 95.52 mm/s peak radial velocity,
  and the smallest confirmed RMS radial displacement is 0.10 mm. Nothing
  establishes that Khufu moves anywhere near that.

## Command

```sh
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 512 --cell 2.0 \
    --subap uniform --overlap 0.4 --estimator correlation \
    --n 182 --win 32 --coherence 0.4 --upsample 10 --null-trials 8 \
    --shifts spot_n182_ov40.csv --out spot_n182_ov40
```

`--upsample 10` is the paper's `A = 10`. `--null-trials 8` is a first look; the
shuffle is a *valid* null here, unlike in the phase runs, because a correlation
observable reads each look independently and reordering does not change the
per-look noise. That is the exemption written into
`rs_warn_shuffle_null_on_phase()` today, and this is the case it exempts.

## THE PREDICTION, recorded before the result

**No detection, and the most likely form is an honest refusal rather than a
spurious peak.**

- **P9** `RS_ERR_RANGE` -- no window resolves motion above the quantisation
  floor. `--upsample 10` puts that floor at 3 sigma of 1/10 px, i.e. 0.245 px,
  four times coarser than the 1/40 px used by the patent-chain runs, which found
  a largest excursion anywhere of 0.05 px. `--reference first` accumulates full
  displacement rather than differencing across a lag, so excursions should be
  larger than those runs saw -- but not, on any available evidence, by a factor
  of five.
- **P10** If it does report a peak, it comes from very few of the 961 windows,
  and the fmin ladder relocates it, exactly as in
  `2026-07-30-uniform-phase-khufu`. Watch the eligible count, not the prominence.
- **P11** The 8-shuffle floor behaves itself this time: the shuffled prominences
  should sit at a similar level to the measurement rather than 5x below it, since
  correlation does not have the red-noise asymmetry that made the phase runs beat
  their own shuffles.

**What would falsify it:** a peak in a contiguous patch of at least four windows,
clear of the lowest bins, surviving `--fmin`, above the shuffled floor, with
per-window excursions comfortably above 0.245 px. That would be the first
positive micro-motion result this project has produced, and it would need
repeating at `--upsample 40` and against `--null-static` before being believed.

**A note on what a null here does and does not mean.** Unlike the patent chain's
null, a null at this operating point cannot be dismissed as a structural artefact
of the method -- this configuration demonstrably measures 1-4 Hz vibration on a
corner reflector with millimetre displacements. A null here is evidence about
**Khufu**, or about the difference between a corner reflector and a limestone
massif, and not about the processing. That makes it the most informative null
this project could produce, which is why it was worth stopping a three-hour job
to run it.

The three `--shifts` series are stored gzipped, as in the sibling run: 961
windows by 159-182 looks is 8-13 MB of text each. `gunzip -k` before pointing
`shift_stats.py` or `freq_map.py` at them.

## Result of the first attempt: nothing tracks at all

```
sub-apertures: 182 looks, dt 0.1799 s
  observable band  f_max 2.78 Hz   AT sub-look resolution 5.58 m
sub-pixel refinement: 1/10 px
tracked 961 windows (31 x 31); 0 pass the 0.40 coherence mask
  no window is coherent enough to carry a measurement.
spectra: 92 bins, 0.0305 Hz resolution
mmotion: value out of range (spectrum: no window resolved motion above the
  0.2449 px quantisation floor (3 sigma of 1/10 px))
```

**P9 held exactly**, including the predicted 0.2449 px floor. But it held for a
reason the prediction did not anticipate: not that the excursions were too small
to clear the floor, but that **no window correlated at all**.

From `spot_n182_ov40.csv`, over all 961 windows:

| | value |
|---|---|
| tracking quality, median | 0.0505 |
| tracking quality, maximum | **0.1080** |
| windows at or above 0.4 | **0 of 961** |
| windows at or above 0.2 | 0 of 961 |
| windows at or above 0.1 | 3 of 961 |

**A caution about the shift series in that CSV: it is identically zero in every
window, and that is an artefact, not a measurement.** `src/core/microm.c:588`
zeroes `disp_az`, `disp_rg`, `vel_los` and `disp_los` for any window failing the
coherence mask -- deliberately, so an incoherent window contributes a flat
spectrum instead of a confident-looking peak. Since all 961 failed, all 961 were
zeroed. Nothing about the ground follows from those zeros.

**An unmasked re-run of this configuration was started and abandoned** -- same
parameters with `--coherence 0 --upsample 40 --null-trials 8` -- killed at 36
minutes, during the eighth of nine tracking passes it did not need. No headline
result was produced and none is claimed.

**Its `--shifts` CSV is complete, though**, because `--shifts` is written after
tracking and before the null trials, and 174,902 rows is exactly 961 windows x
182 looks. So the diagnostic the run was for survived the kill:

```
alpha = 0.91%, unmasked:  quality median 0.0506, max 0.1082
                          azimuth excursion median 31.975 px, max 33.000 px
                          (the tracking window is 32 px)
```

The zeros in `spot_n182_ov40.csv` were the mask. Underneath it the tracker was
producing peak-to-peak excursions of very nearly the full window width -- the
same window-scale wander the corrected run below shows, and for the same reason.
Both aperture fractions fail identically once the mask is removed; the mask was
only hiding it at 0.91%.

## The diagnosis: the aperture fraction was wrong, by a factor of five

Reading Lotti et al., *Feasibility of micro-motion from SAR imagery for
vibration-based SHM* (SHMII-13, DOI 10.3217/978-3-99161-057-1-142) -- the same
Umbra campaign as Vattulainen, same shaker and LVDT -- names this failure
directly:

> When the aperture fraction becomes too small, the target may no longer appear
> as a distinct feature, thus degrading the tracking performance. This imposes a
> practical limit on how finely the aperture can be segmented, although this can
> be mitigated by increasing the overlap.

The aperture fraction `alpha = t_sap/t_a` is the parameter this run got wrong:

| | `t_a` | `t_sap` | `alpha` |
|---|---|---|---|
| Lotti Test 1 | 6.04 s | 0.46 s | 7.6% |
| Lotti Test 2 | 5.95 s | 0.27 s | 4.5% |
| Lotti Test 3 | 14.70 s | 0.66 s | 4.5% |
| **this run** | 32.869 s | 0.300 s | **0.91%** |

The parameters here were derived by holding the *observation ratio* `eta` at 0.6
and letting `alpha` fall where it would. It fell five to eight times below the
validated range, each sub-look resolving 5.58 m, and the scene stopped
correlating. That is the mechanism Lotti et al. describe, observed.

**This also refines the criticism of `README.md` recorded in the sibling run.**
The README's numbers were wrong -- the validated work uses a few dozen looks at
23-49% overlap, not thousands at 90-99%. But its *mechanism* was right, and on a
32.869 s dwell high overlap is not optional. Holding `alpha` at 5% means
`t_sap = 1.64 s`, and the only way to still sample near 5 Hz is `Omega ~ 0.88`.
Their short dwells bought adequate `f_s` at moderate overlap; this collect
cannot. The two papers' own mitigation sentence says exactly this.

## The corrected configuration

```sh
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 512 --cell 1.0 \
    --subap uniform --overlap 0.88 --estimator correlation \
    --n 159 --win 32 --coherence 0 --upsample 40 --null-trials 8 \
    --shifts alpha5_n159_ov88.csv --out alpha5_n159_ov88
```

```
alpha 5.0%   t_sap 1.643 s   dt 0.1972 s
f_s 5.07 Hz  f_max 2.54 Hz   sub-look azimuth resolution ~1.06 m
N = 159      n_az 512 >= 2N = 318
```

`--null-trials` is deliberately **absent**, after being specified and abandoned
twice. At `--upsample 40` with `--coherence 0`, every one of the 961 windows is
tracked and sub-pixel refined, and the refinement cost is quadratic in the
refinement window -- so 1/40 px costs roughly sixteen times 1/10 px per call.
`--null-trials 8` then repeats that whole pass eight more times. The unmasked
diagnostic below was killed at 36 minutes without producing output for exactly
this reason, and the first attempt at this run repeated the error. A shuffled
floor is only worth paying for once there is a peak to test; if this run produces
one, the floor can be measured then.

Three deliberate changes from the first attempt:

- **`--cell 1.0`** rather than 2.0, because the sub-looks now resolve 1.06 m. This
  is the first Giza run in this project whose tracking images are properly
  sampled by the grid rather than aliased onto it.
- **`--coherence 0`**, so the series are not zeroed and the quality distribution
  can be inspected. Any threshold can be applied offline from the CSV, which the
  first attempt's output made impossible.
- **`--upsample 40`** rather than the paper's `A = 10`, because their targets
  moved at 1.42-95.52 mm/s and a 1/10 px floor was ample for that. Nothing
  suggests Khufu is in that range.

### Prediction for the corrected run

- **P12** Tracking quality rises substantially -- median well above 0.1, and a
  usable number of windows above 0.4. This is the parameter the diagnosis says
  was wrong, so if quality does not move, the diagnosis is wrong.
- **P13** Still no detection. Quality is necessary, not sufficient; the validated
  targets are corner reflectors moving at millimetres per second, and a limestone
  massif is neither. Expect either `RS_ERR_RANGE` again -- now for the honest
  reason, excursions below the floor rather than nothing to track -- or a peak
  that the fmin ladder relocates.
- **P14** If a peak does appear, the 8-shuffle floor is a valid test of it, unlike
  in the phase runs.

## Result of the corrected run: coherence doubles, and Khufu still does not track

```
sub-apertures: 159 looks, dt 0.1976 s
  observable band  f_max 2.53 Hz   AT sub-look resolution 1.02 m
sub-pixel refinement: 1/40 px
tracked 961 windows (31 x 31); 961 pass the 0.00 coherence mask
spectra: 80 bins, 0.0318 Hz resolution
strongest peak in window 455: 0.064 Hz, prominence 25.4, quality 0.128,
                              peak-to-peak velocity 278.3 mm/s
  526 of 961 windows were eligible for selection
  sub-aperture response 0.9820 (-0.2 dB) at an observation ratio of 0.10
```

### P12 half held: the diagnosis was right about the direction and wrong about the cure

| | first attempt, `alpha` 0.91% | corrected, `alpha` 5.0% |
|---|---|---|
| tracking quality, median | 0.0505 | **0.1226** |
| tracking quality, maximum | 0.1080 | **0.2439** |
| windows at or above 0.4 | 0 of 961 | **0 of 961** |
| windows at or above 0.3 | 0 of 961 | 0 of 961 |
| windows at or above 0.2 | 0 of 961 | 24 of 961 |

Raising the aperture fraction to the validated range roughly **doubled** the
coherence, which is the prediction the diagnosis rested on and it held. It did
not get the scene anywhere near trackable. Nothing reaches the 0.4 default mask
and only 24 windows of 961 clear even 0.2.

### P13 held, and the shift series says why

```
azimuth excursion, median   31.98 px      (the tracking window is 32 px)
azimuth excursion, maximum  33.00 px
windows above the 0.0612 px quantisation floor   961 of 961
windows with exactly zero excursion                0 of 961
```

**The correlation peak is wandering across the entire patch.** A median
peak-to-peak excursion of 31.98 pixels in a 32-pixel window is not a measurement
of motion; it is a correlator that finds no match and puts the peak wherever
noise is highest. Everything downstream of that is arithmetic on noise, including
the reported "0.064 Hz, prominence 25.4, 278.3 mm/s" -- which is, once again, the
second bin of the spectrum.

This is a different failure from the first attempt, and the difference is the
informative part. There, the mask zeroed every window and the series were
identically zero by construction. Here nothing is masked, every window carries a
series, every series clears the quantisation floor, and they are all garbage.
The tracker has stopped declining and started producing.

### WITHDRAWN: this run does not establish anything about Khufu

**The paragraph that stood here claimed the null was attributable to the target
rather than the configuration. A positive control run afterwards shows that it
is not, and the claim is withdrawn rather than edited, because it was wrong in
the direction that flatters the work.**

The control put a known vibrating target through this exact configuration --
`--subap uniform --overlap 0.88 --estimator correlation --n 159 --win 32
--upsample 40`, aperture fraction 5.01% against this run's 5.00% -- on a clean
synthetic collect with no clutter and nothing else in the scene. Target at
0.5 Hz, observation ratio 0.50, inside the 0.39-0.69 range the published
validation works at, injected displacement 4.19 px peak-to-peak against a
0.061 px quantisation floor.

The four windows containing the target returned:

```
  win 180  quality 0.087  p2p 30.90 px  dominant 1.203 Hz
  win 181  quality 0.116  p2p 10.48 px  dominant 1.569 Hz
  win 199  quality 0.078  p2p 32.00 px  dominant 1.046 Hz
  win 200  quality 0.082  p2p 31.93 px  dominant 0.785 Hz
```

None recovers 0.5 Hz. Excursions sit at the full 32-pixel window width and
quality at 0.08-0.12 -- the same signature, and the same numbers, this run
reported over Khufu. **A configuration that cannot find a bright isolated target
vibrating at a known frequency on an empty synthetic scene says nothing about a
pyramid.**

What is NOT established is why. The candidates are the uniform spectral-split
route, the 0.88 overlap, the window size against the 1.02 m sub-look resolution,
and an outright defect in this implementation of the chain. Distinguishing them
is the next work, and until it is done every number above is uninterpretable
rather than negative.

### A second defect in this configuration, found by the same control

The sub-aperture response nulls at integer observation ratios, `eta = t_sap * f`.
This run's `t_sap` is 1.643 s, so within its own 2.53 Hz observable band there
are nulls at **0.609, 1.217, 1.826 and 2.435 Hz** -- and `eta < 1` only below
0.609 Hz. The configuration was therefore blind to displacement above about
0.6 Hz before the scene was even considered, which is not something this run
recorded or the band figure printed alongside it suggests.

That is the same class of defect the patent chain was criticised for in
`2026-07-29-patent-exact-true-khufu`, where every resolvable bin sat at an
integer observation ratio. It is milder here -- four nulls rather than all of
them -- but it was found by accident, in a control whose first attempt injected
1.0 Hz against a 1.002 s sub-aperture and landed exactly on a null.

Two details that were offered as making it hard to blame the processing, and
that the control shows do not carry that weight:

- **The displacement-averaging problem is gone.** Sub-aperture response is 0.9820,
  which is -0.2 dB, at an observation ratio of 0.10. The mechanism that made the
  patent chain structurally blind -- integer observation ratios, exact nulls --
  does not apply anywhere in this configuration's band.
- **The sub-looks are properly sampled for the first time.** 1.02 m azimuth
  resolution on a 1.0 m grid, so the tracking images are neither aliased nor
  smeared, unlike every previous Giza run.

That reasoning was wrong, and the control is why. An engineered corner reflector
on an empty field is exactly what the synthetic control provides -- a single
bright point, no clutter, no speckle to decorrelate -- and the configuration
failed on it too. So "the scene" cannot be the remaining explanation. Whatever
defeats this configuration defeats it on the easiest target there is.

## Was the correlation failure our own interpolator? No.

Range interpolation in `rs_focus_backproject()` was linear between natively
spaced bins. NGA's reference backprojector zero-pads eightfold before doing the
same thing, on the stated grounds that linear interpolation otherwise "causes
cross-range artifacts" (MATLAB_SAR, `Processing/IFP/BP/bpBasic.m`).

That raised a specific hypothesis about this run rather than a general worry.
The interpolation error depends on the fractional bin position; that position
varies with slant range; different sub-apertures use different pulses and so
sit at different fractional positions. So the error differs BETWEEN LOOKS, which
is a mechanism for look-to-look decorrelation -- and look-to-look decorrelation
is exactly what failed here.

Re-run with an 8-tap Hann-windowed sinc, identical in every other parameter:

| | linear | 8-tap sinc |
|---|---|---|
| tracking quality, median | 0.1226 | **0.1211** |
| tracking quality, maximum | 0.2439 | **0.2321** |
| windows at or above 0.2 | 24 | 22 |
| windows at or above 0.4 | 0 | 0 |
| azimuth excursion, median | 31.975 px | **31.975 px** |

**The hypothesis is dead.** Coherence does not improve; it is marginally lower,
within noise. The median excursion is identical to three decimals -- the
correlation peak still wanders the full 32-pixel window.

This strengthens the null rather than qualifying it. The failure to correlate at
Khufu is not an artefact of this project's range interpolation, and one more
candidate explanation is retired. What the wider kernel does buy, measured
separately on a focused patch of the same collect, is about 1.4 dB of amplitude
and no change in image contrast -- so it is real, modest, and irrelevant to any
conclusion here. It stays off by default.

## The CCD locator on the same stack

Run after `rs_ccd_locate()` was implemented and wired to `mmotion` as
`--ccd-out`, on the identical configuration, so the locator and the tracker see
the same 159 sub-apertures:

```sh
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 512 --cell 1.0 \
    --subap uniform --overlap 0.88 --estimator correlation \
    --n 159 --win 32 --coherence 0 --upsample 40 \
    --ccd-out khufu
```

```
CCD locator: 5x5 window over 157 sub-aperture triples
  statistic  min 7.352  mean 26.608  max 3828.335   (1.000 is the no-change value)
```

### The statistic is saturated across the whole scene

| | synthetic fixture | **Khufu** |
|---|---|---|
| minimum over the map | 1.000 | **7.352** |
| mean | 1.359 | 26.608 |
| median | -- | 14.43 |
| maximum | 35.47 (the vibrating target) | 3828.3 |

**No pixel anywhere on the Khufu map is near the no-change value.** The lowest
point of the entire scene is 7.35, where the synthetic fixture's background sits
at 1.0 and only the deliberately vibrating target rises above it.

The reading is that real distributed clutter violates the detector's null
hypothesis everywhere. `Sigma_X = gamma * Sigma_Y` says the covariance structure
of consecutive sub-aperture pairs differs only by a scale factor; sub-look
speckle decorrelates between looks even at 88% bandwidth overlap, and
decorrelation is a change in structure, not in scale. So every window "changes",
and the statistic measures how much clutter decorrelated rather than what moved.

This is a property of the scene, not of the implementation -- the same code puts
a static point target at 1.049 on the fixture -- and it is why the source paper's
missing detection threshold matters more on rubble than on a shaker in an open
field. It is also worth noting that the paper's statistic maps are normalised to
their own maximum, which would render a saturated background invisible in
exactly this way.

### The extremes are stripes, not targets

The top 0.1% of pixels are not scattered and not clustered on a structure. They
lie in a **five-column band, 464-468, spanning rows 87-463** -- nearly the full
height of the image -- with a second weaker band at columns 46-47:

```
column median, typical            14.41
columns 464/465/466/467/468      498 / 601 / 652 / 632 / 576
columns 463 and 469               40 and 82        (abrupt on both sides)
```

Image column is the grid's Y, so a column-aligned band is a feature at constant
slant range spanning the whole scene. Nothing about a vibrating point target
produces that; a processing artefact or a linear ground feature at constant
range does. **The map's maximum is not a detection**, and reporting 3828 without
saying where it sits would have been the single most misleading number this run
could produce.

### Khufu itself

| region | median | mean |
|---|---|---|
| Khufu, 120 x 120 m about the grid origin | 17.53 | 18.61 |
| wider plateau, 360 x 360 m | 15.09 | 15.98 |
| off-target corners | 15.73 | 60.41 |

The pyramid runs about **16% above** the surrounding plateau, in a map whose
interquartile range is 11.9 to 18.0 and whose tail reaches into the thousands.
That is inside the ordinary variation of the scene and is not a detection of
anything. It is also the direction a slightly different backscatter statistic
would push it, with no motion involved at all.

### What this establishes, and what it needs next

The locator runs on real data, on the same stack as the tracker, and produces a
map where the tracker could not produce a measurement -- which was the point of
implementing it, since it never tracks anything and so is not blocked by the
correlation failure documented above.

**It has not detected micro-motion at Khufu, and it cannot say it has not.** A
map with no floor is not evidence, which the command prints for itself every
time it runs. The floor is `--null-static` through this identical chain: the same
geometry, the same 159 sub-apertures, the same 5x5 window, over a simulated scene
where nothing moves. Whatever value that produces is the number 14.43 has to
clear, and on the evidence here it will not be anywhere near 1.0.

That run is the obvious next one and was not done today.

### What else was not done

`--null-trials` was dropped, so there is no shuffled floor here. It would have
been valid for this observable -- correlation, not phase -- but a floor is a test
of a peak, and with excursions at the window scale there is no peak to test. It
is the obvious addition if any future configuration produces windows that
actually correlate.
