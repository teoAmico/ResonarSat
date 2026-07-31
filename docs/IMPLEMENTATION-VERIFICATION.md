# Independent implementation verification

First audited: 2026-07-29, against the tree at the initial commit
Revised: 2026-07-31, against master. The original findings stand; several are
sharpened below with results obtained since, and three defects in this
implementation — not in the patent — were found and corrected in the interval.

## Question

Does ResonarSat correctly implement WO 2024/008365 A1, and do independent
micro-motion papers validate the physical and signal-processing assumptions on
which that implementation depends?

## Short answer

The project implements a coherent, conventional interpretation of the patent's
computational chain closely enough to test it, including the two Doppler
filters, rigid master/slave sweep, sub-pixel tracking, construction of the
complex shift vector, and steering-matrix pseudoinverse. It does **not**
implement every equation literally as printed in the rendered patent. Most
importantly, Eq. 22 contains an extra `2πt` in its steering exponent. The
default code does not apply it; a fixed caller-supplied interpretation is
available only as experimental `--eq22-literal-t`. The printed equation and
its matrix dimensions are internally inconsistent, so there is no unique
literal implementation. Where
the patent is silent or ambiguous, the program necessarily makes assumptions.
The most important are the value of `B_shift` and the interpretation of Eq. 22.

Independent literature supports the first half of the chain: appropriately
sampled SAR data can recover surface-target vibration frequency, and sometimes
time history and amplitude, using sub-aperture pixel tracking or phase. It does
not validate the second half: replacing the electromagnetic wavelength with an
assumed acoustic wavelength and interpreting an along-orbit sub-aperture
sequence as an elevation/depth aperture. No independent paper located for this
audit demonstrates subsurface depth recovery from a single SAR acquisition.

The present tests therefore establish **software conformance and limited
synthetic micro-motion recovery**, not physical validation of the patent's
depth product.

## Sources used

### Primary specification

The publications themselves are not redistributed here; fetch them from the
publisher.

- WO 2024/008365 A1,
  https://patents.google.com/patent/WO2024008365A1/en
- Biondi and Malanga, *Synthetic Aperture Radar Doppler Tomography Reveals
  Details of Undiscovered High-Resolution Internal Structure of the Great
  Pyramid of Giza* (2022), Remote Sensing 14(20):5231,
  https://doi.org/10.3390/rs14205231
- Biondi, *Scanning Inside Volcanoes with Synthetic Aperture Radar Echography
  Tomographic Doppler Imaging* (2022), Remote Sensing 14(15):3828,
  https://doi.org/10.3390/rs14153828

### Micro-motion cross-checks

- Biondi et al., *Micro-Motion Estimation of Maritime Targets Using Pixel
  Tracking in COSMO-SkyMed SAR Data—An Operative Assessment* (2019),
  https://doi.org/10.3390/rs11141637
- Biondi et al., *Monitoring of Critical Infrastructures by Micromotion
  Estimation: The Mosul Dam Destabilization* (2020),
  IEEE JSTARS 13:6337-6351, https://arxiv.org/abs/2007.05326
- Clemente et al., *On Micro-Motion Extraction from High Resolution X-band SAR
  Products* (2024), https://doi.org/10.1109/IGARSS53475.2024.10642825
- Rollo et al., *Micro-motion Extraction from Spotlight SAR Using a Modified
  Backprojection Approach* (2024),
  https://strathprints.strath.ac.uk/90732/
- Vattulainen et al., *Assessment of Spaceborne SAR Micro-Motion Measurement
  for Vibration-Based SHM* (2026),
  https://doi.org/10.1109/ACCESS.2026.3652346

The last paper supplies the strongest independent metrological evidence used
here: synchronous ground truth, 1–4 Hz motions, and radial RMS displacements
from 10.43 mm down to 0.10 mm.

## Patent-to-code conformance

| Patent operation | Implementation | Assessment |
|---|---|---|
| Block 1: complex SLC or raw SAR input | CPHD and SICD readers | Implemented. Capella CI4 and Umbra CF8 CPHD are supported. |
| Block 2: DFT2 of SLC | Azimuth FFT for every range column | Equivalent for the patent's azimuth-only Doppler filtering because the full range band is retained. It is not literally a 2-D transform, but the omitted range FFT/filter/IFFT pair would cancel. |
| Blocks 3–4: master and slave Doppler band-pass filters | `RS_SUBAP_PAPER`, optional paired slave stack | Implemented. Processed width is `B_CD-B_DL`; default `B_DL/B_CD=0.5`. |
| Rigid sweep in steps `(B_CD-B_DL)/N_D` | Paper-mode band centres | Implemented from the printed formula. |
| Free separation `B_shift` | `b_shift_hz`; zero derives it from the sweep step | Exposed, but **not specified by the source**. The default is an implementation assumption, not a reproduction of a published numeric value. |
| Blocks 5–6: inverse transforms | Inverse azimuth FFT | Equivalent under the same separability condition as Block 2. |
| Block 7: sub-pixel pixel tracking | Normalised cross-correlation with local DFT peak refinement | Implemented. The patent gives no reproducible coregistrator details, so the exact estimator cannot be matched. |
| Master/slave tracking at each sweep position | `RS_MICROM_REF_PAIR` | Implemented and forced by `--patent-exact`. |
| Eq. 20–21 complex vector from two shift components | Azimuth metres + `j` × range metres | Reasonable implementation of an underspecified equation. |
| Eq. 22 steering matrix | Default `exp(j Kz_i z_j)`; experimental fixed-`t` `exp(j 2π Kz_i t z_j)` | The source does not define how the acquisition-time variable becomes the unindexed `t` in every matrix element. No unique literal discretisation exists. |
| `Kz=4πB_perp/(λ r sinθ)` | `rs_tomo_wavenumber()` | Algebraically implemented. The physical meaning of the supplied baseline remains unverified. |
| Eq. 23 dimensions | `Y(k×1)=A(k×F)h(F×1)` | The code uses the only dimensionally consistent orientation. The patent declares `A∈C^(k×F)`, draws an `F×k` matrix, and declares `h∈C^(1×F)`. |
| Eq. 24 pseudoinverse | Shape-aware Moore–Penrose solve | Implemented for the corrected matrix orientation; `--patent-exact` disables regularisation. |
| Acoustic wavelength and quoted resolution | Paper and patent conventions both available | Implemented, including their factor-of-two disagreement. |
| Geocoding | CSV depth-point output | Partially implemented; not a standard 3-D GIS raster/volume. |

### Exact-mode additions and departures

The historically named `--patent-exact` flag forces the patent-chain choices:

- paper spectral decomposition;
- rectangular sub-aperture filters;
- rigid master/slave pair tracking;
- complex coregistrator shifts as `Y`;
- the patent `v/f` wavelength convention;
- an unregularised pseudoinverse;
- no mean removal, taper, or coherence rejection.

Normal/default operation adds tapering, mean removal, coherence filtering, and
regularisation. Those are sensible numerical controls, but a product made with
them is not the patent-chain interpretation. Even with the flag, the result
must not be labelled “Eqs. 21–24 as written,” because the code resolves the
source's Eq. 22–23 inconsistencies rather than reproducing them.

## Findings from independent micro-motion literature

### Supported

1. **A complex SAR acquisition contains target-motion information.**
   Motion along the radar line of sight changes Doppler/phase and commonly
   appears as azimuth displacement or paired echoes.

2. **Sub-aperture processing can recover a temporal sequence.**
   The 2019, 2020, 2024, and 2026 work all use aperture segmentation,
   pulse-resolved backprojection, or both.

3. **There is an unavoidable time/resolution trade-off.**
   Shorter sub-apertures improve temporal reach but degrade azimuth resolution.
   Heavy overlap can increase the sample rate without shortening each window.

4. **SPOT/correlation can recover vibration in a suitable regime.**
   The 2026 assessment validates SPOT against synchronous instrumentation.
   Patch size, upsampling, target velocity, clutter, and resolution materially
   affect the result.

5. **Phase is a different observable, not merely another coregistrator.**
   Modified backprojection can expose phase pulse by pulse and is more
   sensitive to small motion, but clutter and phase unwrapping are practical
   limitations.

6. **Frequency is generally more robust than amplitude/time history.**
   This supports the project's policy of treating amplitude conservatively.

### Not supported by those papers

1. **The patent's rigid `B_shift` pair is not the independently validated SPOT
   configuration.** Earlier SPOT papers describe adjacent/sliding temporal
   sub-apertures. The 2026 assessment uses an overlapping sequence reconstructed
   from an SLC. The Capella modified-backprojection experiment reads target
   phase from a pulse-contribution cube. None independently validates the
   patent's rigid master/slave pair plus depth inversion.

2. **The Capella 2 Hz experiment does not validate this patent chain.**
   It used CPHD, modified backprojection, and direct phase from a controlled
   corner reflector. It did not use the patent's master/slave pixel-offset
   vector or produce a subsurface tomogram.

3. **No independent depth mapping is supplied.**
   The micro-motion papers measure motion of visible surface targets. They do
   not establish that surface vibration phase can be focused into subsurface
   reflectivity or depth using Eq. 22–24.

## Critical specification problems

### 1. `B_shift` is not reproducibly specified

The patent says that `B_shift` selects a desired mechanical frequency and that
larger separation observes lower frequency, but gives neither a numeric value
nor a selection equation. The implementation derives it from the sweep step
when the user does not provide it. This is transparent and testable, but cannot
be called uniquely correct.

There is also a geometric conflict: with `N_D` master windows of width
`B_CD-B_DL` spanning the available band, only one sweep-step of spectral
headroom remains for the final slave. A larger freely selectable `B_shift`
places the last slave outside the measured Doppler support.

The patent's one stated property implies the rule it never gives. The pair
response `|2 sin(pi f dt)|` peaks at `f*dt = 1/2`, so the lag that maximises
sensitivity to frequency `f` is `B_shift_opt = B_CD/(2*f*t_dwell)` -- inversely
proportional to `f`, exactly as the patent describes. Imposing the geometric cap
above gives `f >= N_D/(2*t_sap)`: **the lowest frequency `B_shift` can be tuned
to is `N_D` times the highest the layout can observe**, independent of dwell and
collect. Reading `B_shift` as independent of the sweep, which the patent's
wording allows, lets the sweep span shrink to buy separation; the two then
compete for the same headroom and the requirement becomes a factor `3*(1-L)/L`
above the band edge, which is 3 at the stated `B_DL = B_CD/2`. Neither reading
reaches an observable frequency. See `docs/CORE-QUESTIONS.md` question 6.

### 2. The default patent pair is a temporal difference

When `B_shift` equals the sweep step, `slave(i)` and `master(i+1)` are the same
spectral window. Pair tracking therefore measures a fixed-lag first difference
without accumulation. That observable has response

`|2 sin(π f Δt_pair)|`

and is weak at low frequency. The code models this fact, but the patent does
not describe the required response correction.

Quantified at the band edge, with `B_shift` pinned to the printed step: the
differential costs 20 dB at `N_D = 32`, 32 dB at 128 and 44 dB at 512 --
worsening with the very parameter one would raise for finer time sampling. For
contrast, Multiple Aperture InSAR (Bechor & Zebker 2006) performs the same
two-sub-aperture operation at a normalised squint of 0.5, half the total Doppler
band, which is 128 times the paper's separation at `N_D = 128`.

### 3. The printed paper layout creates a measurement null

With `B_DL=B_CD/2` and the printed step,

`N_D Δt = t_subaperture`

and the spectrum spacing is `1/t_subaperture`. A correlation tracker measures
position averaged across the sub-aperture; a rectangular sub-aperture has sinc
zeros at integer multiples of that spacing. Thus the exact patent layout places
its resolvable bins at the averaging response's nulls.

This is now exact rather than qualitative. Since the bin spacing is `1/t_sap`
while the averaging response reaches only `1/(2*t_sap)`, **the first bin sits at
exactly twice the band edge, for every `N_D` and every `B_DL`** -- no parameter
choice moves it, and it falls precisely on the first averaging null. Measured on
a synthetic collect at 20 s dwell: 0.100 Hz reported against a 0.050 Hz band,
sub-aperture response -240 dB, observation ratio exactly 1.00.

The current Giza run observes this predicted failure: 225/225 windows remain
below the tracking quantisation floor.

### 4. `B_perp` is not provided by a single along-track aperture

Eq. 22 requires an orthogonal/elevation baseline. The implementation, following
the project's reading of the sources, substitutes progressive along-orbit
sub-aperture phase-centre separation. An along-track baseline is not generally
an elevation baseline. Independent SAR tomography normally obtains elevation
wavenumber diversity from genuinely different look geometries/baselines.

The standard decomposition makes this sharper than "an along-track baseline is
not generally an elevation baseline". Bähr (*Orbital Effects in Spaceborne SAR
Interferometry*, DGK Reihe C 719, KIT 2013, Sect. 3.1) gives the cross-track
form used for height as valid only "assuming `B_a ~ 0`", `B_a` being the
along-track component, and warns it otherwise biases the height ambiguity. A
sub-aperture separation is *entirely* `B_a`.

The consequence is not that `B_perp` vanishes. Under the zero-Doppler condition
the line of sight is perpendicular to the velocity, so the general scalar form
`B_perp = sqrt(|B|^2 - B_par^2)` returns `B_par = 0` and `B_perp = |B|` -- the
whole separation, finite and ordinary-looking, from a vector that produces no
parallax across the line of sight. `K_z` accepts the number and returns a normal
depth axis. On the Giza collect the along-orbit arc is 238.7 km against a
repeat-pass baseline spread of order 0.5 km, so `delta_T = lambda*R/(2*A)` reads
0.051 m where a genuine baseline gives 24.6 m.

Moreira et al. draw the same distinction from the other side: across-track
separation gives elevation sensitivity, along-track separation gives a temporal
lag and hence velocity sensitivity. They are different observables.

This is the main physical issue in the depth stage, not a numerical coding
detail.

### 5. The depth scale is imposed, not observed

The Giza paper assumes 12.5 kHz and the patent assumes 22 kHz, while the
spaceborne aperture records only seconds of data and the implemented
micro-motion products measure frequencies of order hertz. The code correctly
labels this an assumed investigation frequency. Consequently the reported
depth axis changes when the user changes `v` or `f`; it is not calibrated by a
measured 12.5/22 kHz signal.

### 6. The acoustic substitution lacks an independently validated forward model

The conventional SAR tomographic wavenumber contains the electromagnetic radar
wavelength. The patent substitutes an acoustic wavelength to obtain metre-scale
depth resolution, without independently demonstrating the measurement equation
that connects a buried acoustic source to the phase of surface radar
scatterers. Correctly evaluating that substituted equation does not validate
the substitution.

Taken with problem 4, `K_z` carries **two** substitutions at once: an
along-track separation standing in for an elevation-sensitive baseline, and a
mechanical wavelength standing in for the radar one. Either alone removes the
expression's derivation. `docs/CORE-QUESTIONS.md` question 4 now leads with the
wavelength substitution rather than with the `v/f` versus `v/(2f)` discrepancy,
which is the smaller half of it.

### 7. Eq. 22 across the source family

The patent PDF, printed page 14 (PDF page 16), shows steering terms of the form

`exp(j 2π k_z,i t z_j)`.

This was checked visually in the rendered PDFs, not inferred from OCR. The
notation has the following traceable lineage:

| Source | Equation | Exponent | Adjacent definition and consequence |
|---|---:|---|---|
| Vesuvius, *Remote Sensing* 14, 3828 (2022), PDF p. 7 | 10 | `j 2π k_z t z` | Defines `k_z=4πB_perp/(λr sinθ)` and later uses `δ=λR/(2A)`, with no `t`. |
| Giza arXiv preprint 2208.00811v1, PDF p. 11 | 22 | `j 2π k_z t z` | Repeats the same matrix, dimensions, `k_z`, and `t`-free resolution. |
| Giza, *Remote Sensing* 14, 5231 (2022), PDF p. 15 | 22 | `j 2π k_z t z` | Retains the expression through peer-reviewed publication. |
| WO 2024/008365 A1, PDF p. 16 | 22 | `j 2π k_z t z` | Repeats the same expression and matrix conflict nearly verbatim. |

The earlier geometry sections of both 2022 papers explicitly define `t` as
the **acquisition-time variable**, with `t=0` and `t=T` at the acquisition
limits. They do not redefine it before the steering equation. Equation 21 then
calls `Y` a vector of `k` samples of the time-domain oscillator. Nevertheless,
Eq. 22 writes the same unindexed `t` in every element rather than `t_i`, gives
no sample times or time step, and its resolution derivation contains no time
factor. The patent removes the earlier symbol list but does not add a
definition local to Eq. 22.

This establishes that the term is repeated source notation, not a Markdown,
OCR, or patent-typesetting error. Repetition does not make its numerical
meaning determinate.

#### Candidate interpretations

1. **`t` is physical acquisition time.** Then the exponent is not
   dimensionless: `(rad/m)·s·m` retains seconds. A sampled model would also
   require `t_i`, and changing the arbitrary time origin would change the
   inferred depth. This cannot be reconciled with the printed resolution.
2. **`t` is the sampling interval.** This still retains seconds in the
   exponent, is never stated, and would make resolution depend on the chosen
   sub-aperture sampling interval.
3. **`t` is a dimensionless sample index or normalized time.** This repairs
   units but is not defined and would duplicate sample dependence already
   represented by `k_z,i`/`B_perp,i`.
4. **`k_z` is a cycles-per-metre frequency.** Then a leading `2π` would be
   normal, but it conflicts with the immediately printed
   `k_z=4πB_perp/(λr sinθ)`, which is already an angular wavenumber. It still
   does not explain `t`.
5. **`2πt` is an editing residue.** Removing it produces the conventional
   TomoSAR steering law, matches the printed `k_z` units, and recovers the
   separately stated resolution scaling. This is the most internally coherent
   interpretation, but the sources never explicitly issue that correction.

The default code therefore evaluates `exp(j k_z,i z_j)`. Experimental
`--eq22-literal-t VALUE` instead evaluates the printed algebra using one fixed
user-supplied scalar: `exp(j 2π k_z,i VALUE z_j)`. That mode is useful for
sensitivity testing, but it is only the fixed-`t` interpretation. It does not
claim to implement a hypothetical `t_i` model, because the sources supply
neither that equation nor the mapping from sub-aperture samples to `t_i`.

The matrix dimensions also conflict:

- the prose declares `A∈C^(k×F)`;
- the displayed matrix has `F` depth rows and `k` baseline columns;
- Eq. 23 requires `A` to be `k×F` if `Y` is `k×1`;
- the patent declares `h∈C^(1×F)`, while Eq. 23 needs an `F×1` column.

ResonarSat chooses `A[i,j]=exp(j k_z,i z_j)`, `A∈C^(k×F)`, and
`h∈C^(F×1)`. That is a conventional and dimensionally valid repair, but not a
literal transcription. The author would need to define whether `t` is scalar
or indexed, give its units and sampling rule, reconcile it with `k_z` and the
resolution formula, and correct the matrix orientation before a uniquely “as
written” mode could exist.

### 8. Eqs. 1–20 are not all executable stages in the code

Equations 1–6 describe an idealised sinc/rect point-target model. ResonarSat
uses vendor data plus backprojection or spectral SLC filtering rather than
evaluating those formulas sample by sample. Equations 7–12 derive constant
range/azimuth velocity and acceleration perturbations; the tracker measures
image offsets without estimating all three printed epsilon parameters.
Equations 13–20 give a spring/nonlinear oscillator model; the program does not
fit its stiffness, damping, forcing, or oscillator parameters. It begins the
depth transform from the measured shift sequence corresponding to Eq. 21.

Accordingly, “implements the whole patent equations as written” is false.
“Implements the patent's Figure 0.5 processing chain under a documented
interpretation of Eqs. 21–24” is the defensible description.

## Defects found in this implementation, not in the patent

Distinguished from the specification problems above because these were ours and
are corrected. All three produced complete, plausible output while wrong, which
is the failure mode the whole project exists to guard against.

1. **The observable band was computed from the sub-aperture step, not from
   `t_sap`.** Each sub-aperture averages the motion over its own duration, so
   the series carries nothing above `1/(2*t_sap)` however finely overlap samples
   it; quoting the step's `1/(2*dt)` overstates the reach by `1/(1-overlap)`, a
   factor of 100 at 0.99 overlap. Every frequency configuration A of
   `runs/giza/2026-07-30-uniform-phase-khufu` reported lies above the band it
   could carry. This also restores `eta < 0.5` as a derived bound, since
   `eta = f*t_sap`.

2. **The coherence mask could be vacuous.** Overlapping sub-looks share spectral
   content by construction, so the pair-averaged estimator has a floor no scene
   can push it below. At 0.99 overlap that floor is 0.574 against the 0.4 mask
   the run actually used: it passed every window it was given.

3. **The sensitivity figure reported the sub-pixel interpolation limit**, not the
   excursion at which the tracker returns the target's frequency rather than its
   own artefact — optimistic by a factor of 57.

All three are now computed ahead of processing by `resonarsat validate`, whose
checks are pinned by `tests/test_validate.c` against configurations whose
outcome is already known from measurement. How the tracker's floor scales with
cell size, sub-look resolution and window size in pixels remains **unresolved**:
three candidate laws have been refuted and those factors are confounded in every
run made so far.

## Executed verification

Both Release and AddressSanitizer/UndefinedBehaviorSanitizer configurations
pass all 23 registered tests.

### How recovery is judged

A per-point match -- `|recovered - injected|` inside a tolerance at a single
frequency -- is too weak to establish anything, and was the source of five
conclusions drawn and withdrawn on 2026-07-31. A chain emitting a FIXED spurious
frequency passes it wherever that value happens to fall near the injection: one
configuration reported 1.569 Hz for every injection from 0.2 to 1.4 Hz and
scored "recovered" at several of them.

`rs_track_fit()` in `tests/rs_test.h` replaces it. The injected frequency is
swept and the reported value fitted against it, requiring slope near 1 -- a
fixed artefact gives 0, and no single point can tell the two apart -- and rms
under half a bin. `test_nullmotion` asserts both, plus a negative control that
must reject a constant series.

That criterion sweeps frequency, not speckle. A five-point sweep passing on one
clutter realisation was later found to pass a configuration that fails on
another, so a claim about a *configuration* needs the sweep repeated over
independent `--seed` realisations and the verdicts pooled.

The tables below therefore stand as measurements, while the "pass" verdicts in
the first of them rest on the older per-point test.

Measured synthetic correlation results at 20 mm amplitude:

| Injected | Recovered | Result |
|---:|---:|---|
| 0.30 Hz | 0.302 Hz | pass |
| 0.50 Hz | 0.503 Hz | pass |
| 0.70 Hz | 0.704 Hz | pass |
| 0.90 Hz | 1.811 Hz | false second harmonic; look count below required bound |
| 1.10 Hz | 2.213 Hz | false second harmonic; look count below required bound |
| 1.30 Hz | 1.308 Hz | recovered despite being below recommended look count |

At `test_nullmotion`'s operating point -- 128 looks, zero overlap, the settings
the phase-ambiguity condition demands -- the same sweep recovers 0.302, 0.504,
0.706 and 0.907 Hz for 0.3 to 0.9 Hz, four of four above the separately measured
static false-positive floor, with **slope 1.008 and rms 0.0052 Hz against a
0.0252 Hz bound**. The negative control, a constant 1.569 Hz, gives slope 0.000
and rms 0.9945 Hz. The false second harmonics in the table above occur at the
lower look counts, which the table marks as under-served.

The synthetic phase case appeared to recover 0.40 Hz as 0.407 Hz. **That is
withdrawn.** Swept rather than sampled, the phase estimator reports 0.407 Hz for
every injection from 0.2 to 0.7 Hz -- slope 0.000, rms 0.1762 Hz against a
0.0254 Hz bound -- and a target with no motion at all reports the same 0.407 Hz
at a HIGHER prominence, 12.5 against 8.0-10.8. It is bin 8 exactly. The single
point passed only because the fixed artefact happens to sit near 0.4 Hz, which
is the fifth instance of the failure `rs_track_fit()` was written to catch.
`test_tracking.c` now records it as a negative.

The patent master/slave pair did not recover its injected frequencies:

| Injected | Reported |
|---:|---:|
| 0.50 Hz | 0.100 Hz |
| 1.00 Hz | 3.300 Hz |

These are deliberate regression assertions in `test_tracking`; they prevent a
software change from silently converting a known negative result into a claim
of reproduction.

## Correctness verdict

| Claim | Verdict |
|---|---|
| The code can read and focus supported Capella/Umbra CPHD | Verified in code/tests and on the Giza Capella product |
| The code implements general sub-aperture micro-motion estimators | Verified synthetically |
| Correlation SPOT works over all advertised settings | Not verified; demonstrated parameter-dependent failures |
| Direct phase reproduces the Capella modified-BPA method | No; the observable is related, but the published pulse-contribution cube is not reproduced exactly |
| `--patent-exact` follows the patent's Figure 0.5 operator chain | Substantially yes, subject to the unspecified `B_shift` and unavailable original coregistrator |
| `--patent-exact` implements Eqs. 21–24 literally as rendered | Not by default; Eq. 22's `2πt` is available only through experimental `--eq22-literal-t`, while Eq. 22–23 dimensions remain repaired |
| The code implements Eqs. 1–20 as executable models | No; several are background derivations and the spring model is not fitted |
| The patent pair recovers known synthetic vibration | No, at the tested operating points |
| The phase estimator recovers known synthetic vibration | No, at twelve operating points; a static scene returns the same frequency |
| Either estimator recovers known synthetic vibration on distributed texture | No; slope swings sign across clutter seeds and static controls match the moving reports |
| Eq. 22–24 recover known synthetic physical depth from SAR motion | Not demonstrated by the current simulator; existing depth fixtures construct `Y` from the inversion model itself |
| Independent papers validate the surface micro-motion stage | Yes |
| Independent papers validate single-pass subsurface tomography | No evidence found |
| A ResonarSat depth peak is currently a measured subsurface depth | No |

### Eq. 22 comparison with conventional TomoSAR

Independent TomoSAR literature uses one of two equivalent conventions:
`exp(j*kz*z)` when `kz = 4*pi*B_perp/(lambda*r*sin(theta))` is an angular
spatial frequency in radians per metre, or `exp(j*2*pi*xi*z)` when `xi` is a
spatial frequency in cycles per metre. It does not multiply the first form by
another `2*pi`, and no independent source found introduces the patent's
additional `t`. Thus `--eq22-literal-t` is a reproduction/sensitivity switch,
not a validated correction to the conventional steering model.

For example, Cazcarra-Bes defines `kz` with `4*pi` in Eq. 2.6 and the steering
vector as `exp(j*kz*z)` in Eq. 2.9
([DLR/ETH dissertation](https://elib.dlr.de/133525/1/phd_thesis_victor_cazcarra_bes.pdf)).
Pardini et al. independently use the same `exp(j*kappa_z*z)` convention
([Remote Sensing 2021](https://elib.dlr.de/146465/1/Pardini-2021_remotesensing-13-02255-v2.pdf)).

## Required tests before a stronger claim

1. **Independent micro-motion reproduction**

   Acquire the controlled Capella/Umbra corner-reflector dataset used in a
   published experiment, freeze processing parameters before inspecting the
   answer, and reproduce frequency, amplitude, and timing against synchronous
   ground truth.

2. **Real-clutter injection test**

   Inject a known phase-history motion into a copy of a real static CPHD scene.
   Compare moving, static, time-shuffled, and wrong-ROI controls through exactly
   the same pipeline.

3. **Estimator-specific validation**

   Validate correlation, phase, and rigid-pair modes independently. A successful
   phase result cannot validate a pixel-offset or pair result.

4. **Forward-model depth simulation**

   Build a physical model that maps buried wave sources through surface
   displacement to radar phase history. Generate data without using the
   inversion steering matrix, then test blind depth recovery. Constructing `Y`
   directly from `A h` tests numerical inversion only and is circular as a
   physical validation.

5. **Depth nulls**

   Run static surface clutter, surface vibration with no buried source, changed
   assumed velocity/frequency, randomized sub-aperture order, and unrelated
   acquisition geometry. A depth feature must disappear or move only where the
   forward model predicts.

6. **Independent-baseline reference**

   Process a conventional multi-baseline SAR tomography dataset with known
   elevation targets. This validates the steering/inversion machinery separately
   from the acoustic hypothesis.

**Status at 2026-07-31: test 3 is executed and returns a negative; the other
five remain open.**

**Test 3, executed.** The phase estimator was swept over frequency and scanned
over twelve operating points -- 32, 64, 128 and 256 looks at 0.00, 0.50 and 0.75
overlap. **None tracks.** The best slope is +0.266 with 0.31 Hz rms, and in every
one of the twelve the static control lands on the same frequency the moving cases
report: the artefact's value moves with the configuration and not with the scene.
So the phase route fails the same way the correlation route does, and the
reasoning below -- which made test 3 the recommended next step -- has been
answered rather than deferred:

> Test 3 is the cheapest and is the recommended next. Every tracking failure
> recorded on 2026-07-31 is the *correlation* estimator; the phase route has a
> Cramér–Rao floor roughly 160x lower and has barely been exercised. Running
> `rs_track_fit()` over frequency and seed on the phase estimator would establish
> whether the difficulty is correlation-specific, and either outcome bounds the
> capability question more than further correlation sweeps would.

The difficulty is **not** correlation-specific.

**Nor is it a fixture artefact.** The remaining explanation was distributed
texture -- the fixtures above are isolated points on empty scenes, which
`microm.c` warns is the case correlation scores badly on. Swept on
`sim_cphd --clutter-vib`, where 400 Rayleigh scatterers vibrate coherently so the
patch moves as a whole, **neither estimator tracks either**. Five frequencies
from 0.3 to 1.1 Hz at a constant 10 px displacement, three clutter seeds, each
with a static control on the identical clutter:

| estimator | slope by seed | rms (Hz) |
|---|---|---|
| correlation | -2.217, +2.344, +1.487 | 0.79-1.29 |
| phase | -0.252, +0.000, -1.135 | 1.10-1.30 |

A tracking chain gives slope near 1 and rms under half a bin, which is 0.025 Hz
here. Correlation's slope swings sign across seeds; phase on one seed returns a
fixed 1.764 Hz for every injection. Every static control falls inside the spread
of the moving reports, and doubling the displacement changes nothing. The
operating point was admissible on this project's own arithmetic -- 10 px sits
mid-window between a 4 px textured floor and a 24.8 px ambiguity ceiling, well
inside a 3.23 Hz band at an observation ratio of 0.33.

**So the synthetic evidence for the tracker rests on exactly one configuration:**
correlation, an isolated point target, an empty scene, 128 looks and zero
overlap. That is the easiest case that exists and the furthest from a real
structure. Full numbers in item 2 of `FOLLOW-UPS.md`.

Test 1 is blocked externally rather than by effort: `docs/DATASETS.md` records
that no corner-reflector collect with synchronous ground truth exists in any
open archive, and every validated result in this literature used one.

Until those tests pass, surface vibration spectra may be reported with their
null controls and operating bounds. Patent-derived depth products should be
reported as experimental transforms under assumed `(v,f)` parameters, not as
validated measurements of internal structure.

