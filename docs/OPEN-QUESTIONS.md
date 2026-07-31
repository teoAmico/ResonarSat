# Open questions -- the detailed list


## Purpose

These questions collect the remaining ambiguities encountered while
implementing and independently checking the processing described in:

- *Scanning Inside Volcanoes with Synthetic Aperture Radar Echography
  Tomographic Doppler Imaging*;
- *Synthetic Aperture Radar Doppler Tomography Reveals Details of Undiscovered
  High-Resolution Internal Structure of the Great Pyramid of Giza*;
- WO 2024/008365 A1;
- the related SAR micro-motion publications.

They are requests for technical clarification, not assumptions about the
authors' implementation. The most useful response would be one complete
numerical example of Eqs. 21–24 containing the actual observation vector,
times, baselines, wavenumbers, depth grid, steering matrix, wavelength
convention, and resulting depth profile.

## 1. Equation 22 and the steering matrix

1. In `A(Kz,z)`, what exactly does `t` represent in
   `exp(j*2*pi*Kz*t*z)`?
2. What are the units of `t`?
3. Is `t` acquisition time, a sampling interval, normalized time, or a sample
   index?
4. Is `t` one scalar for the whole matrix, or should observation `i` use
   `t_i`?
5. If `t_i` is intended, how is it calculated from pulse time, sub-aperture
   centre time, Doppler centroid, or orbital position?
6. Is the time origin the aperture start, aperture centre, or zero-Doppler
   time?
7. If acquisition time is intended, how is inferred depth made invariant to
   the arbitrary time origin?
8. Is the leading `2*pi` intentional?
9. Is `Kz` measured in radians per metre or cycles per metre?
10. Given `Kz=4*pi*B_perp/(lambda*r_i*sin(theta))`, is the additional `2*pi`
    in the exponent redundant?
11. Is the intended conventional expression instead `exp(j*Kz_i*z_j)`?
12. Could `2*pi*t` be an editing residue repeated across the Vesuvius, Giza,
    and patent documents?
13. Was Eq. 22 evaluated literally for the published results?
14. What actual values of `t`, `Kz`, and `z` were used in one published
    experiment?
15. Can the source code or pseudocode that constructs `A` be provided?

## 2. Matrix dimensions and inversion

16. Is `A` intended to be `k x F` or `F x k`?
17. Do its rows represent sub-aperture observations or depth cells?
18. Is `Y` a `k x 1` column vector?
19. Is `h(z)` an `F x 1` column vector rather than the printed `1 x F` row?
20. Is the intended forward model `Y(k x 1) = A(k x F) h(F x 1)`?
21. Was the pseudoinverse computed using SVD, QR, normal equations, or another
    algorithm?
22. Was regularisation used?
23. If so, what type and parameter value were used?
24. Were small singular values truncated?
25. Was `Y` mean-centred, detrended, normalized, tapered, or
    coherence-weighted?
26. Was a window applied across sub-apertures or depth cells?
27. How was the condition number or inversion stability assessed?

## 3. Meaning and calculation of the baseline

28. What precisely is `B_perp,i` in a single SAR acquisition?
29. Is it a true perpendicular/elevation baseline or an along-track
    sub-aperture phase-centre separation?
30. If it comes from along-track motion, how is it projected perpendicular to
    the line of sight?
31. Can the complete vector expression used to calculate `B_perp,i` be
    provided?
32. Is Earth curvature included?
33. Are satellite and target positions evaluated in ECEF coordinates?
34. Is the baseline measured relative to the first look, central look, or
    master sub-aperture?
35. Are pixel-specific slant range and baseline values used?
36. Are actual orbit state vectors used, or is the baseline approximated from
    platform velocity and elapsed time?
37. How does one orbital trajectory provide elevation-wavenumber diversity?
38. What forward model demonstrates that the geometry resolves depth rather
    than surface motion projected into Doppler?
39. Has the inversion been tested on a conventional multi-baseline TomoSAR
    dataset with known elevations?

## 4. Propagation velocity and investigation frequency

40. How was propagation velocity `v` selected for each experiment?
41. Was it measured locally, taken from geological literature, or fitted using
    known structures?
42. Is one constant velocity used throughout each reconstructed volume?
43. How are layering, refraction, anisotropy, and spatially varying velocity
    handled?
44. Does `v` denote P-wave, S-wave, Rayleigh-wave, or another velocity?
45. What exactly is the “investigation frequency” `f`?
46. Is `f` a measured seismic vibration frequency or a processing parameter
    selected from SAR Doppler bandwidth?
47. How can an acquisition lasting seconds observe a physical vibration at
    12.5 or 22 kHz?
48. Are 12.5 and 22 kHz acoustic frequencies, SAR Doppler frequencies, or
    spatial-frequency parameters?
49. What mathematical relationship connects SAR Doppler frequency to
    subsurface acoustic frequency?
50. Why does the Giza preprint calculate `6000/12500 = 0.48 m`, while the
    published paper reports approximately `0.24 m`?
51. Did the published Giza calculation intentionally use `lambda=v/(2f)`?
52. Which convention is authoritative: `lambda=v/f` or `lambda=v/(2f)`?
53. Why does the patent state `v/f` while the published Giza number corresponds
    to `v/(2f)`?
54. Were the published depth axes independently calibrated, or calculated
    entirely from the selected `v` and `f`?

## 5. Physical forward model

55. What is the complete forward model connecting a buried acoustic source to
    a complex SAR surface pixel?
56. Does the radar observe surface displacement, velocity, acceleration,
    strain, or a Doppler-centroid perturbation?
57. How is subsurface propagation delay encoded in surface SAR phase or pixel
    displacement?
58. Where does source depth enter the measured signal before applying `A`?
59. Does the method require coherent acoustic phase propagation to the surface?
60. How are unknown source phase and excitation time handled?
61. How are multiple paths and reflections handled?
62. How are attenuation and geometric spreading included?
63. How are two depths producing the same surface vibration frequency
    distinguished?
64. How are depth and material velocity distinguished?
65. How are subsurface signals separated from wind, traffic, thermal motion,
    atmosphere, orbit error, and registration bias?
66. What null hypothesis was used to show that reconstructed structures were
    not processing artefacts?

## 6. Doppler sub-aperture construction

67. What exact master and slave Doppler-band definitions were used?
68. What were `B_CD`, `B_DL`, `N_D`, and `B_shift` in each experiment?
69. How was `B_shift` selected?
70. Was `B_shift` held fixed while the master/slave pair was swept?
71. Does every `Y_i` come from a master/slave pair at the same sweep position?
72. Were sub-apertures overlapping or non-overlapping?
73. What was the effective time interval between observations?
74. What was the temporal impulse response of each sub-aperture?
75. How was temporal aliasing handled?
76. How can frequencies beyond the pulse-history observation bandwidth be
    interpreted as measured physical vibrations?
77. Were the sub-aperture measurements treated as statistically independent?
78. What Doppler-filter window shape was used?
79. Were range-frequency sub-apertures used, or only azimuth Doppler bands?
80. Was processing performed from raw echoes, phase history, or focused SLC?

## 7. Pixel tracking and construction of `Y`

81. Which coregistration or pixel-tracking algorithm was used?
82. What search and correlation window sizes were used?
83. What sub-pixel estimator was used?
84. What displacement precision was measured on stationary controls?
85. How were range and azimuth pixel shifts converted to metres?
86. Is `Y_i = delta_x_i + j*delta_y_i` the intended construction?
87. Which shift component is real and which is imaginary?
88. Are shifts relative to one fixed master or to each paired master look?
89. Was interferometric phase used directly, or only intensity correlation?
90. How was phase unwrapping handled?
91. Were low-coherence observations rejected?
92. Was the constant displacement component removed?
93. Was bulk scene or orbital motion removed before forming `Y`?
94. What displacement and velocity noise floors were measured?

## 8. Resolution and unambiguous depth

95. How is `delta_T=lambda*R/(2*A)` derived from Eq. 22?
96. Why does that resolution contain no `t` if Eq. 22 contains `2*pi*t`?
97. Is `A` physical orbital distance, projected baseline span, or an equivalent
    aperture derived from Doppler bandwidth?
98. Why may the full along-track aperture be used as an elevation aperture?
99. What is the formula for maximum unambiguous depth?
100. How is depth aliasing detected?
101. What depth-grid spacing and extent were used in each publication?
102. Were displayed depth cells finer than the theoretical resolution?
103. What is the measured or simulated depth point-spread function?
104. Was depth resolution validated using a target at known depth?

## 9. Validation of the reconstructions

105. Which reconstructed features had independently known positions before
     processing?
106. Which discoveries were evaluated blindly rather than selected after
     viewing the tomograms?
107. What excavation, borehole, seismic, muographic, radar, or architectural
     ground truth was used?
108. What coordinate errors were obtained for known structures?
109. Were control regions without expected underground structures processed?
110. Were stationary real or synthetic controls passed through the identical
     chain?
111. Were time-order permutation tests performed?
112. Were phase-randomized or Doppler-band-shuffled null tests performed?
113. Were sensitivity sweeps over `v`, `f`, `B_perp`, and `B_shift` performed?
114. Do structures remain fixed when sub-aperture count and width change?
115. Do they remain fixed across acquisition time and viewing direction?
116. Do ascending and descending passes recover the same structures?
117. Were independent acquisitions processed without registration to expected
     structures?
118. What false-positive rate was measured?
119. What uncertainty accompanies each recovered depth?
120. Are displayed tomograms magnitude, power, normalized magnitude, or
     logarithmic magnitude?
121. What filtering, thresholding, interpolation, and enhancement were applied
     before publication?

## 10. Commercial X-band SAR products

122. Which Capella product level and format were used?
123. Did the Capella data contain pulse-level phase history or only focused
     SLC samples?
124. Which metadata fields were essential for timing and geometry?
125. Were Doppler-centroid, weighting, deskew, timing, and orbit parameters all
     available?
126. How were proprietary focusing choices or missing processing details
     handled?
127. Can vendor autofocus alter or suppress the micro-motion signature?
128. How do multilooking, quantization, radiometric processing, and phase
     flattening affect it?
129. Which Capella acquisition mode and dwell duration were used?
130. What prevented identical processing across Capella, Umbra, and
     TerraSAR-X?
131. Did the Capella experiment estimate frequency, time history, displacement,
     or tomographic depth?
132. Does the 2024 X-band paper support only surface micro-motion extraction,
     or also the depth-tomography stage?
133. What minimum signal access and metadata must a commercial product provide?

## 11. Reproducibility and errata

134. Can the original processing code be released or inspected?
135. Can one complete example provide `Y_i`, `t_i`, `B_perp,i`, `Kz_i`, the
     depth grid, `A`, and `h(z)`?
136. Can the exact parameter files for the Vesuvius, Giza, and patent examples
     be provided?
137. Can a small source-data crop and expected output be released as a
     regression fixture?
138. Which published equations directly describe the implementation, and which
     are conceptual?
139. Were undocumented corrections applied in the processing code?
140. Are there errata for Eq. 22, matrix orientation, wavelength convention, or
     the Giza arithmetic?
141. Can an authoritative corrected form of Eqs. 21–24 be provided with units,
     indices, and matrix dimensions?

## Minimum reproducibility package requested

A compact package answering the following would resolve most of the questions:

1. one measured `Y` vector;
2. the corresponding sub-aperture times and phase-centre coordinates;
3. the calculation of every `B_perp,i` and `Kz_i`;
4. the propagation velocity and investigation-frequency source;
5. the exact depth grid and steering matrix;
6. pseudoinverse and conditioning settings;
7. the resulting complex `h(z)` before display processing;
8. an independently known depth or a negative-control result.

Until those items are available, ResonarSat treats the conventional
`exp(j*Kz*z)` formulation as the default interpretation and the printed
`2*pi*t` factor as an explicitly labelled experimental sensitivity test.

---

## 11. Carrier frequency drift within a single dwell

Not a question for the sources -- a mechanism this project has never examined,
recorded so it is not mistaken for one that has been.

Bähr, *Orbital Effects in Spaceborne SAR Interferometry* (DGK Reihe C 719, KIT
2013), separates three error families that this project has been treating as
one. Section 3.4.2 classes a biased pulse repetition frequency as a **clock**
error rather than a timing error, and section 3.4.3 finds that "as long as
coregistration is implemented by amplitude cross-correlation, the
interferometric phase measurement is completely insensitive to errors in
`f_PRF` and `f_RSR`", whereas carrier frequency biases "can indeed produce
significant artefacts".

`RS_VALIDATE_PRF_STABILITY` measures the first of these. Its concern is
legitimate but different from Bähr's: this project does not coregister two
acquisitions, it uses pulse times to place sub-aperture centres on the time axis
of a spectrum, so a PRF error distorts that axis directly. The check stands.

**The error Bähr identifies as the one that corrupts phase is unmeasured here.**
He reports short-term ERS-1 carrier drifts of up to 82 Hz/s "that started and
stopped abruptly and lasted some tens of seconds", citing Massonnet and Vadon
(1995), and notes that such errors "were neither expected nor are they easy to
validate". Tens of seconds is the timescale of a long-dwell spotlight collect --
Giza is 32.869 s.

### It does not cancel

A pulse sent at slow time `t` carries `f0 + df(t)`, so a target at range `R`
returns phase `-4*pi*(f0+df(t))*R/c`. Comparing sub-look `m` against look 0:

```
dphi_m = -4*pi*R*[df(t_m) - df(t_0)] / c
```

Both looks come from one acquisition, but from different *instants* of it, and
the observable is exactly the drift difference across the lag. The common part
cancels; the part that accumulated between the two looks does not. As an
apparent line-of-sight displacement, for drift rate `f'` over lag `dt`:

```
d = lambda * R * f' * dt / c  =  8.193e-5 * f'[Hz/s] * dt[s]   metres   (Giza)
```

| drift | dt = 0.1 s | dt = 1 s | dt = 10 s |
|---|---|---|---|
| 1 Hz/s | 0.008 mm | 0.082 mm | 0.82 mm |
| 4 Hz/s | 0.033 mm | 0.328 mm | 3.28 mm |
| 82 Hz/s | 0.672 mm | 6.72 mm | 67.2 mm |

The phase CRLB floor at a 0.4 coherence is 0.329 mm per look, reached by a drift
of **4.02 Hz/s at a one-second lag** -- a twentieth of the short-term rate Bähr
reports. So the mechanism is not marginal at documented drift levels.

A *constant* drift rate gives a linear phase ramp and is removed by the
least-squares detrend in `rs_spectrum_compute()`. The dangerous case is the one
Bähr actually describes -- drifts "that started and stopped abruptly" -- which
leave a piecewise-linear residual that detrending does not remove and whose
energy sits in the lowest bins.

### It is not observable in this product's metadata

The Giza CPHD's PVP block carries `SC0` (word 29), `FX1` (21) and `FX2` (22) per
pulse. Read across all 335,149 vectors:

```
SC0  9.300000e9 Hz    min = max, spread 0.0000 Hz
FX1  9.300000e9 Hz    min = max, spread 0.0000 Hz
FX2  9.900000e9 Hz    min = max, spread 0.0000 Hz
```

Bit-identical, at exactly round values. That is a **declared nominal constant,
not a measurement** -- no oscillator produces the same value 335,149 times. The
fields are populated and carry no information about what the carrier actually
did.

(The PVP parse is verified against two independent quantities: `TxTime` spans
32.8686 s against the documented 32.869 s dwell, and `FX2 - FX1` = 600 MHz gives
`c/2B` = 0.2498 m range bins, the spacing `CLAUDE.md` records.)

### Where that leaves it

The mechanism is real, does not cancel, and is large enough at documented drift
rates to dominate the phase floor. It cannot be checked from this collect's
metadata, and estimating it from the data alone would mean separating a
low-frequency phase trend from exactly the low-frequency signal the method is
looking for -- which is the same separation problem the whole project keeps
meeting.

Bähr's remark that such errors "were neither expected nor are they easy to
validate" is borne out here. Nothing asserts this affects any run; it is
recorded as a mechanism that cannot currently be ruled out, and a reason not to
read a low-bin peak as a measurement.
