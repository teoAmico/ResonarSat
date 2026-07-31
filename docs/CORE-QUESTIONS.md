# Core questions

Six questions for the author of the method, Filippo Biondi, on the formulas as
printed in WO 2024/008365 A1, *Scanning Inside Volcanoes…* (RS 14(15):3828) and
*…Internal Structure of the Great Pyramid of Giza* (RS 14(20):5231).

They are requests to clarify the equations, not challenges to the results. Each
is a place where an independent implementation had to choose between readings,
and the choice changes the output.

---

## 1. Equation 22: what is `t`?

The steering matrix is printed as

```
A(K_z, z) = [ 1, exp(j*2*pi*K_z2*t*z_0), ..., exp(j*2*pi*K_z(k-1)*t*z_0) ]
            [ ...                                                        ]
```

`t` is not defined in the surrounding text, and `delta_T = lambda*R/(2*A)`,
derived from the same matrix, contains no `t`.

- Is `t` one scalar for the whole matrix, or per-observation `t_i`?
- `K_z = 4*pi*B_perp/(lambda*r_i*sin(theta))` is already in radians per metre.
  Is the leading `2*pi` intentional, or redundant with it?
- Is the intended form simply `exp(j*K_z,i*z_j)`?

---

## 2. Equation 23: what are the dimensions of `A`?

Eq. 21 gives `Y` as `k x 1`. Eq. 22 prints `A` with `k` columns and `F` rows.
Eq. 24 prints `h(z)` as a row.

`Y = A h(z)` is conformable only if `A` is `k x F` and `h(z)` is `F x 1`. Is
that the intended orientation, with rows as sub-aperture observations and
columns as depth cells?

---

## 3. Is `B_perp` an orthogonal baseline or the orbit aperture?

Eq. 22 defines `K_z = 4*pi*B_perp/(lambda*r_i*sin(theta))`, with `B_perp` "the
`i`-th **orthogonal baseline**".

One paragraph later the resolution is `delta_T = lambda*R/(2*A)`, with `A` "the
**orbit aperture** considered in the tomographic synthesis … proportional to the
Doppler bandwidth used to synthesize the sub-apertures".

An orbit aperture is along-track; an orthogonal baseline is perpendicular to the
line of sight.

- Which of the two enters `K_z`?
- If it is the along-track sub-aperture phase-centre separation, what projection
  converts it to a perpendicular component?

**What the standard decomposition says.** Bähr, *Orbital Effects in Spaceborne
SAR Interferometry* (DGK Reihe C 719 / KIT 2013, refereed by Hanssen), §3.1
defines the baseline as the difference of two sensor positions, `B = x_S - x_M`,
and decomposes it two ways. The general scalar form is

```
B_par  = <B, r_M>                    r_M the unit line-of-sight
B_perp = sqrt(|B|^2 - B_par^2)
```

and the form actually used for height is the projection onto the cross-track
plane, `B_perp ~ B_h cos(theta) + B_v sin(theta)`, stated explicitly **"assuming
B_a ~ 0"** -- `B_a` being the along-track component. Bähr notes the
simplification "can entail more or less significant biases when used for the
computation of the height ambiguity".

That assumption is not merely unmet for sub-aperture diversity; it is inverted.
A sub-aperture phase-centre separation is *entirely* along-track: `B = B_a e_a`
with no horizontal or vertical part. Under the zero-Doppler condition the line
of sight is perpendicular to the velocity, so `<e_a, r_M> = 0` and the scalar
form returns

```
B_par  = 0
B_perp = |B| = B_a
```

The scalar `B_perp` is therefore **not zero** -- it equals the whole along-track
separation. But the vector it measures points along track, so it produces no
height parallax, while `K_z` and `delta_T` will accept the number and return a
finite, ordinary-looking depth axis from it. This is the failure mode the rest of
this project keeps meeting: a plausible result from a configuration that cannot
produce it.

The independent literature points the same way. Azimuth sub-aperture pairs are
used in Multiple-Aperture InSAR (Bechor & Zebker 2006) to measure **along-track
displacement**, which is what an along-track baseline is sensitive to. We are
not aware of a use of azimuth sub-aperture separation as a height baseline.

- Is `B_perp` in Eq. 22 the scalar `sqrt(|B|^2 - B_par^2)`, or the cross-track
  projection that the height ambiguity conventionally uses?
- If the former, what supplies the height sensitivity, given that the separation
  is along-track and produces no parallax across the line of sight?

---

## 4. Why does a mechanical wavelength appear where the radar one belongs?

The factor-of-two question below is the smaller half of this. The prior question
is which wave `lambda` refers to at all.

`K_z = 4*pi*B_perp/(lambda*r_i*sin(theta))` is the standard TomoSAR vertical
wavenumber, and in every independent derivation of it `lambda` is the **radar
carrier wavelength** `lambda_c = c/f_c` -- Zhu & Bamler (*Very High Resolution
Spaceborne SAR Tomography in Urban Environment*, TGRS 48(12), 2010); Kim,
Villano & Krieger (*Volume Structure Retrieval Using Drone-Based SAR
Interferometry*, RS 16(8), 2024), who write the same expression directly as
`4*pi*B_perp*f_c/(c*R*sin(theta))`.

The patent instead supplies a **mechanical** wavelength: an elastic wave speed
and a vibration frequency, `6600/22000 = 0.30 m`. The arithmetic is a correct
`lambda_s = v_s/f_v`; the question is why `lambda_s` stands in a formula derived
with `lambda_c`. The micro-Doppler literature keeps them strictly apart --
Clemente & Soraghan carry radar phase through `lambda_c` and the mechanical
motion separately as `A_v*cos(omega_v*eta)`, coupled but never interchangeable.

**This compounds with question 3.** `K_z` then carries two substitutions at
once: an along-track separation standing in for an elevation-sensitive
baseline, and a mechanical wavelength standing in for the radar one. Either
alone removes the expression's derivation.

- What derivation licenses `lambda_s` in place of `lambda_c` here? At minimum it
  would need the coupling from a subsurface elastic mode to surface
  displacement, from that displacement to SAR phase, and from mechanical
  frequency to a spatial depth wavenumber.
- Which convention is authoritative for the wavelength itself? The patent states
  `lambda = v/f` and works `6600/22000 = 0.30 m`; the Giza preprint computes
  `6000/12500 = 0.48 m`; the published Giza paper reports approximately `0.24 m`,
  which is `v/(2f)`. The answer is a factor of two on every depth.
- Is `f` a measured vibration frequency or a processing parameter? A
  sub-aperture sequence over a 33 s dwell samples surface motion at a few hertz,
  not at 12.5 or 22 kHz.

## 5. Does `N_D * dt = t_sap` follow as intended?

With the master–slave pair stepped by `(B_CD - B_DL)/N_D` across a span of
`B_CD - B_DL`, the record length is

```
N_D * dt = t_dwell * (B_CD - B_DL) / B_CD = t_sap
```

for every `N_D` and every `B_DL`. Then `df = 1/t_sap`, and spectral bin `k` sits
at `f_k = k/t_sap` — an integer observation ratio, where a sub-aperture
integrates a whole number of cycles of the motion.

- Is that the intended construction?
- If so, what makes a displacement-averaging observable non-zero at exactly
  those frequencies?

The same identity also fixes where the observable band ends, and the two land in
a way worth asking about directly. A sub-aperture averages the motion over its
own duration, so the series carries nothing above `1/(2*t_sap)` however finely
the steps sample it. Since `df = 1/t_sap`, that band edge sits at exactly half
the first bin:

```
band edge  = 1/(2*t_sap)        (observation ratio 0.5)
first bin  = 1/t_sap            (observation ratio 1.0)
```

The ratio is **2 for every `N_D` and every `B_DL`** -- it does not depend on the
dwell, the collect or the held-out fraction, so no choice of parameter moves it.
The lowest frequency the layout can report is therefore twice the highest one it
can carry, and it falls precisely on the first averaging null.

Measured on a synthetic collect (20 s dwell, `B_DL = B_CD/2`, `N_D = 128`): the
chain reports 0.100 Hz against a 0.050 Hz band, with the sub-aperture response
at -240 dB and the observation ratio at exactly 1.00.

- Is the intended sweep span `B_CD - B_DL`, as printed, or something wider?
- If as printed, what is meant to be read from a spectrum whose every bin sits
  on an averaging null?

---

## 6. What sets `B_shift`, and can it reach the observable band?

WO 2024/008365 A1 [0004] holds the two bands "rigidly held at a distance
`B_shift`", chosen to select "the precise vibrational frequency one wishes to
observe", and adds that the higher it is, the lower the frequency observed. No
value or selection rule is given, here or in the claims; a best value is stated
only for the held-out band, `B_CL = B_CD/2`.

The stated property implies a rule. Each pair sample is a displacement
difference across a lag `dt = B_shift * t_dwell / B_CD`, so its response is
`|2 sin(pi f dt)|`, greatest at `f*dt = 1/2`:

```
B_shift_opt(f) = B_CD / (2 * f * t_dwell)
```

which is inversely proportional to `f`, as the patent says. The difficulty is
that the sweep geometry in the same document appears to cap it. `N_D` masters of
width `B_CD - B_DL` stepping by `(B_CD - B_DL)/N_D` leave one step of headroom,
so `B_shift <= (B_CD - B_DL)/N_D`. Combining:

```
B_shift_opt(f) <= step   <=>   f >= N_D / (2 * t_sap)   =   N_D * (band edge)
```

**The lowest frequency `B_shift` can be tuned to is `N_D` times the highest
frequency the layout can observe** -- again independent of dwell and collect.
Measured: an implementation refuses `B_shift` above 0.0621 Hz for `N_D = 128` on
a 15.90 Hz Doppler band, where tuning to that collect's 0.050 Hz band edge would
need 7.95 Hz, a factor of 128.

Reading `B_shift` as independent of the sweep -- which the patent's wording
allows, since it advances the rigid pair in `N_D` steps separately -- lets the
sweep span shrink to buy separation. The two then compete for the same `B_CD *
B_DL/B_CD` of headroom: span sets the frequency resolution, separation sets the
lag. Asking for both a resolvable bin and the optimal lag gives

```
f >= 1.5 / (L * t_dwell),    L = B_DL/B_CD
```

which exceeds the band edge `1/(2*t_dwell*(1-L))` by `3*(1-L)/L` -- a factor of
3 at the stated `L = 1/2`, and below 1 only for `L > 3/4`.

- Is there a selection rule for `B_shift` that was omitted?
- Is `B_shift` intended to be bounded by the sweep step, or may the sweep span
  shrink as `B_shift` grows?
- Under either reading, which frequencies is the parameter meant to select, given
  that both readings place them above the band the sub-apertures carry?
