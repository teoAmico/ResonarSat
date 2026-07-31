# Core questions

Six questions for the author of the method, Filippo Biondi, on the formulas as
printed in WO 2024/008365 A1, *Scanning Inside Volcanoes…* (RS 14(15):3828) and
*…Internal Structure of the Great Pyramid of Giza* (RS 14(20):5231).

They are requests to clarify the equations, not challenges to the results. Each
is a place where an independent implementation had to choose between readings,
and the choice changes the output.

Equation numbers follow the patent. For reference, the depth stage as printed is

$$
\mathbf{Y}=\mathbf{A}(K_z,z)\,\mathbf{h}(z)
\qquad\text{(Eq. 23)},
\qquad
\hat{\mathbf{h}}(z)=\mathbf{A}^{+}\mathbf{Y}
\qquad\text{(Eq. 24)},
$$

with $\mathbf{Y}$ the complex sub-pixel shift vector of Eq. 21, whose $i$-th
entry combines the two tracked components as

$$
Y_i=\Delta x_{\mathrm{az},i}+j\,\Delta x_{\mathrm{rg},i}.
$$

---

## 1. Equation 22: what is `t`?

The steering matrix is printed with elements

$$
A_{ij}=\exp\!\left(j\,2\pi\,K_{z,i}\,t\,z_j\right),
\qquad
K_{z,i}=\frac{4\pi B_{\perp,i}}{\lambda\,r_i\sin\theta}.
$$

`t` is not defined in the surrounding text, and the depth resolution derived from
the same matrix,

$$
\delta_T=\frac{\lambda R}{2A},
$$

contains no `t`.

- Is `t` one scalar for the whole matrix, or per-observation `t_i`?
- $K_z$ as defined already carries $4\pi$, so it is an angular spatial frequency
  in radians per metre. Is the leading $2\pi$ intentional, or redundant with it?
- Is the intended form simply

$$
A_{ij}=\exp\!\left(j\,K_{z,i}\,z_j\right)?
$$

---

## 2. Equation 23: what are the dimensions of `A`?

Eq. 21 gives `Y` as $k\times1$. Eq. 22 prints `A` with $k$ columns and $F$ rows.
Eq. 24 prints `h(z)` as a row, $1\times F$.

$\mathbf{Y}=\mathbf{A}\mathbf{h}$ is conformable only in one orientation:

$$
\underbrace{\mathbf{Y}}_{k\times1}
=\underbrace{\mathbf{A}}_{k\times F}\;
\underbrace{\mathbf{h}}_{F\times1}.
$$

Is that the intended orientation, with rows as sub-aperture observations and
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

$$
B_\parallel=\langle\vec B,\hat r_M\rangle,
\qquad
B_\perp=\sqrt{|\vec B|^2-B_\parallel^2},
$$

with $\hat r_M$ the unit line of sight

and the form actually used for height is the projection onto the cross-track
plane,

$$
B_\perp\approx B_h\cos\theta+B_v\sin\theta,
$$

stated explicitly **"assuming
B_a ~ 0"** -- `B_a` being the along-track component. Bähr notes the
simplification "can entail more or less significant biases when used for the
computation of the height ambiguity".

That assumption is not merely unmet for sub-aperture diversity; it is inverted.
A sub-aperture phase-centre separation is *entirely* along-track: `B = B_a e_a`
with no horizontal or vertical part. Under the zero-Doppler condition the line
of sight is perpendicular to the velocity, so `<e_a, r_M> = 0` and the scalar
form returns

$$
B_\parallel=0,
\qquad
B_\perp=|\vec B|=B_a .
$$

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

The standard TomoSAR vertical wavenumber

$$
K_z=\frac{4\pi B_\perp}{\lambda\,r_i\sin\theta}
$$

uses, in every independent derivation of it, the **radar carrier wavelength**
$\lambda_c=c/f_c$ -- Zhu & Bamler (*Very High Resolution
Spaceborne SAR Tomography in Urban Environment*, TGRS 48(12), 2010); Kim,
Villano & Krieger (*Volume Structure Retrieval Using Drone-Based SAR
Interferometry*, RS 16(8), 2024), who write the same expression directly in
terms of the radar frequency:

$$
K_z=\frac{4\pi B_\perp f_c}{c\,R\sin\theta}.
$$

The patent instead supplies a **mechanical** wavelength: an elastic wave speed
and a vibration frequency, $6600/22000=0.30$ m. The arithmetic is a correct
mechanical wavelength,

$$
\lambda_s=\frac{v_s}{f_v},
$$

but the question is why $\lambda_s$ stands in a formula derived
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

$$
N_D\,\Delta t=t_{\mathrm{dwell}}\frac{B_{CD}-B_{DL}}{B_{CD}}=t_{\mathrm{sap}}
$$

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

$$
f_{\mathrm{band}}=\frac{1}{2t_{\mathrm{sap}}}\;(\eta=0.5),
\qquad
f_1=\Delta f=\frac{1}{t_{\mathrm{sap}}}\;(\eta=1.0),
\qquad
\frac{f_1}{f_{\mathrm{band}}}=2 .
$$

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
difference across a fixed lag

$$
\Delta t_{\mathrm{pair}}=\frac{B_{\mathrm{shift}}\,t_{\mathrm{dwell}}}{B_{CD}},
$$

so the observable's response is

$$
|H(f)|=\left|2\sin\!\left(\pi f\,\Delta t_{\mathrm{pair}}\right)\right|,
$$

greatest at $f\,\Delta t_{\mathrm{pair}}=1/2$, which gives

$$
B_{\mathrm{shift}}^{\mathrm{opt}}(f)=\frac{B_{CD}}{2f\,t_{\mathrm{dwell}}},
$$

which is inversely proportional to `f`, as the patent says. The difficulty is
that the sweep geometry in the same document appears to cap it. `N_D` masters of
width `B_CD - B_DL` stepping by `(B_CD - B_DL)/N_D` leave one step of headroom,
so `B_shift <= (B_CD - B_DL)/N_D`. Combining:

$$
B_{\mathrm{shift}}^{\mathrm{opt}}(f)\le\frac{B_{CD}-B_{DL}}{N_D}
\iff
f\ge\frac{N_D}{2t_{\mathrm{sap}}}=N_D\,f_{\mathrm{band}} .
$$

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

$$
f\ge\frac{1.5}{L\,t_{\mathrm{dwell}}},
\qquad L=\frac{B_{DL}}{B_{CD}},
$$

which exceeds the band edge $f_{\mathrm{band}}=1/[2t_{\mathrm{dwell}}(1-L)]$ by

$$
\frac{3(1-L)}{L},
$$

a factor of 3 at the stated $L=1/2$, and below 1 only for $L>3/4$.

- Is there a selection rule for `B_shift` that was omitted?
- Is `B_shift` intended to be bounded by the sweep step, or may the sweep span
  shrink as `B_shift` grows?
- Under either reading, which frequencies is the parameter meant to select, given
  that both readings place them above the band the sub-apertures carry?
