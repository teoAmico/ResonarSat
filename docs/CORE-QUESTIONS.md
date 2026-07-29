# Core questions

Five questions for the author of the method, Filippo Biondi, on the formulas as
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

---

## 4. Is `lambda = v/f` or `v/(2f)`, and what is `f`?

The patent states `lambda = v/f` and works it as `6600/22000 = 0.30 m`. The Giza
preprint computes `6000/12500 = 0.48 m`. The published Giza paper reports
approximately `0.24 m`, which is `v/(2f)`.

- Which convention is authoritative? The answer is a factor of two on every
  depth.
- Is `f` a measured vibration frequency or a processing parameter? A
  sub-aperture sequence over a 33 s dwell samples surface motion at a few hertz,
  not at 12.5 or 22 kHz.

---

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
