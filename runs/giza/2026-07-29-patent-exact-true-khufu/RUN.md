# Run: 2026-07-29 giza / patent-exact-true-khufu

**Question this run is meant to answer:** Does the patent's own chain still
measure nothing when pointed at the ACTUAL Great Pyramid, now that the range
mirror is corrected?

- started:    2026-07-29T19:35:50Z
- host:       Darwin arm64 (8 cores, 25.8 GB RAM)

## Why this run exists

This reader mirrored Capella imagery in range about the scene reference point, so
a grid placed by coordinates landed on ground that was not the target -- about a
kilometre out at Khufu. A null measured there is not a statement about the
pyramid; it has to be measured at the pyramid to be one.

The mirror is now fixed: Capella declares `SGN = +1` while shipping the opposite
convention, and the reader inverts the FX-to-delay transform for their products,
keyed on `CollectorName`.

## Collect

```
CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd  ->  data/giza.cphd
  Capella open data, CC BY 4.0
  dwell 32.869 s   335,141 pulses x 25,073 range bins   PRF 10,196.35 Hz
  carrier 9.3000 GHz   wavelength 0.0322 m
```

**Khufu by coordinates, not by a hand-derived offset.** `--at
29.979175,31.134186` resolves to `--offset -152,-552`, and a focused image there
shows the Great Pyramid centred with the mastaba field's tomb rows around it. A
hand-typed offset is not good enough for this target; the coordinate resolution
is what puts the grid on the monument.

## Command

```sh
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 256 --cell 2.0 \
    --subap paper --reference pair --estimator correlation \
    --n 128 --win 32 --coherence 0 --no-detrend --upsample 40 \
    --shifts runs/giza/2026-07-29-patent-exact-true-khufu/khufu_n128.csv \
    --out runs/giza/2026-07-29-patent-exact-true-khufu/khufu_n128
```

`--rbins 4096` gives +/-512 m of slant range, enough for Khufu's 552 m
cross-track offset plus the grid's own half-extent. This is the same chain
`--patent-exact` forces in `tomo`: the paper front end, the swept master-slave
pair, and the coregistrator shift observable.

## THE PREDICTION, recorded before the result

**The same null, and for a reason that does not depend on pointing.** In paper
mode the sweep step is `(B_CD - B_DL)/N_D`, so `N_D * dt = t_sap` identically,
`df = 1/t_sap`, and every resolvable bin sits at an INTEGER observation ratio --
a frequency at which a displacement-averaging observable has no response.
`--patent-exact` also forces a rectangular sub-aperture filter, which makes
those nulls exact rather than merely deep.

That argument is analytic. It says the chain cannot report a vibration at any
frequency it can resolve, wherever it is aimed. So:

- **P1** `RS_ERR_RANGE`, or a peak one or two bins above DC.
- **P2** The tracked master-slave offset is at or near zero in essentially every
  window.
- **P3** `df = 0.0608 Hz` and `f_max = 3.89 Hz`, unchanged, since both follow
  from the dwell and `N_D` alone.

**What would falsify it:** a peak well clear of the lowest bins, in a contiguous
patch of windows, with a non-zero excursion. That would mean the earlier null
was an artefact of pointing at the wrong ground after all, and the structural
argument is wrong.

## Result

### The null holds at the Great Pyramid itself

```
sub-apertures: 128 looks, dt 0.1284 s
  observable band  f_max 3.89 Hz   AT sub-look resolution 0.10 m
tracked 225 windows (15 x 15); 225 pass the 0.00 coherence mask
spectra: 65 bins, 0.0608 Hz resolution
mmotion: value out of range (spectrum: no window resolved motion above the
  0.061225 px quantisation floor)
```

All three predictions held:

| | predicted | measured |
|---|---|---|
| **P1** | `RS_ERR_RANGE`, or a peak 1-2 bins above DC | **`RS_ERR_RANGE`** |
| **P2** | offset at or near zero in essentially every window | **220 of 225 exactly zero; 225 of 225 below the floor** |
| **P3** | `df = 0.0608 Hz`, `f_max = 3.89 Hz` | **0.0608 Hz, 3.89 Hz** |

The largest excursion anywhere is 2.00 quantisation steps -- 0.05 px -- against a
floor of 2.449 steps. Nothing clears it, so the chain reports no frequency at
all rather than a fabricated one.

### Pointing correctly changed nothing for the conclusion

**Five of the 225 windows move at all**, by one or two quantisation steps, so the
Great Pyramid is marginally less inert than bare desert. That is the whole of the
difference, and it was predicted rather than hoped for.

The reason is analytic and has nothing to do with aim: the paper sweep makes
`N_D * dt = t_sap` identically, so `df = 1/t_sap` and every resolvable bin sits
at an integer observation ratio, where a displacement-averaging observable has no
response. `--patent-exact`'s rectangular filter makes those nulls exact. The
chain cannot report a vibration at any frequency it can resolve, wherever it is
aimed -- and now that has been measured on the monument the claim is about.

### Status of the falsification test

The prediction named what would overturn it: a peak clear of the lowest bins, in
a contiguous patch of windows, with non-zero excursion. None of the three
occurred. There is no peak, no patch, and no excursion above the floor.

**This is the version of the result that should be cited**, because it is the one
measured on the Great Pyramid itself.
