# Run: 2026-07-30 giza / no-optimize-khufu

**Question this run is meant to answer:** Does the exhaustive correlation peak
search of `--no-optimize` change the null measured at Khufu by the patent chain
on 2026-07-29?

- git commit: `2581c91`
- started:    2026-07-30T13:43:14Z
- host:       Darwin arm64 (8 cores, 25.8 GB RAM), OpenMP enabled

## Why this run exists

`runs/giza/2026-07-29-patent-exact-true-khufu` measured a complete null at the
Great Pyramid with the patent chain: 220 of 225 windows exactly zero, none above
the quantisation floor, `RS_ERR_RANGE` rather than a reported frequency. That
result rests on a tracked shift, and the tracked shift rests on a peak search
that does **not** examine the whole correlation surface -- it refines within one
pixel of the integer peak.

`--no-optimize` replaces that with the global maximum of the entire zero-padded
surface. Per the flag's own notice this is *the only change in the mode that can
move a number*; the backprojection half is bitwise identical and the serial
tracking loop cannot matter because windows are independent. So this run isolates
exactly one question: was the null an artefact of where the peak search was
allowed to look?

This is a second measurement by a slower route, not a more trustworthy one.
Neither run passes a null test on its own.

## Collect

```
CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd  ->  data/giza.cphd
  Capella open data, CC BY 4.0
  dwell 32.869 s   335,141 pulses x 25,073 range bins   PRF 10,196.35 Hz
  carrier 9.3000 GHz   wavelength 0.0322 m
```

Same collect, same grid: `--at 29.979175,31.134186` -> `--offset -152,-552`,
Khufu.

## Commands

Byte-for-byte the 2026-07-29 command with `--no-optimize` added and the output
paths changed. Nothing else differs; a comparison is only worth making if it is
the only difference.

```sh
resonarsat mmotion --cphd data/giza.cphd --rbins 4096 \
    --at 29.979175,31.134186 --size 256 --cell 2.0 \
    --subap paper --reference pair --estimator correlation \
    --n 128 --win 32 --coherence 0 --no-detrend --upsample 40 \
    --no-optimize \
    --shifts runs/giza/2026-07-30-no-optimize-khufu/khufu_n128.csv \
    --out runs/giza/2026-07-30-no-optimize-khufu/khufu_n128
```

The tomographic half, likewise the 2026-07-29 cube's settings plus the flag.
Every value is read back off that cube's `.meta` sidecar rather than typed from
memory:

```sh
resonarsat tomo --cphd data/giza.cphd --at 29.979175,31.134186 \
    --size 256 --grid-cell 2.0 --n 128 --rbins 4096 \
    --velocity 6600 --frequency 22000 \
    --cell 1.3 --depth 37.8 \
    --patent-exact --no-optimize \
    --out runs/giza/2026-07-30-no-optimize-khufu/khufu_tomo.f32 \
    --geocode runs/giza/2026-07-30-no-optimize-khufu/khufu_geocoded.csv
```

**`--depth 37.8`, not the `--depth 60` in CLAUDE.md's worked example.** The
first attempt here used 60 and was refused:

```
tomo: value out of range (tomo: depth extent 60 m exceeds the 37.8252 m this
geometry represents without ambiguity at 128 sub-apertures
(v=6600, f=22000, A=238734, R=754216). Beyond it the profile folds.
Use --depth 37.8252 or fewer, or raise --n to about 204)
```

CLAUDE.md is not wrong; its example pairs `--depth 60` with `--n 512`, where 60 m
is inside the limit. The unambiguous extent scales with `--n`, and this run uses
`--n 128` to match 2026-07-29, which puts the ceiling at 37.8252 m. The two
values cannot be mixed between examples. The guard caught it, which is the
behaviour worth recording: a folded profile would have looked like a perfectly
ordinary cube.

## THE PREDICTION, recorded before the result

**The null holds, and the shifts come back all but identical.** The 2026-07-29
null was argued analytically: in paper mode `N_D * dt = t_sap` identically, so
`df = 1/t_sap` and every resolvable bin sits at an integer observation ratio,
where a displacement-averaging observable has no response. Nothing about that
argument involves the peak search, so widening the search cannot rescue a
frequency from it.

- **P1** `RS_ERR_RANGE` again -- no window resolves motion above the
  quantisation floor.
- **P2** `df = 0.0608 Hz`, `f_max = 3.89 Hz`, 65 bins, 225 windows -- all
  unchanged, since they follow from the dwell and `N_D` alone and the flag
  touches neither.
- **P3** The per-look shifts agree with 2026-07-29 in the large majority of the
  28,800 rows. Where they disagree it should be in windows whose correlation
  surface is flat, and the exhaustive search should give the *larger* excursion,
  not the smaller -- a distant sidelobe can win a global maximum that a local
  refinement would never have visited. Disagreement in that direction is
  expected and is not evidence of motion.

**What would falsify it:** the exhaustive search finding a peak clear of the
lowest bins, in a contiguous patch of windows, with non-zero excursion. That
would mean the null was an artefact of the restricted search rather than of the
chain's structure -- and it would call the 2026-07-29 result into question, not
confirm it.

## Result

### Every numeric output is byte-identical to 2026-07-29

Not "agrees within tolerance". Identical, by `cmp`:

| output | comparison |
|---|---|
| `khufu_n128.csv` (28,800 rows: 225 windows x 128 looks) | **byte-identical** |
| `khufu_tomo.f32` (15 x 15 x 30 float32 cube) | **byte-identical** |
| `khufu_geocoded.csv` (204 geocoded windows) | **byte-identical** |
| `khufu_tomogram_annotated.png` | **byte-identical** |

The only differences anywhere in the two runs are the provenance annotations the
flag is supposed to add, which is correct behaviour rather than a discrepancy:

```
hdr:   axes window_az, window_rg, depth [UNOPTIMIZED]
meta:  arithmetic_mode   [UNOPTIMIZED] exhaustive correlation peak search,
                         serial execution
meta:  measurement_chain [UNOPTIMIZED] subap=paper estimator=correlation ...
```

A cube from this mode cannot be mistaken for one from the optimised path, even
detached from this directory.

### The micro-motion half

```
sub-apertures: 128 looks, dt 0.1284 s
  observable band  f_max 3.89 Hz   AT sub-look resolution 0.10 m
sub-pixel refinement: 1/40 px (default 1/10 azimuth, 1/20 range)
tracked 225 windows (15 x 15); 225 pass the 0.00 coherence mask
spectra: 65 bins, 0.0608 Hz resolution
mmotion: value out of range (spectrum: no window resolved motion above the
  0.061225 px quantisation floor (3 sigma of 1/40 px))
```

All three predictions held:

| | predicted | measured |
|---|---|---|
| **P1** | `RS_ERR_RANGE`, no window above the floor | **`RS_ERR_RANGE`, 0 of 225 above the floor** |
| **P2** | `df = 0.0608 Hz`, `f_max = 3.89 Hz`, 65 bins, 225 windows | **all four unchanged** |
| **P3** | shifts agree in the large majority of 28,800 rows | **agree in 28,800 of 28,800** |

P3 was stated with a clause describing where disagreement was expected and which
direction it should take. There was no disagreement to characterise.

### What this establishes, and what it does not

**Establishes:** the null is not an artefact of where the peak search was allowed
to look. The exhaustive global maximum over the whole zero-padded surface
selected exactly the same peak as the +/-1 px refinement, in all 225 windows at
all 128 looks. This retires the peak search as a candidate explanation.

That the agreement is *exact* rather than close is the informative part, and it
is what the local-refinement assumption predicts when it holds. The patches are
mean-removed but untapered, and `--patent-exact` forces rectangular sub-aperture
filters, so each peak carries sinc sidelobes at -13.3 dB -- which can never
outrank their own main lobe. A global maximum can therefore only differ from the
local one where the surface has **two comparable main lobes**, i.e. no dominant
match. Nowhere on this scene does it. See `docs/SUBPIXEL-PEAK-LITERATURE.md`
(local, untracked) for the supporting literature.

**Does not establish:** anything independent about the null itself. Per the
flag's own notice, backprojection under `--no-optimize` is bitwise identical by
construction, so the two runs share every stage upstream of the correlation. This
is one measurement checked along one axis, not two measurements agreeing. The
2026-07-29 conclusion stands on its own analytic argument -- paper mode forces
`N_D * dt = t_sap`, so every resolvable bin sits at an integer observation ratio
where a displacement-averaging observable has no response -- and nothing here
strengthens or weakens that.

**Not a passed null test.** Neither run has one. The standing project position is
unchanged: no micro-motion measurement made with this implementation has passed
its null test, and none of this should be presented as a demonstrated
sensitivity.

### Incidental finding: the depth-ambiguity guard fired

Reconstructing the tomo command by taking `--depth 60` from CLAUDE.md's example
while keeping `--n 128` from the 2026-07-29 cube produced a refusal, not a cube.
The unambiguous extent scales with `--n`; the two examples cannot be mixed. The
correct value came from reading the existing cube's `.meta` sidecar
(`depth_max_m 37.8`, `%g`-printed, so exactly 37.8) rather than from the docs.
Worth noting because a folded profile would have looked like a perfectly
ordinary cube.
