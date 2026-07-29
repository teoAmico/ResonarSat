# Contributing to ResonarSat

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

No external dependencies are required — the FFT is vendored. Before opening a
pull request, also run the sanitiser build, which is what CI gates on:

```sh
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=ASAN
cmake --build build-asan --parallel
ctest --test-dir build-asan --output-on-failure
```

## Every function needs a block comment

Not a style preference — a hard requirement, enforced by CI
(`scripts/check_docs.py`), including for `static` helpers.

The house style is Redis's: a `/* ... */` block immediately above the function
**definition in the `.c` file** (not the header, so the contract lives in one
place and cannot drift), prose in complete sentences, no Doxygen tags.

State **units and conventions on every argument**: Hz against rad/s, slant
against ground range, samples against metres, azimuth line rate against transmit
PRF, electromagnetic against acoustic wavelength. Nearly every defect found in
this codebase and in the literature it reproduces has been a units or convention
error, not a logic error. Also cover what the function returns, the error
contract, memory ownership, preconditions, cost when it is not obvious, and any
assumption or known limitation. The last of those is the most valuable and the
most often skipped.

When implementing a published relation, cite it by section, equation or figure
number so a reader can check the line against the source. Cite published,
citable sources only — no private prior work, no local filesystem paths.

## Errors

Every fallible function returns `resonarsat_status_t` and describes the failure
through `rs_set_error()`. Readers parse untrusted external binary formats: a
truncated, corrupt or hostile file must produce a described status, never a
crash, a partial write into a caller's struct, or a silently zero-filled result.
Validate any field that determines an allocation or a file offset *before* using
it.

## Changes to the tomography stage

For changes to Phase 5 (`src/core/tomo.c` and its interfaces), please include
two results in the pull-request description:

1. **A synthetic depth test** showing how the inversion recovers an injected
   structure.
2. **A parameter-sensitivity sweep** showing how recovered depths move as
   the assumed velocity and investigation frequency vary.

These results make tomography changes easier to review and compare. The depth
axis depends on assumed constants, so a synthetic recovery test checks the
implementation while the sweep makes the effect of those assumptions visible.
Plots, tables or the relevant test output are all suitable.

The current interface follows two conventions:

- `--velocity` and `--frequency` are required rather than defaulted, so the
  chosen depth scale is explicit.
- Model A keeps the measured vibration frequency separate from its kHz
  investigation-frequency parameter. They represent different quantities and
  substituting one for the other changes the depth resolution by orders of
  magnitude.

Proposals to change either convention are welcome. Include the source or
derivation for the new interpretation, update the command documentation, and
extend the synthetic and sensitivity tests so reviewers can evaluate the change
on the same evidence as the existing implementation.

## Reporting results

When reporting ResonarSat results, state the assumed constants used for the
depth axis, the wavelength convention, and the along-track tomographic
baseline. The README and each tomogram's metadata sidecar contain these values
to make that reporting straightforward.
