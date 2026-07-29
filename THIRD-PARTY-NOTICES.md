# Third-party notices

ResonarSat itself is MIT licensed (see [`LICENSE`](LICENSE)). It vendors and links
the components below. All are permissively licensed and compatible with
distributing ResonarSat under MIT, including in commercial and closed-source
products.

Keep this file current. It is the first thing a corporate legal review asks for,
and reconstructing it after the fact is far more work than maintaining it.

## Vendored

### PocketFFT

- **Location:** `third_party/pocketfft/`
- **Upstream:** https://github.com/mreineck/pocketfft
- **License:** BSD-3-Clause — full text in `third_party/pocketfft/LICENSE.md`
- **Copyright:** © 2010-2019 Max-Planck-Society
- **Used for:** every Fourier transform in the project, behind the `rs_fft`
  wrapper in `src/core/fft.c`.

FFTW was evaluated and **rejected**: its standard build is GPL licensed, and
linking it into an MIT-licensed binary would make the distributed work
GPL-encumbered while `LICENSE` promised otherwise. That is a licensing defect
and a misrepresentation to downstream users, so the FFT backend is BSD by
requirement rather than by preference. See §2 and §6 of the implementation plan.

## Build and runtime dependencies

None required. The project builds with a C11 compiler and CMake alone.

OpenMP is used when the toolchain provides it and is detected optionally; every
parallel region is guarded and the code is correct single-threaded. Apple clang
ships without it by default, which is supported.

## Planned, not yet integrated

These will be needed by readers not yet implemented. Their licences were checked
before selection so that adding them cannot change ResonarSat's licensing
position:

| Component | Purpose | License |
|---|---|---|
| GDAL | Sentinel-1 GeoTIFF read, GeoTIFF write | MIT-style |
| HDF5 | COSMO-SkyMed SCS products | BSD-style |
| libxml2 | Sentinel-1 and SICD metadata | MIT |
| OpenBLAS / LAPACKE | Model C inversion, sparse spectral estimation | BSD |
| libpng | PNG quicklooks (PGM is used today) | libpng license |
| NITRO | SICD and CPHD NITF containers | verify per fork — NITF library licensing has varied |
