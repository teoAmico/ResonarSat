# SAR Data Sources

This document lists SAR datasets relevant to micro-motion and Doppler
sub-aperture analysis. Phase-preserving products are required:

- **CPHD** — compensated phase history
- **SICD/SLC** — focused single-look complex imagery

Detected products such as GEC, GEO, GRD, DGM and SIDD preserve amplitude but
are not substitutes for complex data.

## Available X-band open data

| provider | open complex products | access | licence | Giza coverage |
|---|---|---|---|---|
| **Capella** | CPHD, SICD, SLC | Public AWS bucket, no account | CC BY 4.0 | **Yes** |
| **Umbra** | CPHD, SICD | Public AWS bucket, no account | Creative Commons | None located |
| **ICEYE** | Primarily SLC COG; selected SICD/CPHD assets | Public AWS bucket, no account | CC BY 4.0 | None as of 2026-07-28 |

### Capella Open Data

Capella's public catalogue contains X-band Spotlight, Sliding Spotlight and
Stripmap collections. Available product families include CPHD, SICD, SLC, SIDD,
GEC and GEO.

Resources:

- AWS catalogue: https://registry.opendata.aws/capella_opendata/
- STAC catalogue: https://capella-open-data.s3.us-west-2.amazonaws.com/stac/catalog.json
- Access guide: https://support.capellaspace.com/how-do-i-access-capellas-open-data
- Product formats: https://support.capellaspace.com/what-sar-imagery-products-are-available-with-capella

#### Giza Spotlight CPHD

This public collection covers the Giza pyramid complex:

```text
CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012
Centre:        29.97500 N, 31.13068 E
Mode:          Spotlight
Polarization:  HH
Format:        CPHD 1.1
Dwell:         32.87 s
PRF:           10,196 Hz
Pulses:        335,149
Range bins:    25,073
File size:     39,180,270,944 bytes
Scene size:    approximately 5 km square
```

The scene centre is approximately 576 m from the Great Pyramid. Khufu, Khafre,
Menkaure and the Great Sphinx are inside the footprint.

Download:

```sh
curl -L -C - --retry 10 --retry-all-errors \
  -o data/giza.cphd \
  "https://capella-open-data.s3.us-west-2.amazonaws.com/data/2024/10/4/CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012/CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd"
```

Use `-C -` to resume interrupted transfers and verify the final byte count.

#### Cairo Spotlight SICD

Capella also provides a long-dwell scene over central Cairo:

```text
CAPELLA_C13_SP_SICD_HH_20241123062737_20241123062813
Centre:             30.04426 N, 31.23577 E
Mode:               Spotlight Ultra
Polarization:       HH
Orbit:              Descending
Observation side:   Left
Incidence angle:    43.4 degrees
Dwell:              35.02 s
Resolution:         0.05 m azimuth, 0.27 m range
Format:             SICD
```

Its footprint is 31.1993–31.2722 E and 30.0126–30.0759 N. It does not contain
the pyramids.

Resources:

- Collection: https://capella-open-data.s3.us-west-2.amazonaws.com/stac/capella-open-data-by-capital/capella-open-data-cairo/collection.json
- SICD item: https://capella-open-data.s3.us-west-2.amazonaws.com/stac/capella-open-data-by-datetime/capella-open-data-2024/capella-open-data-2024-11/capella-open-data-2024-11-23/CAPELLA_C13_SP_SICD_HH_20241123062737_20241123062813/CAPELLA_C13_SP_SICD_HH_20241123062737_20241123062813.json

#### Other long-dwell Capella locations

| centre | dwell | approximate CPHD size | location |
|---|---:|---:|---|
| 19.369 N, 155.196 W | 34–39 s | 36–39 GB | Kilauea, Hawaii |
| 45.976 N, 7.659 E | 33.1 s | 40.1 GB | Alps |
| 41.016 N, 28.967 E | 60.0 s | 60.8 GB | Istanbul |
| 39.480 N, 0.400 W | 60.0 s | 60.5 GB | Valencia |
| 10.738 S, 25.381 E | 30.3 s | 31.4 GB | Katanga copper belt |
| 48.115 N, 78.098 W | 29.4 s | 26.2 GB | Abitibi, Quebec |

### Umbra Open Data

Umbra provides X-band Spotlight collections in CPHD, SICD, SIDD and GEC.
The catalogue contains repeating observations of more than twenty locations
and assorted collections over many additional sites.

Open Spotlight CPHD collections inspected for this project generally have
2–8-second dwells.

Resources:

- Open Data Program: https://umbra.space/open-data/
- AWS catalogue: https://registry.opendata.aws/umbra-open-data/
- STAC catalogue: https://stacindex.org/catalogs/umbra-open-sar-data

Notable repeating locations include:

- Bingham Canyon copper mine, Utah
- Komati Power Station, South Africa
- Centerfield, Utah
- University of North Dakota
- Ports, mines and industrial facilities

### ICEYE Open Data

ICEYE's open X-band catalogue normally supplies:

- SLC COG
- GRD COG
- Quicklook imagery
- STAC/JSON metadata
- CSI and SAR-video products for selected Dwell collections

The public footprint catalogue contained no scene intersecting the Great
Pyramid coordinates (`31.1342 E, 29.9792 N`) when checked on 2026-07-28.
The catalogue is actively updated, so this result should be rechecked.

Resources:

- Open Data Initiative: https://www.iceye.com/open-data-initiative
- Open-data documentation: https://sar.iceye.com/6.0.6/opendata/opendata/
- STAC collection: https://iceye-open-data-catalog.s3.amazonaws.com/collections/iceye-sar.json
- Footprint summary: https://iceye-open-data-catalog.s3.amazonaws.com/stac-items/summary/footprint-summary-points.geojson
- Imaging modes: https://sar.iceye.com/latest/productspecification/imagingmodes/
- Product formats: https://sar.iceye.com/latest/productspecification/dataproducts/

## Corner reflectors: searched for, and not in any open archive

The one thing this project cannot get from an image alone is a target whose
motion is independently known. Every validated micro-motion result in the
literature uses a corner reflector, so the obvious question is whether any open
collect contains one.

**ESA's EDAP assessment of Capella** (EDAP.REP.072, Telespazio for ESA, 2022)
confirms the data exists. It characterises the Capella impulse response over the
**Rosamond calibration site in California**, which carries 38 permanently
installed corner reflectors -- 23 with a 2.4 m leg, 5 with 4.8 m, 10 with 0.7 m
-- and it names 29 acquisitions used, by full product ID and timestamp, such as
`CAPELLA_C03_SP_SLC_HH_20210409054258_20210409054302`.

None of it is public. Searched exhaustively on 2026-07-31:

| archive | items scanned | over the Rosamond box |
|---|---|---|
| Capella open data, CPHD | 1174 | **0** |
| Capella open data, SICD | 1504 | **0** |
| Umbra open data, all acquisitions | 2809 | **0** |

The box was generous, 34.65-35.05 N by 118.35-117.85 W. The named acquisitions
are not in the bucket either -- `data/2021/4/9/` and `data/2021/4/15/` are empty
-- which matches the assessment's own statement that the products were "provided
directly from the Capella team".

Rosamond is public infrastructure and its reflector positions are published, so
the obstacle is coverage rather than secrecy: neither open archive happens to
contain a collect over it. Getting one means a tasking request or an approved
data request, not a download.

### What the assessment is still worth

Three measured numbers that apply to the Capella products this project does use:

- **Azimuth sidelobes are high by design.** "Azimuth window is a rectangular
  window for SP and SS products resulting in a PSLR of about -13 dB", and "the
  azimuth antenna pattern is not compensated in the CAPELLA products". Measured
  on a spotlight acquisition: azimuth PSLR **-13.23 dB**, range PSLR -31.3 dB.
  That is a 5% amplitude sidelobe in azimuth on every Capella product, and it
  bears directly on the correlation-peak hopping documented in
  `runs/giza/2026-07-30-validated-spot-khufu/POSITIVE-CONTROL.md`.
- **Geolocation is good to about 1.5 m.** Absolute localisation error for
  spotlight: mean -0.73 m azimuth, 0.38 m range, standard deviations 1.5-1.6 m,
  never worse than 7.6 m. So a grid landing tens or hundreds of metres from its
  coordinates is this project's transform, not the product -- which is what the
  Capella SGN mirror in `readers/cphd.c` turned out to be.
- **Absolute radiometric calibration is suspect.** "The measured NESZ level is
  lower than the values reported in the products' annotations, suggesting a
  possible residual error in the absolute radiometric calibration." One more
  reason amplitudes from these products are labelled qualitative here.

## X-band data requiring approval or purchase

| provider | complex products | access |
|---|---|---|
| **COSMO-SkyMed** | SCS Level 1A | ASI institutional licence, ESA proposal or e-GEOS purchase |
| **TerraSAR-X** | Complex Spotlight and Staring Spotlight | DLR science proposal or commercial purchase |
| **Synspective StriX** | CEOS or SICD SLC | Commercial or research agreement |
| **PAZ** | Complex Spotlight and Stripmap | Commercial/research access |
| **KOMPSAT-5** | X-band complex products | Institutional or commercial access |
| **QPS-SAR** | X-band complex products | Commercial access |

### COSMO-SkyMed

The 2022 Giza paper identifies these COSMO-SkyMed Second Generation
acquisitions:

| date | geometry | reported beam/angle | polarization | use |
|---|---|---:|---|---|
| 28 October 2021 | right-descending | 06 | HH | external |
| 13 November 2021 | right-descending | 06 | HH | external |
| 27 October 2021 | right-descending | 08 | HH | external |
| 12 November 2021 | right-descending | 08 | HH | external |
| 24 July 2021 | right-ascending | 39 | HH | external |
| 9 August 2021 | right-ascending | 39 | HH | external |
| 25 February 2022 | left-descending | 46 | VV | internal |
| 16 November 2021 | right-descending | 48 | HH | internal |
| 22 February 2022 | right-descending | 48 | VV | internal |
| 16 February 2022 | right-ascending | 48 | VV | internal |

The paper's `46` and `48` labels should be confirmed against the catalogue;
official CSG Spotlight beam identifiers normally use names such as
`S2A-###`, `S2B-###` and `S2C-###`.

Access:

- ASI institutional portal: https://www.asi.it/en/earth-science/cosmo-skymed/
- ESA collection description: https://eocat.esa.int/eo-catalogue/collections/COSMO-SkyMed.full.archive.and.tasking
- ESA archive/tasking page: https://earth.esa.int/eogateway/catalog/cosmo-skymed-full-archive-and-tasking
- Archaeological research-access precedent: https://doi.org/10.3390/rs11111326
- Commercial catalogue: https://www.cleos.earth/

The COSMO-SkyMed entries exposed through Copernicus Data Space are specialized
sea-ice collections and do not provide general Giza coverage:

https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/missions/Cosmo-Sky-Med

### TerraSAR-X

DLR accepts proposals for scientific use of TerraSAR-X. Available modes include
Stripmap, ScanSAR, Spotlight and high-resolution Spotlight.

https://www.dlr.de/en/research-and-transfer/projects-and-missions/terrasar-x/data-access-and-products

### Synspective StriX

StriX operates at 9.65 GHz. SLC products are available in CEOS or SICD for
Stripmap, Sliding Spotlight and Staring Spotlight.

https://synspective.com/data/synspective-sar-data/

### PAZ

The immediately downloadable PAZ collection in Copernicus Data Space is a
sea-ice dataset dominated by detected GEC/MGD products. It is not a general
open complex archive.

https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/missions/paz

## Other open SAR data

These sources are not X-band but provide complex data useful for comparison.

### Sentinel-1

- Band: C-band
- Product: SLC
- Coverage: global archive
- Access: open

Resources:

- ASF Vertex: https://search.asf.alaska.edu
- Copernicus Data Space: https://dataspace.copernicus.eu

**IW cannot do the single-pass work, and the reason is dwell rather than
anything fixable in processing.** It is listed here for multi-pass comparison
only, and this section says why so that nobody spends a week rediscovering it.

TOPS steers the beam through each burst, so a target is illuminated for a
fraction of the burst. Working back from the product's own azimuth resolution --
`dx = lambda*R/(2*v*T)` with `dx` about 22 m at C band, `R` about 850 km --
gives an illumination time of **0.141 s** and therefore an observed aperture of
about **1.07 km**. The Giza Capella spotlight collect gives **32.87 s** and
**238.7 km**, a factor of 223.

Depth resolution is `lambda_ac*R/(2*A)` (`rs_tomo_resolution()`), so it inherits
that factor directly. At the same assumed constants for both -- `v = 1500 m/s`,
`f = 500 Hz` -- and computed with this project's own code:

| | illumination | aperture | **depth resolution** |
|---|---|---|---|
| Sentinel-1 IW | 0.141 s | 1.07 km | **594.6 m** |
| Capella spotlight (Giza) | 32.87 s | 238.7 km | **2.37 m** |

A resolution cell of 595 m cannot say anything about structures tens of metres
across, and no choice of `v` or `f` repairs it: both scale the two columns
identically, so the ratio is fixed by geometry. The unambiguous depth compounds
the point -- at 16 sub-apertures Sentinel-1's is 5.5 km, so the entire depth axis
falls inside one resolution cell and a half.

**Note what is NOT the problem.** The steering matrix is equally well conditioned
in both cases -- `rs_tomo_conditioning()` returns essentially the same condition
number and rank for either geometry, because conditioning depends on the depth
cell size relative to the resolution, not on the resolution itself. So a
Sentinel-1 inversion will run, converge, and produce a confident, well-posed,
structured profile. It will simply be resolving in units of 600 m. That is worth
stating because a failed inversion announces itself and a well-conditioned
meaningless one does not.

An independent attempt on exactly this data
(github.com/mfwarren/Pyramid, Sentinel-1A over Giza) reached the same conclusion
from the other direction, reporting an EM-style vertical resolution of about
4.6 km and acoustic-style figures of 108-215 m for its own assumed constants. Its
code is unlicensed, so nothing here derives from it; the numbers above are this
project's own.

### UAVSAR

- Band: L-band
- Platform: airborne
- Products: SLC, repeat-pass stacks and interferometric products
- Access: open

Resources:

- UAVSAR portal: https://uavsar.jpl.nasa.gov
- Delta-X SLC stack: https://daac.ornl.gov/DELTAX/guides/DeltaX_L1_UAVSAR_SLC_Stack.html
- NASA Earthdata entry: https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-deltax-l1-uavsar-slc-stack-1984-1.1
- Product formats: https://uavsar.jpl.nasa.gov/science/documents/rpi-format.html

## COSMO-SkyMed sample products

ESA publishes sample products suitable for inspecting the mission formats:

- First generation: https://earth.esa.int/eogateway/missions/cosmo-skymed/sample-data
- Second generation: https://earth.esa.int/eogateway/missions/cosmo-skymed-second-generation/sample-data
