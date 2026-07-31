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
| **Umbra** | CPHD, SICD, SIDD, CSI — **1,413 of 1,556 items carry CPHD** | Public AWS bucket, no account | Creative Commons | None located |
| **ICEYE** | Primarily SLC COG; **6 of 374 items carry CPHD + SICD** | Public AWS bucket, no account | CC BY 4.0 | None as of 2026-07-28 |

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
Range bins:    29,160        (Data/Channel/NumSamples)
Image area:    25,073 x 25,073 at 0.1994 m  (SceneCoordinates/IAXExtent)
File size:     39,180,270,944 bytes
Scene size:    approximately 5 km square
```

An earlier version of this block gave the range-bin count as 25,073, which is
the **image-area line count**, not the signal array width. The two sit in
different parts of the CPHD XML and both appear as `NumSamples`. The signal
array is 29,160 wide, confirmed three ways: `Data/Channel/NumSamples`, SarPy's
`data_size`, and the signal block itself -- 39,091,779,360 bytes at 4 bytes per
CI4 sample over 335,149 vectors is exactly 29,160. Caught by cross-reading the
file with SarPy.

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

That table listed the long-dwell collects known at the time. A full scan of the
catalogue on 2026-07-31 found the archive is far larger, and the survey below
replaces it.

**Every CPHD item in the catalogue was fetched and its dwell taken from the
start/end timestamps in the product ID.** Of 1174:

| dwell | collects | unique sites (clustered to 0.2 deg) |
|---|---:|---:|
| >= 20 s | **531** | — |
| >= 25 s | 379 | — |
| >= 30 s | 175 | **93** |
| >= 40 s | 8 | 6 |
| 60 s | 3 | 2 (Istanbul, Valencia) |

The longest are 60.0 s: Istanbul once and Valencia twice. Nothing in the open
CPHD archive exceeds a minute.

Sites at >= 30 s dwell over built infrastructure, which is what the micro-motion
stage needs. Incidence is listed because the tracker's amplitude floor scales as
`1 / cos(incidence)` for vertical motion, so a steep look costs sensitivity:

| dwell | n | inc | site | first ID |
|---:|---:|---:|---|---|
| 60.0 s | 1 | 38.5 | **Istanbul** | `CAPELLA_C11_SP_CPHD_HH_20230907223849_20230907223949` |
| 60.0 s | 2 | 36.1 | **Valencia** | `CAPELLA_C09_SP_CPHD_HH_20240227111009_20240227111109` |
| 41.0 s | 1 | 52.6 | Abu Dhabi | `CAPELLA_C13_SP_CPHD_HH_20241121060730_20241121060811` |
| 41.0 s | 1 | 53.8 | Addis Ababa | `CAPELLA_C14_SP_CPHD_HH_20241121180053_20241121180134` |
| 40.0 s | 3 | 56.1 | **Budapest** (Danube) | `CAPELLA_C15_SP_CPHD_HH_20241115212743_20241115212823` |
| 40.0 s | 1 | 52.5 | Algiers | `CAPELLA_C14_SP_CPHD_HH_20241121212131_20241121212211` |
| 39.0 s | 22 | 55.5 | Kilauea, Hawaii | `CAPELLA_C13_SP_CPHD_HH_20241003024738_20241003024817` |
| 39.0 s | 1 | 48.6 | **Venice** | `CAPELLA_C13_SP_CPHD_HH_20240816153121_20240816153200` |
| 39.0 s | 2 | 56.2 | Dalian | `CAPELLA_C14_SP_CPHD_HH_20240705020414_20240705020453` |
| 39.0 s | 2 | 56.6 | Ningbo | `CAPELLA_C14_SP_CPHD_HH_20240902210947_20240902211026` |
| 38.0 s | 3 | 54.9 | **Paris** | `CAPELLA_C14_SP_CPHD_HH_20240730213056_20240730213134` |
| 38.0 s | 1 | 50.3 | **Bangkok** | `CAPELLA_C14_SP_CPHD_HH_20241121131023_20241121131101` |
| 37.0 s | 1 | 43.3 | **Amsterdam** | `CAPELLA_C15_SP_CPHD_HH_20241122215521_20241122215558` |
| 37.0 s | 2 | 54.1 | Chicago | `CAPELLA_C14_SP_CPHD_HH_20240613040046_20240613040123` |
| 34.0 s | 1 | **37.5** | **Rome** (Tiber) | `CAPELLA_C13_SP_CPHD_HH_20240816102624_20240816102658` |
| 33.0 s | 1 | **37.8** | Bergen County, NJ (Hudson) | `CAPELLA_C14_SP_CPHD_VV_20240514125342_20240514125415` |
| 33.0 s | 1 | 49.4 | Toronto | `CAPELLA_C10_SP_CPHD_HH_20231225194617_20231225194650` |
| 32.0 s | 1 | **28.6** | North Tyneside, UK | `CAPELLA_C11_SP_CPHD_VV_20240514051800_20240514051832` |
| 32.0 s | 1 | 31.0 | Astana | `CAPELLA_C11_SP_CPHD_HH_20241121035026_20241121035058` |
| 30.0 s | 1 | 46.2 | Shanghai (Yangpu) | `CAPELLA_C10_SP_CPHD_HH_20240518171634_20240518171704` |
| 30.0 s | 1 | 37.0 | Ankara | `CAPELLA_C10_SP_CPHD_HH_20241124182803_20241124182833` |
| 33.0 s | 2 | 38.6 | *Giza, for reference* | `CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012` |

Place names are reverse-geocoded scene centres, not verified footprint contents;
only Istanbul below has been checked against what is actually inside the box. The
remaining ~70 sites at >= 30 s are terrain -- Morocco (many), Mexico, the Alaska
North Slope, Congo, Myanmar, Indonesia, assorted islands. Useful for the
tomography side, not for a structure with a mode.

**Reproducing this survey.** The catalogue is plain STAC over HTTPS and needs no
account. `stac/capella-open-data-by-product-type/capella-open-data-cphd/collection.json`
lists every item by relative href; each item is a ~4 KB JSON carrying `bbox`,
`view:incidence_angle` and `sar:instrument_mode`. Dwell is not a field -- take it
from the two timestamps in the product ID. Fetching all 1174 in parallel takes
about 30 seconds.

**Reading a collect's geometry without downloading it.** A CPHD begins with an
ASCII header giving `XML_BLOCK_BYTE_OFFSET` and `XML_BLOCK_SIZE`, and the XML
block that follows carries `ReferenceGeometry` -- slant range, ARP velocity,
incidence angle -- plus `NumVectors`, `TxTime2` and the `FxBand`. Two HTTP range
requests totalling about 12 KB give everything needed to compute a collect's
amplitude window, against a 17-61 GB download. Do this before fetching anything.

#### The SICD collection is larger and the CLI cannot use it

The same scan over the 1504 SICD items finds **633 collects at >= 20 s dwell
across 339 unique sites**, including 60 s scenes over Istanbul and Valencia and
40 s over Budapest. SICD is focused imagery, so a stack comes from spectral
splitting rather than the pulse route -- `rs_subaperture_split()`, which
`subaperture.h` notes is the route the published method actually describes.

**But `--sicd` is accepted only by `resonarsat info`.** `focus`, `mmotion`,
`tomo`, `sweep` and `validate` all read `--cphd` and nothing else, so none of
those 633 collects can currently be measured. `rs_read_sicd()` exists and
`rs_subaperture_split()` takes an `rs_slc_t`, so the gap is CLI wiring rather
than missing science. Closing it widens the usable archive considerably and cuts
per-collect download sizes by an order of magnitude. It also gives up the pulse
route, which is the one the only working positive control used, so a SICD run is
not a substitute for a CPHD run until the split route is shown to track
something.

### Umbra Open Data

Umbra provides X-band Spotlight collections in CPHD, SICD, SIDD and GEC.
The catalogue contains repeating observations of more than twenty locations
and assorted collections over many additional sites.

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

#### Corrected: this is the largest open CPHD archive, not a starter dataset

An earlier version of this section said Umbra's open Spotlight CPHDs "generally
have 2-8-second dwells." That describes the median and misses the entire tail,
and it is why Umbra sat here as a reader-development starter rather than a
primary source. The full STAC tree was walked on 2026-07-31 -- **1,556 items,
all SPOTLIGHT**:

| | |
|---|---|
| **items carrying a `.cphd` asset** | **1,413 of 1,556 (91%)** -- more than Capella's 1,174 |
| dwell | min 1.0 s, median 6.4 s, **max 42.58 s** |
| >= 20 s dwell with CPHD | **40** |
| >= 25 s | 28 |
| >= 30 s | 9 |

Each collect also ships SICD, SIDD, a CSI image and an `ANALYSIS.zip`.

The 40 long-dwell CPHD collects fall over just seven sites:

| collects | dwell | azimuth res | site |
|---:|---|---|---|
| **28** | 20.6-30.0 s | 0.062-0.063 m | **Tippecanoe County, Indiana** (West Lafayette) |
| 8 | 24.2-42.6 s | 0.062-0.50 m | North Slope Borough, Alaska (Kivalina area) |
| 1 | 34.4 s | **0.050 m** | Tucson, Arizona |
| 1 | 25.4 s | 0.125 m | Jeddah, Saudi Arabia |
| 1 | 23.6 s | 0.125 m | Orchard Park, New York |
| 1 | 23.0 s | 0.125 m | Campo, Cameroon |

#### The Tippecanoe series is a control substrate, not a target

Twenty-eight CPHD collects of the same ~6.6 km box (40.4562-40.5154 N,
87.0372-86.9596 W) at 20.6-30.0 s dwell and 6.2-6.3 cm azimuth resolution,
incidence 32-43 degrees, spanning January to August 2025.

Checked against OpenStreetMap, the box holds **no significant structure**: a rail
bridge on the Lafayette Line, a trunk-road bridge on Sagamore Parkway, and five
residential road crossings. Nothing with a mode worth measuring.

**That is what makes it valuable.** The hardest-won conclusion in
`runs/giza/2026-07-30-validated-spot-khufu/POSITIVE-CONTROL.md` is that a claim
about a *configuration* needs the sweep repeated over independent realisations
and the verdicts pooled -- seed 7 scored slope 1.004 where seed 23 reported
6.275 Hz against an injected 1.1 Hz from the same configuration. That criterion
was built for synthetic data and has had no real-data equivalent.

**Twenty-eight acquisitions of the same ordinary ground is that equivalent**, and
it needs no ground truth whatsoever:

- Run one configuration across all 28 and pool the reported dominant
  frequencies. A configuration that answers differently on each pass is
  measuring itself, and that conclusion follows without knowing what is on the
  ground.
- Because the box holds nothing that should vibrate at millimetre scale, a
  confident spectral peak is a **false positive measured on real clutter and
  real geometry** -- a quantity this project has never had. Every false-positive
  figure it currently quotes comes from `rs_simulate_static_like()` or from
  synthetic fixtures.
- The repeats span eight months, so seasonal and atmospheric variation is
  sampled rather than assumed.

An uninteresting scene is the right scene for this. It sits alongside
`--null-static` rather than replacing it: the simulated null shares the real
geometry but not the real clutter, and this shares both.

Its one limitation is dwell class. At 20-30 s it is Giza-like, so it inherits
the aperture-fraction problem discussed under the Istanbul section, and the
`t_sap` a usable window demands will again be far below the published 4.5-7.6%.

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

#### Six open ICEYE CPHD collects, and the best target set found so far

Surveyed 2026-07-31: all **374 items** in the open STAC collection were fetched.
Modes are spotlight 320, stripmap 44, scan 9. Asset sets:

| assets | items |
|---|---:|
| grd, qlk, slc, csi, vid | 258 |
| grd, qlk, slc | 100 |
| grd, qlk only | 9 |
| **grd, qlk, slc, csi, vid + `cphd`, `sicd`, `sidd`** | **6** |

So six collects carry phase history. Geometry below is from each CPHD's own XML
header, read by range request:

| site | centre | product | dwell | `df` | inc | R | size |
|---|---|---|---:|---:|---:|---:|---:|
| **Bratislava** | 48.138 N, 17.105 E | spot-fine | **15.574 s** | 0.0642 Hz | 40.87 | 760.8 km | 16.5 GB |
| Houston (Texas Medical Center) | 29.685 N, 95.410 W | dwell-precise | 15.345 s | 0.0652 Hz | **25.21** | 580.1 km | 19.4 GB |
| Vienna | 48.202 N, 16.333 E | dwell-precise | 15.010 s | 0.0666 Hz | 27.61 | 600.6 km | 15.0 GB |
| Paris (Le Bourget) | 48.960 N, 2.437 E | spot-fine | 14.3 s | 0.070 Hz | 30.25 | -- |
| Mexico City (south) | 19.303 N, 99.150 W | dwell-fine | 17.4 s | 0.057 Hz | 36.94 | -- |
| Vandenberg / Lompoc | 34.633 N, 120.614 W | dwell-precise | 14.7 s | 0.068 Hz | 26.20 | -- |

**ICEYE's STAC `start_datetime`/`end_datetime` IS the aperture time**, unlike
TerraSAR-X's. Verified against the files: STAC 15.574 s against `TxTime2`
15.574 s at Bratislava, and likewise at Vienna and Houston. The two conventions
are indistinguishable in a catalogue listing, so check one file per mission
before trusting the field.

##### Bratislava

Its footprint (17.0619-17.1471 E, 48.1098-48.1669 N) contains **five Danube
crossings of four structural types**, verified against OpenStreetMap:

| structure | position | type |
|---|---|---|
| **Most SNP** (UFO Bridge) | 48.1383, 17.1046 | **`bridge:structure=cable-stayed`**, single pylon |
| Starý most | 48.1371, 17.1171 | truss; tram, cycle and foot |
| Most Apollo | 48.1358, 17.1279 | arch, primary road |
| Prístavný most | 48.1350, 17.1404 | beam; motorway and rail |
| Most Lanfranconi | 48.1412, 17.0752 | beam, motorway |

Four structural types in one scene is an internal control the single-target
scenes do not offer: a beam viaduct and a cable-stayed span should not carry the
same modal signature, and if they do, that is the processing rather than the
ground.

Its window, by the bounds in the Istanbul section: floor `2.960/f` mm, open
below `t_sap` 0.833 s, 4.4x margin below 0.189 s -- so `alpha` 1.22%, 82 looks,
sub-look 8.2 m, and an admissible **2.80-12.32 mm at f <= 1.06 Hz**.

Against the Istanbul 25 s candidate:

| | Bratislava | Istanbul 25 s |
|---|---|---|
| dwell / `df` | 15.6 s / 0.064 Hz | 25.0 s / **0.040 Hz** |
| incidence | 40.87 deg | **19.53 deg** |
| floor | **2.960/f mm** | 2.673/f mm |
| window | 2.80-12.32 mm | 2.24-9.84 mm |
| `alpha` for the 4.4x margin | **1.22%** | 0.67% |
| looks | 82 | **149** |
| size | 16.5 GB | 17.3 GB |
| targets | **5 bridges, 4 types** | 3 crossings |

Istanbul wins on incidence, frequency resolution and look count; Bratislava on
the target set, and on an `alpha` nearly double and therefore closer to the
published range. Istanbul first on the numbers, Bratislava second and arguably
the better scientific run, because five bridges give controls one scene cannot.

##### The SGN override will not fire, and that is a trap

All three ICEYE headers read declare **`SGN = -1`**; Capella declares `+1` and
ships the opposite, which is why `rs_read_cphd()` mirrors it. That override is a
substring match on `CollectorName`, which here reads `ICEYE-X50`, so **it will
not fire.** Whether ICEYE honours its declared sign is untested.

This is exactly the failure mode item 3 of `FOLLOW-UPS.md` records: "if another
vendor ships the same defect it is not caught." The first ICEYE run must focus
and look for range mirroring before any measurement, and the data-driven check
proposed in that item -- compressing a pulse both ways and choosing the direction
whose energy lands inside the declared `TOA1`/`TOA2` support -- would settle it
per product without a vendor name.

##### An aside worth knowing

258 of the 374 items ship `csi-cog` (Colorized Sub-Aperture Image) and
`vid-mp4`. ICEYE productises sub-aperture decomposition -- the same operation
`rs_subaperture_split()` performs. Their CSI is a three-look colorisation rather
than a time series, so it validates nothing here, but it is a useful external
reference for what the decomposition of a real scene should look like.

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

Since no reflector is obtainable, the next question is whether any open collect
holds a target whose ordinary motion is inside what the tracker can measure. The
next section is the first one found.

## Istanbul: a footprint whose targets could be inside the instrument's window

No corner reflector exists in any open archive, so the remaining question is
whether any open collect contains a target whose *ordinary* motion falls inside
what the tracker can measure. This section records the first one found, and the
arithmetic that qualifies it.

### Why the window, and not the dwell, is the selection criterion

`src/core/validate.c` bounds the tracked azimuth excursion from both sides. At a
0.4 m cell, the cell size where the floor was measured:

```
floor   (vertical amplitude)  a_min = 3.5 * cell * V / (R * 2*pi*f * cos(inc))
ceiling                       a_max = 0.75 * lambda * R / (2*V*t_sap) / (R*2*pi*f*cos(inc)/V)
window width                  a_max / a_min  =  0.75 * lambda * R / (2*V*t_sap) / (3.5*cell)
```

**Which wavelength.** Every window figure in this document is computed at the
**band-centre** wavelength, `c/((FxMin+FxMax)/2)`, read from the product header.
`readers/cphd.c:455` deliberately uses something else -- the per-pulse `SC0`,
the first sample's frequency -- because the FX-to-delay transform leaves a
residual `exp(-j*2*pi*SC0*tau)` that the backprojector must undo, and using the
band centre there would leave about 3.4 radians of phase error per metre of
offset. Both are right for their own purpose, and cross-checking against SarPy
(which reports the band) is what surfaced the difference.

The consequence is that `resonarsat validate` will report ceilings slightly
higher than the tables here, since `res_sap` scales with lambda: **+2.7% for the
Istanbul collects and +6.5% for the ICEYE ones**, whose transmit bands are wider.
Nothing below changes -- no collect crosses a threshold -- but treat the quoted
ceilings as conservative by that much, and take the authoritative figure from
`validate` run on the file.

Three consequences drive everything below:

- **The window's width depends only on `t_sap`.** Frequency slides the window up
  and down; it cannot open it. So the selection question is how short the
  sub-aperture must be, not what frequency is sought.
- **The window does not depend on dwell at all.** A longer collect buys
  frequency resolution `df = 1/T` and more looks, nothing else. A 60 s collect
  and a 25 s collect of the same geometry measure the same amplitudes.
- **Higher frequency is strictly better,** because the floor goes as `1/f`.

Giza, for comparison, has floor `2.750/f` mm and ceiling `2.491/(t_sap*f)` mm,
its window is open below `t_sap` 0.906 s, and reaches the 4.4x margin that is the
only ratio ever demonstrated to work below `t_sap` 0.206 s -- an aperture
fraction of 0.63% against the 4.5-7.6% the published campaigns use. See
`runs/giza/2026-07-30-validated-spot-khufu/POSITIVE-CONTROL.md` for where the
7 px floor and the 4.4x margin come from, and for the fact that a ratio of 1.2
was measured *closed*.

### The two Istanbul collects

Every CPHD item in the catalogue was scanned; **exactly two are over Istanbul.**
Geometry below is from each product's own XML header, read by HTTP range request
without downloading the signal block.

| | 60 s collect | 25 s collect |
|---|---|---|
| ID | `CAPELLA_C11_SP_CPHD_HH_20230907223849_20230907223949` | `CAPELLA_C09_SP_CPHD_HH_20230321101754_20230321101819` |
| size | 60,792,882,864 B (60.8 GB) | **17,312,033,344 B (17.3 GB)** |
| acquired | 2023-09-07 22:39 UTC = **01:39 local** | 2023-03-21 10:18 UTC = **13:18 local** |
| dwell / PRF | 60.00 s / 10,152 Hz | 25.00 s / 9,977 Hz |
| `df` | 0.0167 Hz | 0.0400 Hz |
| slant range | 811.0 km | 636.5 km |
| platform speed | 7230 m/s | 7197 m/s |
| incidence | 38.55 deg | **19.53 deg** |
| wavelength | 0.0311 m | 0.0311 m |
| CPHD version | 1.1.0 | 1.0.1 |
| bbox | 28.9257-29.0076 E, 40.9854-41.0474 N | 28.9422-29.0120 E, 40.9790-41.0319 N |

`readers/cphd.c` checks only the `CPHD/` magic and not the version, so the 1.0.1
product should read, but that is untested.

### Their windows

| | floor | ceiling | open below | 4.4x margin below |
|---|---|---|---|---|
| 25 s | `2.673/f` mm | `1.967/(t_sap*f)` mm | `t_sap` 0.736 s | 0.167 s |
| 60 s | `2.540/f` mm | `2.371/(t_sap*f)` mm | `t_sap` 0.933 s | 0.212 s |
| *Giza* | *`2.750/f` mm* | *`2.491/(t_sap*f)` mm* | *0.906 s* | *0.206 s* |

Concrete operating points at the conservative `eta = 0.2`:

```
25 s   t_sap 0.167 s  (alpha 0.67%, 149 looks, sub-look  8.2 m)  f<=1.20 Hz  ->  2.24 -  9.84 mm
25 s   t_sap 0.084 s  (alpha 0.33%, 298 looks, sub-look 16.4 m)  f<=2.39 Hz  ->  1.12 -  9.84 mm
60 s   t_sap 0.212 s  (alpha 0.35%, 282 looks, sub-look  8.2 m)  f<=0.94 Hz  ->  2.69 - 11.85 mm
60 s   t_sap 0.106 s  (alpha 0.18%, 565 looks, sub-look 16.4 m)  f<=1.89 Hz  ->  1.35 - 11.85 mm
```

The window is essentially identical to Giza's -- same sensor class, similar
orbit. **The instrument's reach does not change with the collect.** What changes
is whether the target lives inside it.

### What is actually in the footprints

Verified against OpenStreetMap rather than inferred from the scene centre: 23
named bridge features in the 60 s box, 10 in the 25 s. Present in **both**:

| structure | position | type |
|---|---|---|
| **M2 Yenikapi-Haciosman metro bridge** over the Golden Horn | 41.0227 N, 28.9667 E | `bridge:structure=cable-stayed`, `railway=subway` |
| **Galata Koprusu** | 41.0200 N, 28.9733 E | arch/bascule, road + T1 tram |
| **Ataturk Koprusu** | 41.0242 N, 28.9651 E | beam viaduct, primary road |
| **Marmaray** Sirkeci-Kazlicesme rail line | 41.0048 N, 28.9536 E | 16 / 11 bridge segments |
| **Bozdogan Su Kemeri** (Valens Aqueduct) | 41.0162 N, 28.9552 E | Roman masonry arch |

In the 60 s box only: **Halic Koprusu** (O-1 motorway viaduct, 8 segments,
41.0436 N, 28.9420 E), plus T4 tram and Metrobus busway bridges.

The Valens Aqueduct matters as much as the movers: a massive masonry arch in the
same scene, at the same range, through the same processing, is a *static
reference* the null tests can use without simulating anything.

### Which collect, and why the shorter one

**The 25 s daytime collect first**, despite half the dwell:

1. **13:18 local.** Traffic, metro and trams running. The 60 s collect is 01:39
   local -- the M2 metro is not operating and road traffic is minimal. Traffic
   is the principal excitation for these decks, and a 60 s dwell over an
   unexcited bridge measures nothing.
2. **Incidence 19.5 deg against 38.5 deg**, so cos 0.943 against 0.782 -- about
   20% better projection of vertical modes onto the line of sight.
3. **17.3 GB against 60.8 GB.**
4. The extra dwell buys only `df`, and 0.040 Hz is ample to place a mode in a
   1-3 Hz band.

The cost is losing Halic Koprusu, the motorway bridge. The cable-stayed metro
bridge, Galata and Ataturk are all retained.

### What this does and does not establish

**Does not:** close test 1 of `IMPLEMENTATION-VERIFICATION.md`. There is still no
accelerometer, no reflector, no synchronous truth. Published modal frequencies
for these bridges, if any exist, would be weak external truth -- more than Khufu
offers, and less than validation.

**Does:** put a target class whose plausible amplitudes overlap the instrument's
window in front of this software for the first time. The relevant comparison is
the Trento validation in `MODIFIED-BACKPROJECTION.md`: a corner reflector driven
at **1.5 cm at 2 Hz**, which sits at the top edge of the window computed above.
That is the point. The window is engineered-target territory, and a bridge deck
under traffic is the nearest thing to it that occurs naturally.

Four caveats to carry into any run:

- **2.2-9.8 mm is still large for ambient deck response.** Traffic-induced
  mid-span deflection on a medium-span bridge is typically sub-millimetre to a
  few millimetres. Plausible, not comfortable. This is the first target where the
  two ranges touch at all, against Khufu's three-orders-of-magnitude gap.
- **The required `alpha` of 0.33-0.67% is an order of magnitude below the
  published 4.5-7.6%,** with 8-16 m sub-looks. `validate` warns there for a
  reason: the target may stop being a distinct feature in each sub-look. Nothing
  in the literature or in this project's tests has been run in that regime.
- **A bridge over water is a bright structure on a dark background.** High
  contrast, but `microm.c` warns that isolated targets on empty scenes score low
  even when tracking perfectly. A deck is extended rather than a point, so
  windows on it should carry structure -- an assumption, not a measurement.
- **Ships in the Marmara** are visible in the 25 s scene. They translate rather
  than oscillate, and are worth using as a "does the chain see motion at all"
  sanity check, but they are not vibration.

### Order of operations

1. `validate` the geometry before committing to the download -- the header read
   above already gives every input it needs.
2. Seed a run directory with `tools/new-run.sh` and write the question into
   `RUN.md` **before** fetching anything: *does a bridge deck in the Istanbul
   daytime collect produce a frequency inside the configuration's own amplitude
   window?*
3. Focus and look, per the grid-placement rule -- a misplaced grid over a river
   is as invisible as a misplaced grid over desert.
4. Run with `--null-static`, and use the Valens Aqueduct as the in-scene static
   reference.

## X-band data requiring approval or purchase

| provider | complex products | access |
|---|---|---|
| **COSMO-SkyMed** | SCS Level 1A | ASI institutional licence, ESA proposal or e-GEOS purchase |
| **TerraSAR-X** | Complex Spotlight and Staring Spotlight | DLR science proposal or commercial purchase |
| **Synspective StriX** | CEOS or SICD SLC | Commercial or research agreement |
| **PAZ** | Complex Spotlight and Stripmap | Commercial/research access |
| **KOMPSAT-5** | X-band complex products | Institutional or commercial access |
| **QPS-SAR** | X-band complex products | Commercial access |
| **ASNARO-2** (Tellus) | L1.1 SLC | Free with account, **but cannot leave their cloud** |

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
Stripmap, ScanSAR, Spotlight, high-resolution Spotlight and Staring Spotlight.

https://www.dlr.de/en/research-and-transfer/projects-and-missions/terrasar-x/data-access-and-products

#### The DLR Geohazard Supersites archive

Surveyed 2026-07-31 via the open STAC API at
`geoservice.dlr.de/eoc/ogc/stac/v1`. Of the **71 collections** that portal
serves, exactly one carries SAR complex data. Everything else is derived --
TanDEM-X DEMs, forest/non-forest, World Settlement Footprint, snow, atmospheric
L3, and `SAR4TEC` (Sentinel-1-derived tectonic products, CC BY 4.0 but not
phase-preserving).

The exception is `SUPERSITES`, "TerraSAR-X Geohazard Supersites Products",
licence `proprietary`, temporal extent from 2008-05-02. **10,656 items.** Of
6,000 sampled:

| | |
|---|---|
| product type | **94.5% SSC** (single-look slant-range complex), 5.2% GEC, 0.3% MGD |
| sensor mode | SL 2753, SM 2544, **ST 494**, HS 209 |
| sites | LatinAmerica 2034, SE-Asia 1495, Ecuador 705, Iceland 587, NorthPacific 403, **Vesuv 289**, NewZealand 182, **Marmara 177**, Hawaii 68, Kahramanmaras 38 |

So this is genuine SLC-grade complex data, unlike the rest of the portal.

**The catalogue is open; the data is not.** STAC browsing needs no account, and
each item carries a direct `download.geoservice.dlr.de` href. Every one of them
302s to `sso.eoc.dlr.de/eoc/auth/login` -- tested on all four modes. Access is by
Principal Investigator on a Supersite region under an accepted TerraSAR-X Science
proposal, not registration.

#### `start_datetime` is the scene transit, not the aperture time

The trap worth carrying away from this archive. Reported durations by mode:

| mode | reported duration | x ground velocity ~6.7 km/s | matches |
|---|---|---|---|
| ST | 0.389-0.441 s | 2.6 km | ST scene azimuth extent |
| HS | 0.7-0.8 s | 4.7-5.4 km | HS scene |
| SL | 1.4-1.5 s | 9.4-10 km | SL scene |
| SM | 3.0-8.14 s | 20-55 km | stripmap scene |

All four agree, so `start_datetime`/`end_datetime` is the **scene ground-extent
transit**. Read as a dwell it inverts the ranking exactly backwards: staring
spotlight, the only mode with a usable illumination time, reports the *shortest*
number in the catalogue, and stripmap the longest. This is the Sentinel-1 IW
error in a new disguise -- see the IW section below, where the same distinction
between scene time and target illumination time decides the mission.

Target illumination follows from the mode's azimuth resolution,
`t = lambda*R / (2*V*res_az)`, at R ~ 650 km, V ~ 7600 m/s, lambda 0.0311 m:

| mode | azimuth resolution | illumination time |
|---|---|---|
| **ST** staring spotlight | 0.24 m | **~5.5 s** |
| HS high-resolution spotlight | 1.1 m | ~1.2 s |
| SL spotlight | 2.0 m | ~0.7 s |
| SM stripmap | 3.3 m | ~0.4 s |

Only ST is in contention, and it is 6x shorter than the 32.9 s Giza collect and
11x shorter than the 60 s Istanbul one.

#### What ST would have bought, and why it still fails

At TerraSAR-X geometry on a 0.4 m cell, by the bounds in the Istanbul section
above: floor `3.18/f` mm, ceiling `2.27/(t_sap*f)` mm, window open below `t_sap`
0.71 s and reaching the 4.4x margin below `t_sap` 0.162 s. At `eta` 0.2 that is
`f <= 1.23 Hz` and an admissible **2.6-11.4 mm** -- the same window as Capella,
which is one more confirmation that the window is a property of the sensor class
rather than of the collect.

Two differences from Capella, one in each direction:

- **`alpha` = 0.162/5.5 = 2.9%**, against 0.63% at Giza and 0.67% at Istanbul,
  and against the 4.5-7.6% the published campaigns use. Independent confirmation
  from a second sensor that **`t_sap` in seconds is the transferable quantity and
  `alpha` is not**: Giza demands an absurd aperture fraction only because its
  dwell is long.
- **`df` = 0.18 Hz**, against 0.040 Hz for the Istanbul 25 s collect. A 1.5 Hz
  mode lands in bin 8. Coarse, but not disqualifying on its own.

**The archive fails on coverage, not on physics.** Staring spotlight exists only
over LatinAmerica (192), SE-Asia (126), NorthPacific (95) and NewZealand (81) --
geohazard regions, no infrastructure. The two sites that look promising by name
both disappoint:

- **Vesuv: 289 items, all SM.** Worth stating plainly, because `docs/` carries
  Biondi's *Scanning inside volcanoes with SAR echography*: this archive holds
  289 TerraSAR-X SSC scenes over that paper's exact target, at ~0.4 s
  illumination, which cannot support the micro-motion stage at any
  configuration.
- **Marmara: 18 HS + 159 SM**, no ST.

So: gated, and behind the gate the one useful mode does not point at a useful
target. Not a candidate.

#### The one angle worth keeping

The corner-reflector section above concludes that obtaining a Rosamond collect
"means a tasking request or an approved data request, not a download." DLR's
TerraSAR-X Science Service System (https://sss.terrasar-x.dlr.de/) is exactly
that mechanism, and a science proposal can request **new tasking** rather than
archive access alone.

TerraSAR-X ST over Rosamond, or over a bridge with an instrumented deck, is a
concrete route to the collect-with-ground-truth that test 1 of
[`IMPLEMENTATION-VERIFICATION.md`](IMPLEMENTATION-VERIFICATION.md) needs and that
no open archive contains. It is a months-long administrative path rather than a
download, so it changes nothing about what to do next -- but it is the only route
to closing that test found so far, and 2.9% is an aperture fraction the published
literature would recognise.

### Synspective StriX

StriX operates at 9.65 GHz. SLC products are available in CEOS or SICD for
Stripmap, Sliding Spotlight and Staring Spotlight.

https://synspective.com/data/synspective-sar-data/

### PAZ

The immediately downloadable PAZ collection in Copernicus Data Space is a
sea-ice dataset dominated by detected GEC/MGD products. It is not a general
open complex archive.

https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/missions/paz

### ASNARO-2, via Tellus

Checked 2026-07-31. Two datasets, and the pair is instructive because each fails
for a different reason.

**The urban-area stripmap dataset is detected.** Level 1.5 and 2.1, so multilook
detected and orthorectified amplitude with no phase -- excluded by the rule at
the top of this document. 2 m resolution, 12 km swath, six Japanese cities
(Sapporo, Toyama, Nagoya, Kumamoto, Kagoshima, Ube, Nakatsu), 2021-2022,
purchase required.
https://www.tellusxdp.com/en-us/catalog/data/asnaro-2_x-band_sar_images_urban_area_dataset_stripmap.html

**The L1.1 dataset is genuinely complex and still unusable.** It is SLC
"including phase information", SPOTLIGHT-1/2 at 1.0 m, STRIPMAP at 2.0 m,
SCANSAR at 16 m, CEOS format, Japan/Asia/Oceania, 2018-2021, and **free with a
Tellus account**. The blocker is not the licence fee but the environment: "this
can only be used in the Tellus environment", and the terms forbid copying the
original out, "including just a part of, or formatted data that can be restored
back to the original". ResonarSat is a local binary with its own readers, so
data that cannot leave a hosted notebook is not usable by it at all.
https://www.tellusxdp.com/en-us/catalog/data/asnaro-2_l1_1.html

**And the dwell would settle it anyway.** By the same derivation used for the
TerraSAR-X modes above, at lambda 0.0312 m (9.6 GHz), R ~ 623 km, V ~ 7610 m/s:

| mode | azimuth resolution | illumination time |
|---|---|---|
| SPOTLIGHT-1/2 | 1.0 m | **~1.3 s** |
| STRIPMAP | 2.0 m | ~0.6 s |
| SCANSAR | 16 m | ~0.08 s |

At 1.28 s the frequency resolution is **0.78 Hz**: bins at 0, 0.78, 1.56 and
2.34 Hz, four of them across the entire band a structure's modes live in. And
the 4.4x window margin needs `t_sap` <= 0.156 s, which over a 1.28 s dwell is
**eight looks** at zero overlap. That is not a spectrum.

#### What the dwell comparison across four sensors shows

| source | usable dwell | `df` | `alpha` for the 4.4x margin |
|---|---:|---:|---:|
| Capella Istanbul 60 s | 60.0 s | 0.017 Hz | 0.35% |
| Capella Giza | 32.9 s | 0.030 Hz | 0.63% |
| TerraSAR-X ST | ~5.5 s | 0.18 Hz | **2.9%** |
| ASNARO-2 SPOTLIGHT | ~1.3 s | 0.78 Hz | **12.2%** |
| Sentinel-1 IW | 0.141 s | 7.1 Hz | -- |
| *published campaigns* | -- | -- | *4.5-7.6%* |

ASNARO-2 overshoots the published aperture fraction, Capella undershoots it by
an order of magnitude, and TerraSAR-X ST lands beside it. Three sensors now
agree that **`t_sap` in seconds is the invariant and `alpha` is not.**

Run backwards, `t_sap` of 0.16-0.21 s at `alpha` 4.5-7.6% implies the published
campaigns worked with **effective dwells of about three seconds**, which is what
`MODIFIED-BACKPROJECTION.md` records the Trento experiment doing: displacement
reconstructed from the first three seconds of a roughly 20 s dwell, with
spectrogram clarity degrading after about 11 s. Whether a long collect should be
deliberately truncated rather than used whole is item 4 of `FOLLOW-UPS.md`.

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

### AFRL Gotcha, and the only open dataset that addresses a required test

The **Gotcha Volumetric SAR Data Set, Version 1.0** is X-band **phase history**
at 640 MHz bandwidth, with **full 360-degree azimuth coverage at eight different
elevation angles**, fully polarimetric, over a scene of civilian vehicles and
calibration targets. Airborne, circular collection geometry.

https://www.sdms.afrl.af.mil/index.php?collection=gotcha

**Why it matters here.** Test 6 of
[`IMPLEMENTATION-VERIFICATION.md`](IMPLEMENTATION-VERIFICATION.md) asks for an
independent-baseline reference:

> Process a conventional multi-baseline SAR tomography dataset with known
> elevation targets. This validates the steering/inversion machinery separately
> from the acoustic hypothesis.

Eight elevation passes over vehicles of known height is exactly that -- a real
elevation baseline with known three-dimensional structure, entirely separate
from the acoustic hypothesis the depth stage rests on. Every other source in
this document was assessed against the *micro-motion* stage; this is the only
open dataset found that speaks to the *tomography* stage, and the only one that
addresses any of the six required tests directly.

It also carries **calibration targets at known positions**. Not a moving
reflector, so it cannot touch test 1, but it would let geolocation and impulse
response be checked against ground truth for the first time -- which the failed
Rosamond search (see the corner-reflector section) left with no substitute.

**What it would cost.** Three things are unresolved and none is trivial:

- **Circular airborne geometry, not a satellite pass.** `rs_focus_backproject()`
  takes arbitrary platform positions and times, so the focusing is not the
  obstacle; the grid, the geodesy and every "along-track" assumption in the
  tomography stage are. `IMPLEMENTATION-VERIFICATION.md` records that the
  patent's baseline is along-track rather than elevation, and a circular collect
  is the case where those genuinely differ.
- **The scene is static.** It says nothing about micro-motion, and a depth
  product from it would be testing the inversion machinery alone. That is the
  point of test 6, but it must not be reported as anything else.
- **Access is by request** through the AFRL SDMS public site. Whether it is
  obtainable without US affiliation has **not** been verified, and should be
  before anything is planned around it.

## COSMO-SkyMed sample products

ESA publishes sample products suitable for inspecting the mission formats:

- First generation: https://earth.esa.int/eogateway/missions/cosmo-skymed/sample-data
- Second generation: https://earth.esa.int/eogateway/missions/cosmo-skymed-second-generation/sample-data
