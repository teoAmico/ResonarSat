# Run outputs

Everything this software produces from a real collect lands here: rasters, depth
cubes and their sidecars, geocoded CSVs, and the manifest that says what produced
them. One directory per scene, and **one directory per run inside it**. Nothing
is written to a scene directory directly.

```
runs/
  README.md
  giza/
    YYYY-MM-DD-khufu-paper-phase/
      RUN.md                  # what was run, verbatim, and why
      giza_sar.png
      giza_tomogram.png
      giza_tomo.f32{,.hdr,.meta}
```

`tools/new-run.sh <scene> <suffix> "<question>"` creates the directory and seeds
the manifest.

## Why per-run directories

A depth product is not interpretable without knowing how it was made. The same
collect and the same target give different answers depending on the sub-aperture
route, the estimator, the look count, the coherence threshold and the assumed
constants — and two figures that look alike may have come from entirely
different processing.

Sidecars carry a `measurement_chain` line recording all of it, which covers the
machine-readable half. This directory rule covers the other half: a reader who
opens a figure finds, in the same folder, the commands that produced it and the
question it was meant to answer.

## Naming

`YYYY-MM-DD-<target>-<what-distinguishes-this-run>`

The suffix should name whatever a *different* run would change. `paper-phase`
and `pulse-correlation` are useful; `run2` and `final` are not, because they
tell a later reader nothing about why two directories differ.

## RUN.md

Every run directory carries one, containing:

- the exact command lines, copy-pasteable, including every option;
- the collect, with its size and dwell;
- the git commit the binary was built from;
- what question the run was meant to answer;
- and, once results exist, what it actually showed — including if that was
  nothing.

A run that produced a null result keeps its directory. Deleting the runs that
did not work is how a body of evidence quietly becomes a highlight reel.
