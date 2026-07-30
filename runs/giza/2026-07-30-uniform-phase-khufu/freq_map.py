#!/usr/bin/env python3
"""Per-window dominant frequency map from a --shifts CSV.

Mirrors rs_spectrum_compute_opts(): linear detrend, periodogram of the
line-of-sight DISPLACEMENT (the phase estimator's observable), dominant bin
above an optional floor. Exists to answer a question the CLI does not print:
whether the winning frequency occupies a contiguous patch of windows, or a
different window each time the search is nudged.
"""
import csv
import sys
import numpy as np

path = sys.argv[1]
fmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

with open(path) as f:
    lines = f.readlines()
dt = None
for ln in lines[:4]:
    if 'dt_s=' in ln:
        for tok in ln.replace('#', '').split():
            if tok.startswith('dt_s='):
                dt = float(tok.split('=')[1])
hdr = next(i for i, ln in enumerate(lines) if ln.startswith('window,'))
rows = list(csv.DictReader(lines[hdr:]))

wins = {}
for r in rows:
    wins.setdefault(int(r['window']), []).append(float(r['disp_los_m']))
n_win = len(wins)
side = int(round(n_win ** 0.5))
n_looks = len(wins[0])

freqs = np.fft.rfftfreq(n_looks, dt)
lo = np.searchsorted(freqs, max(fmin, 1e-12))

dom = np.zeros(n_win)
prom = np.zeros(n_win)
for w in sorted(wins):
    y = np.array(wins[w])
    t = np.arange(n_looks)
    y = y - np.polyval(np.polyfit(t, y, 1), t)          # linear detrend
    p = np.abs(np.fft.rfft(y)) ** 2
    p[0] = 0.0
    band = p[lo:]
    k = int(np.argmax(band)) + lo
    dom[w] = freqs[k]
    prom[w] = p[k] / np.mean(band) if np.mean(band) > 0 else 0.0

g = dom.reshape(side, side)
pg = prom.reshape(side, side)
print(f"{path}   fmin {fmin} Hz   df {freqs[1]:.4f} Hz   f_max {freqs[-1]:.2f} Hz")
print(f"\ndominant frequency per window (Hz), {side}x{side} grid:")
for r in g:
    print("  " + " ".join(f"{v:6.2f}" for v in r))
print(f"\nprominence per window (peak / mean of band):")
for r in pg:
    print("  " + " ".join(f"{v:6.1f}" for v in r))

order = np.argsort(prom)[::-1][:8]
print("\ntop 8 windows by prominence:")
for w in order:
    print(f"  window {w:4d} at grid ({w//side},{w%side})   "
          f"{dom[w]:7.3f} Hz   prominence {prom[w]:6.1f}")

# Do the top windows touch each other?
top = set(int(w) for w in order[:8])
adj = 0
for w in top:
    r, c = w // side, w % side
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        if (r + dr) * side + (c + dc) in top and 0 <= r + dr < side and 0 <= c + dc < side:
            adj += 1
print(f"\nadjacent pairs among the top 8 windows: {adj // 2}")
print(f"spread of their dominant frequencies: "
      f"{dom[order].min():.3f} to {dom[order].max():.3f} Hz")
