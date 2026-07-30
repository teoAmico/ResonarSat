#!/usr/bin/env python3
"""Summarise a --shifts CSV from resonarsat mmotion.

Checks the phase-fold tell documented in include/resonarsat/microm.h: a
peak-to-peak line-of-sight velocity near lambda/(2*dt) means the unwrapped
series wrapped rather than that the target moved.
"""
import csv
import sys
import statistics as st

LAMBDA = 0.0322  # m, from the collect's carrier

path = sys.argv[1]
dt = None
with open(path) as f:
    lines = f.readlines()

for ln in lines[:4]:
    if ln.startswith('#') and 'dt_s=' in ln:
        for tok in ln.replace('#', '').split():
            if tok.startswith('dt_s='):
                dt = float(tok.split('=')[1])

hdr_i = next(i for i, ln in enumerate(lines) if ln.startswith('window,'))
rows = list(csv.DictReader(lines[hdr_i:]))

wins = {}
for r in rows:
    wins.setdefault(int(r['window']), []).append(r)

vel_pp, disp_pp, qual, phase_step = [], [], [], []
for w, rs in sorted(wins.items()):
    v = [float(r['vel_los_ms']) for r in rs]
    d = [float(r['disp_los_m']) for r in rs]
    p = [float(r['phase_rad']) for r in rs]
    vel_pp.append(max(v) - min(v))
    disp_pp.append(max(d) - min(d))
    qual.append(float(rs[0]['quality']))
    phase_step.append(max(abs(p[i + 1] - p[i]) for i in range(len(p) - 1)))

fold = LAMBDA / (2 * dt) if dt else float('nan')
n = len(wins)
print(f"file            {path}")
print(f"windows         {n}   looks {len(rows)//n}   dt {dt:.6g} s")
print(f"fold ceiling    lambda/(2*dt) = {fold*1e3:.1f} mm/s")
print()
print(f"peak-to-peak LOS velocity   median {st.median(vel_pp)*1e3:9.1f} mm/s"
      f"   max {max(vel_pp)*1e3:9.1f} mm/s")
print(f"  as a fraction of the fold ceiling: median {st.median(vel_pp)/fold:.3f}"
      f"   max {max(vel_pp)/fold:.3f}")
print(f"  windows at or above the ceiling:   {sum(1 for x in vel_pp if x >= fold)} of {n}")
print(f"  windows above half the ceiling:    {sum(1 for x in vel_pp if x >= 0.5*fold)} of {n}")
print()
print(f"peak-to-peak LOS displacement  median {st.median(disp_pp)*1e3:9.2f} mm"
      f"   max {max(disp_pp)*1e3:9.2f} mm")
print(f"  lambda/2 (one phase cycle) = {LAMBDA/2*1e3:.1f} mm;"
      f" windows exceeding it: {sum(1 for x in disp_pp if x > LAMBDA/2)} of {n}")
print()
print(f"max |phase step| between looks  median {st.median(phase_step):.3f} rad"
      f"   max {max(phase_step):.3f} rad   (pi = 3.142; above it the unwrap is ambiguous)")
print(f"  windows whose largest step exceeds pi: "
      f"{sum(1 for x in phase_step if x > 3.14159)} of {n}")
print()
print(f"tracking quality   median {st.median(qual):.3f}   min {min(qual):.3f}"
      f"   max {max(qual):.3f}")
