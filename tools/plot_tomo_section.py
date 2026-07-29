"""Render a depth section from a tomo --out cube with labelled axes and a legend.

`tomo --section-out` writes the section as bare pixels, one per window by one
per depth cell, which is unreadable at its native size and carries no scale. A
depth product without an axis invites being read as a picture of the ground, so
this draws the axes, the colour scale, and the assumed constants that set the
depth axis, straight from the .meta sidecar.

usage: plot_tomo_section.py CUBE.f32 OUT.png [ROW]

ROW is the range-window index the section is taken along; default is the middle.
Reads CUBE.f32.hdr for dimensions and CUBE.f32.meta for the geometry.
"""
import re
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def font(size):
    for p in ("/System/Library/Fonts/Menlo.ttc",
              "/System/Library/Fonts/Supplemental/Arial.ttf"):
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()


def meta(path, key, cast=str):
    """One value from the .meta sidecar, or None."""
    try:
        for line in open(path):
            m = re.match(rf"{re.escape(key)}\s+(.+)", line.strip())
            if m:
                return cast(m.group(1).strip())
    except Exception:
        pass
    return None


cube, out = sys.argv[1], sys.argv[2]
hdr = {k: v for k, v in
       (l.split(None, 1) for l in open(cube + ".hdr") if len(l.split(None, 1)) == 2)}
n_plane, n_row, n_col = (int(hdr["n_plane"]), int(hdr["n_row"]), int(hdr["n_col"]))
row = int(sys.argv[3]) if len(sys.argv) > 3 else n_row // 2

a = np.fromfile(cube, dtype=np.float32)[:n_plane * n_row * n_col]
a = a.reshape(n_plane, n_row, n_col)
sec = a[:, row, :].astype(float)          # [window_az, depth]

mp = cube + ".meta"
cell = meta(mp, "depth_cell_m", float) or 1.0
dmax = meta(mp, "depth_max_m", float) or n_col * cell
res = meta(mp, "depth_resolution_m", float)
unamb = meta(mp, "depth_unambiguous_m", float)
vel = meta(mp, "assumed_velocity_ms", float)
frq = meta(mp, "assumed_frequency_hz", float)
lam = meta(mp, "acoustic_wavelength_m", float)
conv = meta(mp, "wavelength_convention") or "?"
model = meta(mp, "model") or "?"

# Grid geometry, so the along-track axis can be metres rather than window index.
chain = meta(mp, "measurement_chain") or ""
win = re.search(r"win=(\d+)x", chain)
gcell = re.search(r"cell=([\d.]+)", chain)
win = int(win.group(1)) if win else 32
gcell = float(gcell.group(1)) if gcell else 1.0
stride = win // 2
extent = (n_plane - 1) * stride * gcell    # metres spanned by the section

# dB relative to the section's own peak, which is what a tomogram is read in.
mag = np.abs(sec)
peak = mag.max() if mag.max() > 0 else 1.0
db = 20 * np.log10(np.maximum(mag, 1e-30) / peak)
LO = -30.0
v = np.clip((db - LO) / (0 - LO), 0, 1)

# Jet, to match --palette jet, built explicitly so the legend cannot drift.
def jet(t):
    r = np.clip(1.5 - np.abs(4 * t - 3), 0, 1)
    g = np.clip(1.5 - np.abs(4 * t - 2), 0, 1)
    b = np.clip(1.5 - np.abs(4 * t - 1), 0, 1)
    return np.dstack([r, g, b])

# depth down the page, along-track across
img_arr = (jet(v.T) * 255).astype(np.uint8)      # [depth, window]
PW, PH = 900, 460
panel = Image.fromarray(img_arr).resize((PW, PH), Image.NEAREST)

L, T, R, B = 110, 96, 190, 92
canvas = Image.new("RGB", (L + PW + R, T + PH + B), (255, 255, 255))
canvas.paste(panel, (L, T))
d = ImageDraw.Draw(canvas)
f_ttl, f_ax, f_sm = font(21), font(15), font(12)

d.rectangle([L, T, L + PW, T + PH], outline=(0, 0, 0), width=1)

# ---- depth axis (vertical) ----
step = 5.0 if dmax <= 45 else 10.0
z = 0.0
while z <= dmax + 1e-9:
    y = T + (z / dmax) * PH
    d.line([(L - 6, y), (L, y)], fill=(0, 0, 0), width=1)
    d.text((L - 12 - d.textlength(f"{z:.0f}", font=f_ax), y - 8),
           f"{z:.0f}", fill=(0, 0, 0), font=f_ax)
    z += step
d.text((16, T + PH / 2 - 60), "d\ne\np\nt\nh\n\n(m)", fill=(0, 0, 0), font=f_ax)

# ---- along-track axis (horizontal) ----
for k in range(5):
    x = L + k * PW / 4
    xm = -extent / 2 + k * extent / 4
    d.line([(x, T + PH), (x, T + PH + 6)], fill=(0, 0, 0), width=1)
    lbl = f"{xm:+.0f}"
    d.text((x - d.textlength(lbl, font=f_ax) / 2, T + PH + 10), lbl,
           fill=(0, 0, 0), font=f_ax)
d.text((L + PW / 2 - 150, T + PH + 34),
       f"along-track position across the section (m), window row {row}",
       fill=(0, 0, 0), font=f_ax)

# ---- colour legend ----
CB = L + PW + 26
for i in range(PH):
    t = 1.0 - i / (PH - 1)
    c = tuple(int(x * 255) for x in jet(np.array([[t]]))[0, 0])
    d.line([(CB, T + i), (CB + 26, T + i)], fill=c)
d.rectangle([CB, T, CB + 26, T + PH], outline=(0, 0, 0), width=1)
for frac, lbl in ((0.0, "0"), (0.25, "-7.5"), (0.5, "-15"), (0.75, "-22.5"), (1.0, "-30")):
    y = T + frac * PH
    d.line([(CB + 26, y), (CB + 31, y)], fill=(0, 0, 0), width=1)
    d.text((CB + 35, y - 7), lbl, fill=(0, 0, 0), font=f_sm)
d.text((CB - 4, T - 22), "dB re peak", fill=(0, 0, 0), font=f_sm)

# ---- title and the constants that set the axis ----
d.text((L, 16), "Depth section - patent-exact chain, Great Pyramid (Khufu)",
       fill=(0, 0, 0), font=f_ttl)
sub = (f"model {model}   lambda = {lam:.3g} m ({conv} convention)"
       f"   ASSUMED v = {vel:.0f} m/s, f = {frq:.0f} Hz")
d.text((L, 44), sub, fill=(40, 40, 40), font=f_ax)
line3 = f"depth cell {cell:.2f} m"
if res:
    line3 += f"   resolution {res:.2f} m"
if unamb:
    line3 += f"   unambiguous to {unamb:.1f} m"
d.text((L, 66), line3, fill=(40, 40, 40), font=f_ax)

# ---- the caveat, on the figure rather than only in the sidecar ----
d.text((L, T + PH + 58),
       "NOT A MEASUREMENT OF DEPTH. The tracked observable feeding this "
       "inversion was zero in 220 of 225 windows;",
       fill=(150, 0, 0), font=f_sm)
d.text((L, T + PH + 73),
       "the depth axis is set by the assumed v and f above, neither of which "
       "this pipeline measures. See the .meta sidecar.",
       fill=(150, 0, 0), font=f_sm)

canvas.save(out)
print(f"wrote {out}  ({canvas.width}x{canvas.height})")
print(f"  section: {n_plane} windows x {n_col} depth cells, row {row}")
print(f"  span {extent:.0f} m along-track, 0 to {dmax:.1f} m depth")
print(f"  peak |h| {peak:.4g}, dynamic range shown {LO:.0f} dB")
