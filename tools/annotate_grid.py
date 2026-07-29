"""Overlay a labelled offset grid, and optional targets, on a focus --raw dump.

usage: annotate_grid.py RAW.f32 N CELL OX OY OUT.png [STEP] [TARGETS]
       TARGETS = "name:X:Y,name:X:Y,..."   metres, same frame as --offset

Axes follow focus.c:117-121 -- image ROW is the grid's X (n_x/uIAX) and image
COLUMN is its Y (n_y/uIAY), with the grid centre at --offset. Every gridline is
therefore labelled in the metres that --offset and --at speak, so a feature can
be pointed at and converted back without guessing.

NOTE ON WHAT "X DOWN, Y RIGHT" MEANS HERE. On the Giza collect the product's
image-area axes are LEFT-handed with respect to local up --
(uIAX x uIAY) . up = -0.9999958 -- so this layout is mirrored relative to a
north-up map. Compass intuition about which end of a line is north will be
wrong; read positions off the grid rather than inferring them from orientation.
"""
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

def _font(size):
    for p in ("/System/Library/Fonts/Menlo.ttc",
              "/System/Library/Fonts/Supplemental/Arial.ttf",
              "/Library/Fonts/Arial.ttf"):
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()

def label(d, xy, text, fg, font):
    """Text on a dark plate, so it stays readable over bright speckle."""
    x, y = xy
    try:
        l, t, r, b = d.textbbox((x, y), text, font=font)
    except Exception:
        r, b = x + 8 * len(text), y + 14
        l, t = x, y
    d.rectangle([l - 3, t - 2, r + 3, b + 2], fill=(0, 0, 0))
    d.text((x, y), text, fill=fg, font=font)

path, n, cell, ox, oy, out = (sys.argv[1], int(sys.argv[2]), float(sys.argv[3]),
                              float(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
step = float(sys.argv[7]) if len(sys.argv) > 7 else 100.0

a = np.fromfile(path, dtype=np.float32)[: n * n].reshape(n, n).astype(float)

# Percentile stretch in dB. Referencing the peak wastes the whole ramp when a
# single corner reflector is 40 dB above the scene, which is the usual case
# here and is why the first attempt came out nearly black.
db = 20 * np.log10(np.maximum(a, 1e-30))
lo, hi = np.percentile(db, 5), np.percentile(db, 99.5)
v = np.clip((db - lo) / (hi - lo), 0, 1)
img = Image.fromarray((v * 255).astype(np.uint8)).convert("RGB")

# Upscale so labels are legible at 1:1.
S = max(1, int(round(1200 / n)))
img = img.resize((n * S, n * S), Image.NEAREST)
d = ImageDraw.Draw(img)

def px_of(X, Y):
    """metres -> pixel, inverting focus.c's cell-centre mapping."""
    r = (X - ox) / cell + (n - 1) / 2.0
    c = (Y - oy) / cell + (n - 1) / 2.0
    return c * S, r * S      # PIL is (x=col, y=row)

half = n * cell / 2.0
GRID = (0, 200, 255)
def rng(o):
    lo_ = np.ceil((o - half) / step) * step
    return np.arange(lo_, o + half + 1e-9, step)

fnt = _font(max(13, int(img.width / 70)))
for X in rng(ox):                      # lines of constant X run across columns
    _, y = px_of(X, oy)
    d.line([(0, y), (img.width, y)], fill=GRID, width=1)
    label(d, (5, y + 3), f"X={X:+.0f}", GRID, fnt)
    label(d, (img.width - 105, y + 3), f"X={X:+.0f}", GRID, fnt)
for Y in rng(oy):
    x, _ = px_of(ox, Y)
    d.line([(x, 0), (x, img.height)], fill=GRID, width=1)
    label(d, (x + 4, 5), f"Y={Y:+.0f}", GRID, fnt)
    label(d, (x + 4, img.height - 24), f"Y={Y:+.0f}", GRID, fnt)

TARGETS = {}
if len(sys.argv) > 8 and sys.argv[8].strip():
    for item in sys.argv[8].split(","):
        nm, sx, sy = item.split(":")
        TARGETS[nm] = (float(sx), float(sy))
for name, (X, Y) in TARGETS.items():
    if abs(X - ox) > half or abs(Y - oy) > half:
        continue
    x, y = px_of(X, Y)
    r = 9 * S
    col = (255, 60, 60)
    d.ellipse([x - r, y - r, x + r, y + r], outline=col, width=3)
    d.line([(x - 1.6 * r, y), (x - 0.4 * r, y)], fill=col, width=2)
    d.line([(x + 0.4 * r, y), (x + 1.6 * r, y)], fill=col, width=2)
    d.line([(x, y - 1.6 * r), (x, y - 0.4 * r)], fill=col, width=2)
    d.line([(x, y + 0.4 * r), (x, y + 1.6 * r)], fill=col, width=2)
    label(d, (x + r + 6, y - 10), f"{name}  ({X:+.0f},{Y:+.0f})", col, fnt)

label(d, (6, img.height - 52),
      f"grid {step:.0f} m | centre --offset {ox:+.0f},{oy:+.0f} | cell {cell} m | "
      f"X down (uIAX), Y right (uIAY)", (255, 255, 0), fnt)
img.save(out)
print(f"wrote {out}  ({img.width}x{img.height})")
