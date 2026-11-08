#!/usr/bin/env python
"""Convert CRESCENT movie-tecplot output to Blender-ready heightmap sequences.

Reads the eta-only movie tecplot files of a nested run (zone 1: L0 full domain,
zone 2: L1 restricted to Omega_A) and writes, per frame, one 16-bit grayscale
PNG per level: the time-varying wave height in a form Blender's Displace
modifier reads directly as an image sequence. A meta.json records the common
normalization (a single symmetric eta range across all frames and both levels)
plus the grids and times, so the Blender scene can convert pixel values back to
meters exactly: eta = (pix/65535 - 0.5) * (2 * eta_abs).

Usage:
    python tools/movie/tecplot_to_heightmaps.py RUN_DIR [--prefix jonswap]
        [--out RUN_DIR/heightmaps]
"""

import argparse
import glob
import json
import os
import re

import numpy as np
from PIL import Image


def parse_movie(path):
    """Return [(nx, ny, x0, y0, dx, eta[ny,nx]), ...] for the zones of one file."""
    with open(path) as f:
        lines = f.readlines()
    zone_idx = [i for i, l in enumerate(lines) if l.startswith("ZONE")]
    zones = []
    for zi, ln in enumerate(zone_idx):
        m = re.search(r"I=(\d+), J=(\d+)", lines[ln])
        nx, ny = int(m.group(1)), int(m.group(2))
        stop = zone_idx[zi + 1] if zi + 1 < len(zone_idx) else len(lines)
        block = np.loadtxt(lines[ln + 1:stop])
        assert block.shape[0] == nx * ny, (path, nx * ny, block.shape)
        arr = block.reshape(ny, nx, 3)
        x = arr[0, :, 0]
        y = arr[:, 0, 1]
        zones.append(dict(nx=nx, ny=ny, x0=float(x[0]), y0=float(y[0]),
                          dx=float(x[1] - x[0]), dy=float(y[1] - y[0]),
                          eta=arr[:, :, 2]))
    return zones


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--prefix", default="jonswap")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.run_dir, "heightmaps")
    os.makedirs(os.path.join(out, "L0"), exist_ok=True)
    os.makedirs(os.path.join(out, "L1"), exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.run_dir,
                                          f"{args.prefix}_movie_tecplot*.dat")))
    assert files, "no movie tecplot files found"
    steps = [int(re.search(r"tecplot(\d+)", f).group(1)) for f in files]

    # Pass 1: global symmetric range so all frames share one normalization.
    eta_abs = 0.0
    parsed = []
    for f in files:
        zones = parse_movie(f)
        parsed.append(zones)
        for z in zones:
            eta_abs = max(eta_abs, float(np.abs(z["eta"]).max()))
    eta_abs *= 1.001  # keep the extremes strictly inside the pixel range

    # Pass 2: write 16-bit PNGs, frame numbers starting at 1 for Blender.
    for k, zones in enumerate(parsed):
        for lev, z in enumerate(zones):
            pix = np.clip((z["eta"] / (2 * eta_abs) + 0.5) * 65535.0,
                          0, 65535).astype(np.uint16)
            # image row 0 is the TOP in Blender's UV convention with flipped V;
            # keep y increasing upward by flipping rows here and using default UVs.
            img = Image.fromarray(pix[::-1, :], mode="I;16")
            img.save(os.path.join(out, f"L{lev}", f"h{k + 1:04d}.png"))

    meta = dict(
        eta_abs=eta_abs,
        n_frames=len(parsed),
        steps=steps,
        levels=[
            {**{k2: parsed[0][lev][k2]
                for k2 in ("nx", "ny", "x0", "y0", "dx", "dy")}}
            for lev in range(len(parsed[0]))
        ],
    )
    with open(os.path.join(out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(f"{len(parsed)} frames x {len(parsed[0])} levels -> {out}")
    print(f"eta_abs = {eta_abs:.4f} m ; L0 {meta['levels'][0]['nx']}^2 ; "
          f"L1 {meta['levels'][1]['nx']}^2 at ({meta['levels'][1]['x0']:.2f},"
          f"{meta['levels'][1]['y0']:.2f})")


if __name__ == "__main__":
    main()
