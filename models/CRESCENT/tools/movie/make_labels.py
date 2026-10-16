#!/usr/bin/env python
"""Render the mathtext label PNGs used by movie_scene.py.

Each label's font color matches the boundary it names: white for the L0
domain frame, the frame orange for the L1 window. Transparent background,
serif mathtext, high dpi so the plane texture stays crisp at 4K.

Usage (from the repo root):
    python tools/movie/make_labels.py [--out movie]
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

# (file, text, color, stroke) -- L1 uses the frame-bar orange of
# movie_scene.py; the stroked variant stays legible over sun glitter and is
# used by the plan-view camera
LABELS = [
    ("label_L0.png", r"$L_0$ domain", "white", False),
    ("label_L1.png", r"$L_1$ domain", (1.0, 0.42, 0.05), False),
    ("label_L1_stroke.png", r"$L_1$ domain", (1.0, 0.42, 0.05), True),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="movie")
    args = ap.parse_args()
    for fname, text, color, stroke in LABELS:
        fig = plt.figure(figsize=(4, 1), dpi=400)
        fig.patch.set_alpha(0)
        kw = dict(path_effects=[pe.withStroke(
            linewidth=4, foreground=(0.02, 0.05, 0.08))]) if stroke else {}
        fig.text(0.5, 0.5, text, ha="center", va="center",
                 fontsize=30, color=color, family="serif", **kw)
        path = os.path.join(args.out, fname)
        fig.savefig(path, transparent=True, bbox_inches="tight",
                    pad_inches=0.05)
        plt.close(fig)
        print("wrote", path)


if __name__ == "__main__":
    main()
