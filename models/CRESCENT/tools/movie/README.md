# Making ocean movies from CRESCENT runs with Blender

This folder turns a nested CRESCENT run into a rendered movie: the L0 parent
sea extended periodically to the horizon, the L1 child's own solution inset in
its window, a photographic sky, and formal $L_0$/$L_1$ annotations. Everything
runs headless from the command line; no admin rights and no Blender MCP are
required.

Produced with this pipeline: `movie/crescent_broadband_L0_L1_*.mp4` (broadband
directional JONSWAP sea, 24 peak periods, low grazing camera) and
`movie/crescent_wave_*.mp4` (the class-II crescent instability, 6 carrier
periods, near-real-time, medium vista with sky via `--static-cam
--cam-height 0.44 --cam-back 0.95 --aim-height 0.14`, surface = the literal
simulation field; the `_lowview` and `_planview` files keep the earlier
grazing and 55-degree versions). The `*_clean_*.mp4` files are
annotation-free variants (waves and sky only): the sea pass rendered with
`--l0-only` and encoded directly -- no annotation pass, no compositing.
`--l0-only` shows the continuous parent field with no inset child patch: a
data splice always leaves a faint rim seam (child-parent mismatch plus the
z-fighting lift), and without the frame bars covering it even a
`--uniform-tint` splice stays traceable by eye. The
`crescent_broadband_clean_static_*.mp4` pair holds the camera at the
push-in's first frame (`--static-cam` with `--l0-only`).

## Prerequisites (all user-space)

- **Blender** 4.x/5.x. Copy the app out of the DMG so renders survive
  unmounts: `cp -R /Volumes/Blender/Blender.app ~/Applications/`. The binary
  is `~/Applications/Blender.app/Contents/MacOS/Blender`.
- **Python venv** with numpy, Pillow, matplotlib (`~/venvs/qiskit313` here).
- **Encoder**: `python -m pip install imageio-ffmpeg` ships a static ffmpeg
  binary (`python -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())"`).
- **Sky image**: any wide photo; crop the sky portion (here
  `movie/cloud_upper.jpg`, the upper half of the user's cloud photograph).
- **Labels**: `python tools/movie/make_labels.py` renders the mathtext label
  PNGs ("$L_0$ domain" etc.) on transparent backgrounds. Each label's font
  color matches the boundary it names (white L0 frame, orange L1 window), and
  a dark-stroked L1 variant stays legible over sun glitter in plan views.
  Labels sit *inside* the window they name, never touching its boundary.

## Pipeline overview

```
solver movie-tecplot output      (eta only; zone 1 = L0 full, zone 2 = L1 core)
        |  tecplot_to_heightmaps.py
16-bit PNG heightmaps + meta.json   (the Blender-readable wave height)
        |  movie_scene.py  --render-pass sea          (animation)
        |  movie_scene.py  --render-pass annotations  (animation, transparent)
PNG frame sequences
        |  PIL alpha_composite (per frame)
composited frames
        |  ffmpeg (imageio-ffmpeg static binary)
mp4 (4K / 2K / ...)
```

## Step 1 — produce movie-cadence run data

The solver writes a compact movie tecplot (`<prefix>_movie_tecplot*.dat`,
eta only, L0 full domain + L1 restricted to its accurate core) at every
`crescent.tecplot_int` step. For smooth motion make the frame interval a small
fraction of the peak/carrier period: aim for **T/10 or finer per frame**
(JONSWAP: `tecplot_int = 5` at dt = Tp/50; crescent: `tecplot_int = 2` at
dt = T/100, i.e. T/50 per frame).

The full tecplot files written at the same cadence are large; run a deleter
loop alongside that removes `<prefix>_tecplot*.dat` once the next movie file
exists (see the session history for the exact loop), keeping only the ~MB
movie files. Never reuse an existing run directory.

## Step 2 — convert to heightmaps

```bash
python tools/movie/tecplot_to_heightmaps.py RUN_DIR --prefix <plot_file>
```

Writes `RUN_DIR/heightmaps/{L0,L1}/h0001..NNNN.png` (16-bit grayscale, one
global symmetric normalization across all frames and both levels) and
`meta.json` holding `eta_abs`, frame count, and each level's grid
(`nx, ny, x0, y0, dx, dy` — rectangular cells supported). Exact calibration:
`eta = (pix/65535 - 0.5) * 2 * eta_abs`. 16 bits give ~10 um resolution on a
+-0.33 m field; Blender reads the sequence directly as displacement.

## Step 3 — build and render the scene

`movie_scene.py` builds everything (run once per pass):

```bash
B=~/Applications/Blender.app/Contents/MacOS/Blender
M=RUN_DIR/heightmaps
# sea pass (animation)
$B --background --python tools/movie/movie_scene.py -- \
   --maps $M --sky movie/cloud_upper.jpg --out OUT --anim --render-pass sea \
   --res 3840 2160 [--exag 5.5] [--static-cam] [--save-blend scene.blend]
# annotation pass (transparent overlay)
$B --background --python tools/movie/movie_scene.py -- \
   --maps $M --sky movie/cloud_upper.jpg --out OUT_ANN --anim \
   --render-pass annotations --res 3840 2160 [--exag ...] [--static-cam]
```

With `--static-cam` the annotation pass can be a single `--still 1` instead of
an animation (the overlay never changes).

Scene contents, all derived from `meta.json` so any domain size works:

- **Hi-res sea**: 3x3 tile-aligned periodic images of the L0 field (the HOS
  parent is doubly periodic, so tiling is its exact continuation), ~2300^2
  vertices. Keeps full mesh detail well past the marked domain.
- **Far field**: 9x9 tiles on a lighter mesh with a hole cut under the hi-res
  sea; carries the ocean to the horizon.
- **L1 patch**: the child's own core solution, inset into a hole cut under the
  child window, lifted by 0.006 x (displayed amplitude) + 0.05 x (child cell)
  to avoid z-fighting where it overlaps the hole rim.
- **Annotations**: emissive frame bars (orange = $L_1$ window, white = $L_0$
  cell) and mathtext label planes, placed relative to the domain extent
  horizontally and to the displayed amplitude vertically.
- **Sky**: emissive backdrop plane (the photo the camera sees) + the same
  image as world environment (reflections/ambient).
- **Camera**: slow push-in keyed from a wide vista to low over the child
  window (Bezier ease). Make camera moves large: a 1.2x apparent-size change
  over 10 s is imperceptible; ~1.8x reads clearly. `--static-cam` holds the
  start framing. `--cam-height/--cam-back/--aim-height` raise and pull back
  the static camera while keeping the horizon in frame (raise the aim with
  the camera) for a medium vista. `--cam-elev DEG --cam-dist F` instead
  places a static high-oblique camera (elevation angle above the horizon
  seen from the domain center, slant distance in units of the L0 extent)
  for a plan-style view -- the classic way to show the staggered crescent
  pattern, at the price of the sky leaving the frame.
- **High-view mode** (any `--cam-elev` run) switches the scene dressing: the
  grazing look lives off dark water and specular sky, which reads flat from
  above, so high views use a brighter base, no coat, a raking northwest sun
  (shades the wave fronts and forward-scatters glitter to the camera) and
  less world ambient. `--ripple` additionally adds a drifting geometric
  capillary ripple (7 cm, ~1.5% of the displayed amplitude). It is cosmetic
  sub-grid texture -- OFF by default so the rendered surface is the literal
  simulation field; without it plan views look soft, which is the honest
  spectral truncation of HOS. Disclose it in any caption if used.
- **Exaggeration**: `--exag` scales eta for visibility. Gently sloped seas
  need it (JONSWAP ak~0.06: 5.5); steep seas do not (crescent ak=0.33: use
  ~1, the true geometry is already dramatic). The true amplitude stays in
  `meta.json`.

Regenerate the labels if the wording or colors change:

```bash
python tools/movie/make_labels.py
```

## Step 4 — composite and encode

Composite into a SEPARATE directory -- overwriting the sea frames in place
bakes the overlay in, so any later annotation-only tweak forces a full sea
re-render instead of a cheap re-composite:

```bash
python - <<'EOF'
import glob, os
from PIL import Image
os.makedirs("OUT/composited", exist_ok=True)
sea = sorted(glob.glob("OUT/frames/f*.png"))
ann = sorted(glob.glob("OUT_ANN/frames/f*.png")) or \
      ["OUT_ANN/still_0001.png"] * len(sea)       # static cam: one still
for fs, fa in zip(sea, ann):
    Image.alpha_composite(Image.open(fs).convert("RGBA"),
                          Image.open(fa).convert("RGBA")).convert("RGB") \
         .save("OUT/composited/" + os.path.basename(fs))
EOF
FF=$(python -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())")
$FF -y -framerate 24 -i OUT/composited/f%04d.png -c:v libx264 -pix_fmt yuv420p \
    -crf 17 -preset medium -movflags +faststart out_4k.mp4
$FF -y -framerate 24 -i OUT/composited/f%04d.png -vf scale=2560:1440 ... out_2k.mp4
```

## Hard-won constraints (why the pipeline looks like this)

| symptom | cause | fix baked into the scripts |
|---|---|---|
| every rendered frame identical (movie ~0.1 MB) | image *sequences* on displace textures do not refresh in headless renders | single image datablock + `frame_change_pre` handler that swaps `img.filepath` and reloads per frame |
| dark + bright wavy lines under the marker bars | the bars cast sun shadows and specular reflections on the water; EEVEE honors per-object shadow visibility but ignores *reflection* visibility | two-pass render: sea without annotations, annotations on transparent film, composited (Cycles honors the flags but is ~60x slower on CPU) |
| `enum "FFMPEG" not found` | Blender 5.x removed the built-in movie encoder | render PNG frames, encode with imageio-ffmpeg's static binary |
| near-field "not smooth" outside the domain frame | that strip belonged to the coarse far-field mesh (its vertex spacing low-passes the same height data); the density seam also leaves silhouette specks | hi-res sea spans 3x3 periodic tiles, pushing the seam ~one domain-length away, off-frame |
| sky/clouds invisible, gray gradient instead | backdrop plane beyond the camera's default 1000 m clip | `clip_end` raised (scaled with the domain) |
| a label never visible | placed outside the camera frustum at that depth (perspective misjudged) | verify placement with an annotation-only still; label positions are domain-relative |
| push-in "not visible" | 1.2x apparent-size change over 10 s is below perception against wave motion | measure the effect on encoded frames (track the orange box width); use ~1.8x |
| render job killed mid-chain | long chained background jobs can hit time limits | run sea pass, annotation pass, composite/encode as separate jobs |
| thick slab-like frame bars on a small domain | bar cross-section scaled with wave amplitude, which is large *relative to* a small steep domain | bars are thin squares tied to the domain extent (`t x t`), only their height above water scales with amplitude |
| faint dark outline on the water under the child window | the z-fighting lift of the L1 patch casts a silhouette/shadow line, visible when the lift is large relative to the domain | lift = 0.006 x amplitude + 0.05 x child cell (amplitude piece for the child-parent mismatch, cell piece for interpolation error) |
| `--cycles` flag rejected | Blender parses ambiguous CLI prefixes even after `--` | name custom flags to avoid Blender's (`--use-cycles`) |
| mp4 unplayable in some players | odd chroma/pixel formats | always `-pix_fmt yuv420p`, `-movflags +faststart` |
| giant dark swath across the sea after moving the sun | the emissive sky backdrop plane casts a shadow once the sun sits behind it | `visible_shadow = False` on the backdrop |
| water flat and out-of-focus from a high camera | almost no Fresnel reflection from above, ambient fill compresses shading, and the HOS surface is glassy below the grid scale; a shader Bump node had no visible effect in headless EEVEE Next | high-view mode: raking sun, less ambient, brighter base, and a *geometric* drifting ripple as a second Displace modifier |
| dotted glint line along the child-window rim in plan views | the 3D ripple noise is sampled *after* the wave displacement, so patch and sea decorrelate across the lift and the sea pokes through; open mesh edges also have biased normals that sparkle | z-crushed ripple sampling space (2D noise), a downward skirt on the patch boundary, and 3x wider frame bars in plan views to blanket the seam |
| a label placed inside a small window appears outside it (or slides across it) | the label floats above the water, so its image parallax-shifts relative to the ground; the frame bars are elevated too and shift by their own amount | place labels at the BAR height (co-planar with the frame they sit in): both shift together, so "inside the frame" holds from any camera |

## Performance notes (Apple Silicon, EEVEE)

720p ~1 s/frame; 4K ~6 s/frame with ~8M displaced vertices. A 241-frame 4K
two-pass production is ~40 min end to end. Cycles gives full ray-visibility
control and richer water but ran ~60 s/frame on CPU here; use it only for
single-still checks unless GPU (Metal) rendering is verified.
