# Blender headless scene: CRESCENT nested sea, L0 parent + L1 child.
#
# Builds an ocean from the heightmap sequences written by
# tecplot_to_heightmaps.py. The high-resolution sea spans 3x3 tile-aligned
# periodic images of the L0 field (its exact continuation), a coarser 9x9-tile
# far field carries it to the horizon, and the L1 patch (the child's own
# solution over its accurate core) is inset into a hole cut under the child
# window. The upper half of the user's cloud photograph forms the sky
# (emissive backdrop plane plus world environment for reflections).
#
# The layout is domain-relative: horizontal positions scale with the L0 extent
# and vertical annotation placement scales with the displayed wave amplitude
# (eta_abs * exag), so the same scene works for the 256 m JONSWAP sea and the
# 12.6 m crescent box. Wave heights are exaggerated by --exag (default 5.5 for
# the gently sloped JONSWAP sea; use ~1 for the steep crescent case, whose
# true slopes need no help); the true amplitude is recorded in meta.json.
#
# Learned constraints baked in (do not "simplify" these away):
# - Image *sequences* on displace textures do not refresh in headless renders;
#   a frame_change_pre handler swaps and reloads a single image per frame.
# - EEVEE ignores per-object *reflection* visibility (Cycles honors it), so
#   annotations are rendered as a separate transparent pass and composited.
# - Blender 5.x removed the FFMPEG movie format; animations render PNG frames
#   for external encoding (imageio-ffmpeg works without admin rights).
#
# Run (each pass; see tools/movie/README.md for the full pipeline):
#   Blender --background --python tools/movie/movie_scene.py -- \
#       --maps RUN_DIR/heightmaps --sky movie/cloud_upper.jpg --out OUT_DIR \
#       [--still N | --anim] [--render-pass sea|annotations|all] \
#       [--exag 5.5] [--static-cam | --cam-elev 55 --cam-dist 1.3] \
#       [--res 3840 2160] [--save-blend f.blend]

import argparse
import json
import math
import os
import sys

import bpy

ANNOTATIONS = []    # objects excluded from the sea pass (bars, labels)
HEIGHT_IMAGES = []  # (image datablock, directory) advanced per frame


def parse_args():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    ap = argparse.ArgumentParser()
    ap.add_argument("--maps", required=True)
    ap.add_argument("--sky", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--still", type=int, default=None)
    ap.add_argument("--anim", action="store_true")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--res", type=int, nargs=2, default=[1280, 720])
    ap.add_argument("--exag", type=float, default=5.5,
                    help="vertical exaggeration of eta (use ~1 for steep seas)")
    ap.add_argument("--static-cam", action="store_true",
                    help="hold the camera at the push-in start position")
    ap.add_argument("--cam-elev", type=float, default=None,
                    help="static high-oblique camera: elevation angle (deg) "
                         "above the horizon seen from the domain center; "
                         "overrides the low push-in framing")
    ap.add_argument("--cam-dist", type=float, default=1.4,
                    help="camera slant distance from the domain center in "
                         "units of the L0 x-extent (with --cam-elev)")
    ap.add_argument("--cam-height", type=float, default=0.203,
                    help="camera height as a fraction of the L0 x-extent "
                         "(low/medium framings; default is the low vista)")
    ap.add_argument("--cam-back", type=float, default=0.586,
                    help="camera distance behind the near L0 edge as a "
                         "fraction of the L0 y-extent")
    ap.add_argument("--aim-height", type=float, default=0.0,
                    help="look-at height fraction of the L0 x-extent; raise "
                         "with --cam-height to keep the horizon in frame")
    ap.add_argument("--l0-only", action="store_true",
                    help="do not inset the L1 patch: show the continuous "
                         "parent field only (a data splice always leaves a "
                         "faint seam; use this for annotation-free movies)")
    ap.add_argument("--uniform-tint", action="store_true",
                    help="give the child patch the parent water color "
                         "(for annotation-free renders, where the subtle "
                         "tint cue would read as a ghost rectangle)")
    ap.add_argument("--ripple", action="store_true",
                    help="add the cosmetic sub-grid capillary ripple (plan "
                         "views only; OFF by default so the surface is the "
                         "literal simulation field)")
    ap.add_argument("--save-blend", default=None)
    ap.add_argument("--use-cycles", action="store_true")
    ap.add_argument("--render-pass", choices=["all", "sea", "annotations"],
                    default="all")
    return ap.parse_args(argv)


def clean_scene():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for block in (bpy.data.meshes, bpy.data.materials, bpy.data.textures,
                  bpy.data.images, bpy.data.curves):
        for item in list(block):
            if item.users == 0:
                block.remove(item)


def seq_texture(name, first_png, repeat_x=1, repeat_y=1):
    img = bpy.data.images.load(first_png)
    img.colorspace_settings.name = "Non-Color"
    HEIGHT_IMAGES.append((img, os.path.dirname(first_png)))
    tex = bpy.data.textures.new(name, type="IMAGE")
    tex.image = img
    tex.extension = "REPEAT" if max(repeat_x, repeat_y) > 1 else "EXTEND"
    tex.repeat_x, tex.repeat_y = repeat_x, repeat_y
    tex.use_interpolation = True
    return tex


def advance_heightmaps(scene, _depsgraph=None):
    f = max(1, scene.frame_current)
    for img, d in HEIGHT_IMAGES:
        path = os.path.join(d, f"h{f:04d}.png")
        if os.path.exists(path) and img.filepath != path:
            img.filepath = path
            img.reload()


def water_material(name, tint, rough=0.22, coat=0.4):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*tint, 1.0)
    bsdf.inputs["Roughness"].default_value = rough
    if "Coat Weight" in bsdf.inputs:
        bsdf.inputs["Coat Weight"].default_value = coat
    return mat


def make_grid(name, nvx, nvy, center, size_xy, z, material):
    bpy.ops.mesh.primitive_grid_add(x_subdivisions=nvx, y_subdivisions=nvy,
                                    size=1.0, location=(center[0], center[1], z))
    ob = bpy.context.active_object
    ob.name = name
    ob.scale = (size_xy[0], size_xy[1], 1.0)
    bpy.ops.object.transform_apply(scale=True)
    ob.data.materials.append(material)
    bpy.ops.object.shade_smooth()
    return ob


def add_displace(ob, tex, strength):
    mod = ob.modifiers.new("eta", "DISPLACE")
    mod.texture = tex
    mod.texture_coords = "UV"
    mod.direction = "Z"
    mod.mid_level = 0.5
    mod.strength = strength


def cut_hole(ob, xlo, xhi, ylo, yhi):
    """Delete faces of ob whose centers lie inside the rectangle (world xy)."""
    import bmesh
    bm = bmesh.new()
    bm.from_mesh(ob.data)
    mw = ob.matrix_world
    doomed = [f for f in bm.faces
              if (lambda c: xlo < c.x < xhi and ylo < c.y < yhi)
              (mw @ f.calc_center_median())]
    bmesh.ops.delete(bm, geom=doomed, context="FACES")
    bm.to_mesh(ob.data)
    bm.free()


def add_skirt(ob, depth):
    """Extrude the open boundary downward: open-edge vertex normals average
    one-sided faces and sparkle under sun + ripple; a skirt closes the rim
    (and hides the parent's hole-cut edge in a shadowed slot)."""
    import bmesh
    bm = bmesh.new()
    bm.from_mesh(ob.data)
    edges = [e for e in bm.edges if len(e.link_faces) == 1]
    ret = bmesh.ops.extrude_edge_only(bm, edges=edges)
    for v in [g for g in ret["geom"] if isinstance(g, bmesh.types.BMVert)]:
        v.co.z -= depth
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(ob.data)
    bm.free()


def annotation_only(ob):
    ANNOTATIONS.append(ob)
    for attr in ("visible_shadow", "visible_glossy", "visible_diffuse",
                 "visible_transmission", "visible_volume_scatter"):
        if hasattr(ob, attr):
            setattr(ob, attr, False)


def emissive(name, color, strength):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    nodes.clear()
    em = nodes.new("ShaderNodeEmission")
    em.inputs["Color"].default_value = (*color, 1.0)
    em.inputs["Strength"].default_value = strength
    out = nodes.new("ShaderNodeOutputMaterial")
    mat.node_tree.links.new(em.outputs[0], out.inputs["Surface"])
    return mat


def main():
    a = parse_args()
    meta = json.load(open(os.path.join(a.maps, "meta.json")))
    n = meta["n_frames"]
    eta_abs = meta["eta_abs"]
    L0, L1 = meta["levels"]
    for lev in (L0, L1):
        lev.setdefault("dy", lev["dx"])
    strength = 2.0 * eta_abs * a.exag   # displacement over the full pixel range
    V = eta_abs * a.exag                # displayed wave half-amplitude

    # domain extents and the horizontal layout scale (design was the 256 m sea)
    ex0, ey0 = L0["nx"] * L0["dx"], L0["ny"] * L0["dy"]
    x0lo, y0lo = L0["x0"] - L0["dx"] / 2, L0["y0"] - L0["dy"] / 2
    c0 = (x0lo + ex0 / 2, y0lo + ey0 / 2)
    ex1, ey1 = L1["nx"] * L1["dx"], L1["ny"] * L1["dy"]
    x1lo, y1lo = L1["x0"] - L1["dx"] / 2, L1["y0"] - L1["dy"] / 2
    c1 = (x1lo + ex1 / 2, y1lo + ey1 / 2)
    fx = ex0 / 256.0
    fy = ey0 / 256.0

    clean_scene()
    scn = bpy.context.scene

    # --- water surfaces -----------------------------------------------------
    # The low grazing view lives off dark water + specular sky; a high-oblique
    # view sees almost no Fresnel reflection, so brighten the diffuse base and
    # widen the highlights or the relief reads flat from above.
    if a.cam_elev is not None:
        # same tint both levels: from above a tint step reads as a mismatched
        # rectangle (the orange frame already marks the window); no coat, its
        # sheen mutes the shading that carries the relief
        m_l0 = water_material("water_L0", (0.022, 0.09, 0.155), rough=0.12,
                              coat=0.0)
        m_l1 = water_material("water_L1", (0.022, 0.09, 0.155), rough=0.12,
                              coat=0.0)
    else:
        m_l0 = water_material("water_L0", (0.012, 0.055, 0.10))
        m_l1 = m_l0 if a.uniform_tint else \
            water_material("water_L1", (0.016, 0.066, 0.112))

    # high-resolution sea: 3x3 tile-aligned periodic images of L0
    ob0 = make_grid("L0_sea", 2304, 2304, c0, (3 * ex0, 3 * ey0), 0.0, m_l0)
    add_displace(ob0, seq_texture("etaL0",
                                  os.path.join(a.maps, "L0", "h0001.png"),
                                  repeat_x=3, repeat_y=3), strength)

    # far field: 9x9 periodic tiles on a lighter mesh, hole under the hi-res sea
    obo = make_grid("outer_sea", 1536, 1536, c0, (9 * ex0, 9 * ey0),
                    -0.006 * V, m_l0)
    add_displace(obo, seq_texture("etaL0_tiled",
                                  os.path.join(a.maps, "L0", "h0001.png"),
                                  repeat_x=9, repeat_y=9), strength)
    cut_hole(obo, x0lo - ex0 + L0["dx"], x0lo + 2 * ex0 - L0["dx"],
             y0lo - ey0 + L0["dy"], y0lo + 2 * ey0 - L0["dy"])

    # the child's accurate core, inset into a hole cut under its window;
    # the half-cell margin puts the rim step exactly under the frame bars,
    # hiding the step silhouette that top-down views would otherwise show
    ob1 = None
    if not a.l0_only:
        cut_hole(ob0, x1lo + L1["dx"] / 2, x1lo + ex1 - L1["dx"] / 2,
                 y1lo + L1["dy"] / 2, y1lo + ey1 - L1["dy"] / 2)
        # z-fighting guard over the hole rim: an amplitude piece for the
        # child-parent height mismatch plus a cell-size piece for
        # interpolation error. Larger lifts cast a visible silhouette line
        # on small domains.
        lift = 0.006 * V + 0.05 * L1["dx"]
        ob1 = make_grid("L1_sea", 512, 512, c1, (ex1, ey1), lift, m_l1)
        if a.cam_elev is not None:
            # plan views only: at grazing angles the skirt's exposed top
            # band renders as a dark curtain line along the rim (worse than
            # the sub-pixel open-edge artifact it cures from above)
            add_skirt(ob1, 0.05 * V + 2.0 * lift)
        add_displace(ob1, seq_texture("etaL1",
                                      os.path.join(a.maps, "L1", "h0001.png")),
                     strength)

    if a.cam_elev is not None and a.ripple:
        # cosmetic capillary ripple, real geometry (~1.5% of the wave height):
        # HOS fields are glassy below the grid scale, which reads as
        # out-of-focus from above; fine glitter restores the photographic
        # look. Only the fine meshes resolve it (the far field would alias).
        # A slowly drifting empty advects the pattern so it reads as water.
        bpy.ops.object.empty_add(location=(0, 0, 0))
        drift = bpy.context.active_object
        drift.name = "ripple_drift"
        # huge z-scale makes the 3D noise effectively 2D: the modifier stack
        # samples it after the wave displacement, so patch and sea would
        # otherwise get decorrelated ripples across the rim lift and the sea
        # would poke through as a dotted glint line
        drift.scale = (1.0, 1.0, 1000.0)
        drift.keyframe_insert("location", frame=1)
        drift.location = (-1.2, -0.4, 0)
        drift.keyframe_insert("location", frame=n)
        rip = bpy.data.textures.new("ripple", type="CLOUDS")
        rip.noise_scale = 0.07
        rip.noise_depth = 2
        for ob in ((ob0,) if ob1 is None else (ob0, ob1)):
            mod = ob.modifiers.new("ripple", "DISPLACE")
            mod.texture = rip
            mod.texture_coords = "OBJECT"
            mod.texture_coords_object = drift
            mod.direction = "Z"
            mod.mid_level = 0.5
            mod.strength = 2.0 * 0.015 * V

    # camera position (used by label parallax targeting below)
    if a.cam_elev is not None:
        # static high-oblique view from the -y side of the domain center
        th = math.radians(a.cam_elev)
        R = a.cam_dist * ex0
        cam_start = (c0[0], c0[1] - R * math.cos(th), R * math.sin(th))
    else:
        cam_start = (c0[0], y0lo - a.cam_back * ey0, a.cam_height * ex0)
    cam_end = (c0[0], y0lo - 0.203 * ey0, 0.078 * ex0)

    # --- annotation frames + labels (overlay pass only) ---------------------
    m_frame = emissive("frame", (1.0, 0.42, 0.05), 1.2)
    m_cell = emissive("cell", (0.85, 0.9, 1.0), 0.6)
    t = 0.0023 * ex0
    # plan views look straight down onto the rim seam; a wider bar blankets it
    bw = 3.0 * t if a.cam_elev is not None else t
    zbar = 1.1 * V + t
    for mat, (xl, yl, exr, eyr) in (
        (m_frame, (x1lo, y1lo, ex1, ey1)),
        (m_cell, (x0lo, y0lo, ex0, ey0)),
    ):
        for (cx, cy, sx, sy) in [
            (xl + exr / 2, yl, exr + bw, bw),
            (xl + exr / 2, yl + eyr, exr + bw, bw),
            (xl, yl + eyr / 2, bw, eyr + bw),
            (xl + exr, yl + eyr / 2, bw, eyr + bw),
        ]:
            bpy.ops.mesh.primitive_cube_add(size=1, location=(cx, cy, zbar))
            bar = bpy.context.active_object
            bar.scale = (sx, sy, t)
            bar.data.materials.append(mat)
            annotation_only(bar)

    def label_plane(png, pos, width, emit=1.8):
        img = bpy.data.images.load(os.path.abspath(png))
        aspect = img.size[1] / img.size[0]
        mat = bpy.data.materials.new(os.path.basename(png))
        mat.use_nodes = True
        nodes = mat.node_tree.nodes
        nodes.clear()
        texn = nodes.new("ShaderNodeTexImage")
        texn.image = img
        em = nodes.new("ShaderNodeEmission")
        em.inputs["Strength"].default_value = emit
        tr = nodes.new("ShaderNodeBsdfTransparent")
        mix = nodes.new("ShaderNodeMixShader")
        out = nodes.new("ShaderNodeOutputMaterial")
        lk = mat.node_tree.links
        lk.new(texn.outputs["Color"], em.inputs["Color"])
        lk.new(texn.outputs["Alpha"], mix.inputs["Fac"])
        lk.new(tr.outputs[0], mix.inputs[1])
        lk.new(em.outputs[0], mix.inputs[2])
        lk.new(mix.outputs[0], out.inputs["Surface"])
        if hasattr(mat, "surface_render_method"):
            mat.surface_render_method = "BLENDED"
        elif hasattr(mat, "blend_method"):
            mat.blend_method = "BLEND"
        bpy.ops.mesh.primitive_plane_add(size=1.0, location=pos)
        pl = bpy.context.active_object
        pl.scale = (width, width * aspect, 1.0)
        # tip the label to face the camera (72 deg suits the low framing)
        tilt = 90.0 - a.cam_elev if a.cam_elev is not None else 72.0
        pl.rotation_euler = (math.radians(tilt), 0.0, 0.0)
        pl.data.materials.append(mat)
        annotation_only(pl)

    label_dir = os.path.dirname(os.path.abspath(a.sky))
    zlab = 1.2 * V if a.cam_elev is not None else 1.9 * V
    # L0 label: mid-height of the clear strip between the south L0 and L1
    # edges, east of the child window, so it never crowds the L1 boundary.
    # L1 label: inside the window near its southwest corner (labels never
    # touch the boundary they name; the font color matches the frame).
    label_plane(os.path.join(label_dir, "label_L0.png"),
                (x0lo + 0.70 * ex0, y0lo + 0.5 * (y1lo - y0lo), zlab),
                0.086 * ex0)
    if a.cam_elev is not None:
        # plan view: bottom-center inside the window, shifted north by the
        # label height's parallax so it lands clear of the bar on screen,
        # and brighter to punch through the sun glitter
        ypar = zlab / math.tan(math.radians(a.cam_elev))
        label_plane(os.path.join(label_dir, "label_L1_stroke.png"),
                    (x1lo + 0.5 * ex1, y1lo + 0.10 * ey1 + ypar, zlab),
                    0.095 * ex0, emit=2.6)
    else:
        # at the BAR height the label is co-planar with the frame it sits in,
        # so both parallax-shift together and it reads inside the window from
        # any camera (an elevated label over a small box would slide across it)
        label_plane(os.path.join(label_dir, "label_L1.png"),
                    (x1lo + 0.18 * ex1, y1lo + 0.30 * ey1, zbar + 0.03 * V),
                    0.070 * ex0)

    # --- sky: emissive backdrop plane + world environment -------------------
    sky_img = bpy.data.images.load(os.path.abspath(a.sky))
    m_sky = bpy.data.materials.new("sky")
    m_sky.use_nodes = True
    nodes = m_sky.node_tree.nodes
    nodes.clear()
    texn = nodes.new("ShaderNodeTexImage")
    texn.image = sky_img
    em = nodes.new("ShaderNodeEmission")
    em.inputs["Strength"].default_value = 1.0
    out = nodes.new("ShaderNodeOutputMaterial")
    m_sky.node_tree.links.new(texn.outputs["Color"], em.inputs["Color"])
    m_sky.node_tree.links.new(em.outputs[0], out.inputs["Surface"])

    aspect = sky_img.size[0] / sky_img.size[1]
    sky_w = 6.64 * ex0
    sky_h = sky_w / aspect
    bpy.ops.mesh.primitive_plane_add(
        size=1.0, location=(c0[0], y0lo + 5.47 * ey0, sky_h / 2 - 0.047 * ex0))
    sky = bpy.context.active_object
    sky.rotation_euler = (math.radians(90), 0, 0)
    sky.scale = (sky_w, sky_h, 1)
    sky.data.materials.append(m_sky)
    if hasattr(sky, "visible_shadow"):
        sky.visible_shadow = False   # a backdrop must never shadow the sea

    world = bpy.data.worlds.new("w")
    scn.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    envt = world.node_tree.nodes.new("ShaderNodeTexEnvironment")
    envt.image = sky_img
    world.node_tree.links.new(envt.outputs["Color"], bg.inputs["Color"])
    bg.inputs["Strength"].default_value = 0.18 if a.cam_elev is not None else 0.45

    # --- light + camera -----------------------------------------------------
    bpy.ops.object.light_add(type="SUN",
                             location=(0.23 * ex0, -0.16 * ey0, 0.47 * ex0))
    sun = bpy.context.active_object
    if a.cam_elev is not None:
        # raking light: a low sun draws shading across the wave slopes,
        # which is what makes relief visible from a high viewpoint
        # beam from the northwest: the x-component shades each horseshoe's
        # steep leaning front (waves travel +x) while the southward component
        # forward-scatters sun glitter into the south-side camera
        sun.data.energy = 6.0
        sun.rotation_euler = (math.radians(66), 0.0, math.radians(-142))
    else:
        sun.data.energy = 3.5
        sun.rotation_euler = (math.radians(42), math.radians(-18),
                              math.radians(20))

    bpy.ops.object.empty_add(location=(c0[0], c0[1], a.aim_height * ex0))
    aim = bpy.context.active_object
    bpy.ops.object.camera_add(location=cam_start)
    cam = bpy.context.active_object
    cam.data.lens = 33
    cam.data.clip_end = 40.0 * ex0
    if not a.static_cam and a.cam_elev is None:
        cam.keyframe_insert("location", frame=1)
        cam.location = cam_end
        cam.keyframe_insert("location", frame=n)
        # default Bezier easing gives the slow-in/slow-out push
    tc = cam.constraints.new("TRACK_TO")
    tc.target = aim
    tc.track_axis = "TRACK_NEGATIVE_Z"
    tc.up_axis = "UP_Y"
    scn.camera = cam

    # --- render settings ----------------------------------------------------
    if a.use_cycles:
        scn.render.engine = "CYCLES"
        scn.cycles.samples = 48
        scn.cycles.use_denoising = True
        scn.cycles.device = "GPU"
        prefs = bpy.context.preferences.addons.get("cycles")
        if prefs:
            cp = prefs.preferences
            cp.compute_device_type = "METAL"
            cp.get_devices()
            for d in cp.devices:
                d.use = True
    else:
        for eng in ("BLENDER_EEVEE_NEXT", "BLENDER_EEVEE"):
            try:
                scn.render.engine = eng
                break
            except TypeError:
                continue
        if hasattr(scn, "eevee") and hasattr(scn.eevee, "use_raytracing"):
            scn.eevee.use_raytracing = True
        if a.cam_elev is not None:
            for attr, val in (("shadow_ray_count", 4), ("shadow_step_count", 8)):
                if hasattr(scn.eevee, attr):
                    setattr(scn.eevee, attr, val)
    scn.render.resolution_x, scn.render.resolution_y = a.res
    scn.render.fps = a.fps
    scn.frame_start, scn.frame_end = 1, n

    bpy.app.handlers.frame_change_pre.append(advance_heightmaps)

    if a.render_pass == "sea":
        for ob in ANNOTATIONS:
            ob.hide_render = True
    elif a.render_pass == "annotations":
        for ob in scn.objects:
            if ob not in ANNOTATIONS and ob.type in ("MESH", "FONT"):
                ob.hide_render = True
        scn.render.film_transparent = True

    os.makedirs(a.out, exist_ok=True)
    if a.save_blend:
        bpy.ops.wm.save_as_mainfile(filepath=os.path.abspath(a.save_blend))
        print("BLEND saved:", a.save_blend)
    if a.still is not None:
        scn.frame_set(a.still)
        advance_heightmaps(scn)
        scn.render.image_settings.file_format = "PNG"
        scn.render.filepath = os.path.join(a.out, f"still_{a.still:04d}.png")
        bpy.ops.render.render(write_still=True)
        print("STILL:", scn.render.filepath)
    if a.anim:
        # Blender 5.x has no FFMPEG movie format: render a PNG sequence and
        # encode externally (see tools/movie/README.md).
        try:
            scn.render.image_settings.file_format = "FFMPEG"
            scn.render.ffmpeg.format = "MPEG4"
            scn.render.ffmpeg.codec = "H264"
            scn.render.ffmpeg.constant_rate_factor = "HIGH"
            scn.render.filepath = os.path.join(a.out, "crescent_broadband_")
        except TypeError:
            scn.render.image_settings.file_format = "PNG"
            scn.render.filepath = os.path.join(a.out, "frames", "f")
        bpy.ops.render.render(animation=True)
        print("ANIM done ->", scn.render.filepath)


main()
