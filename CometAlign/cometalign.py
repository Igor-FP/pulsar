#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cometalign - Re-align a star-aligned FITS sequence onto a comet's nucleus.

Input frames are already registered on the stars (shared pixel grid). You mark
the comet nucleus on the FIRST and LAST frames (by capture time) in a pygame
view; every frame is then shifted by a time-linear interpolation of the comet's
apparent motion, so the comet stays fixed and the stars trail. The result is
ready to stack into a comet-locked image.

The two marks may instead be given on the CLI (--start X Y --stop X Y), which
skips the GUI (the user guarantees the coords belong to the earliest and latest
frame by time).
"""

import sys
import os
import csv
import numpy as np
from astropy.io import fits

# Import shared utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../lib")))
import batch_utils


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def usage():
    sys.stderr.write(
        "cometalign - Align a star-aligned FITS sequence on a comet's nucleus.\n"
        "\n"
        "Usage:\n"
        "  cometalign.py input_spec output_spec [options]\n"
        "\n"
        "  input_spec   - star-aligned frames: wildcard (*.fit), numbered\n"
        "                 (img0001.fit), or @list.txt\n"
        "  output_spec  - single file, numbered pattern, or directory\n"
        "\n"
        "Options:\n"
        "  --start X Y  - comet nucleus (px) on the EARLIEST frame\n"
        "  --stop  X Y  - comet nucleus (px) on the LATEST frame\n"
        "                 Give BOTH to skip the interactive GUI (you guarantee\n"
        "                 the coords belong to the first/last frame by time).\n"
        "  --ref F      - frame the comet is parked on: first (default) or last\n"
        "  --mtf M      - initial display midtones (default 0.05, smaller =\n"
        "                 brighter; GUI only)\n"
        "  --fill V     - value for out-of-footprint pixels after the shift\n"
        "                 (default 0 = no-data / black)\n"
        "\n"
        "GUI controls:\n"
        "  Left click      place the comet crosshair on the current frame\n"
        "  Tab             toggle FIRST <-> LAST frame (pan/zoom kept)\n"
        "  Arrows          pan the view (Shift = faster)\n"
        "  + / -           zoom in / out (mouse wheel zooms at the cursor)\n"
        "  Home / End      display brighter / darker (MTF)\n"
        "  Enter           confirm (both crosshairs set) and close\n"
        "  Q / Escape      cancel (writes nothing)\n"
    )
    sys.exit(1)


def _parse_xy(args, idx, name):
    if idx + 2 >= len(args):
        sys.stderr.write(f"Error: {name} requires two numbers: {name} X Y\n")
        sys.exit(1)
    try:
        return float(args[idx + 1]), float(args[idx + 2])
    except ValueError:
        sys.stderr.write(f"Error: {name} X Y must be numbers.\n")
        sys.exit(1)


def parse_args(argv):
    args = argv[1:]
    mtf_m = 0.05
    ref = "first"
    fill = 0.0
    start = None
    stop = None

    positional = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--start":
            start = _parse_xy(args, i, "--start"); i += 3; continue
        if a == "--stop":
            stop = _parse_xy(args, i, "--stop"); i += 3; continue
        if a == "--ref":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --ref needs first|last\n"); sys.exit(1)
            ref = args[i + 1].strip().lower()
            if ref not in ("first", "last"):
                sys.stderr.write("Error: --ref must be 'first' or 'last'\n"); sys.exit(1)
            i += 2; continue
        if a == "--mtf":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --mtf needs a value\n"); sys.exit(1)
            try:
                mtf_m = float(args[i + 1])
            except ValueError:
                sys.stderr.write("Error: --mtf must be a number\n"); sys.exit(1)
            if not (0.0 < mtf_m < 1.0):
                sys.stderr.write("Error: --mtf must be in (0..1)\n"); sys.exit(1)
            i += 2; continue
        if a == "--fill":
            if i + 1 >= len(args):
                sys.stderr.write("Error: --fill needs a value\n"); sys.exit(1)
            try:
                fill = float(args[i + 1])
            except ValueError:
                sys.stderr.write("Error: --fill must be a number\n"); sys.exit(1)
            i += 2; continue
        if a.startswith("--"):
            sys.stderr.write(f"Error: unknown option: {a}\n"); usage()
        positional.append(a); i += 1

    if len(positional) < 2:
        usage()
    if (start is None) != (stop is None):
        sys.stderr.write("Error: give BOTH --start and --stop, or neither.\n")
        sys.exit(1)

    return positional[0], positional[1], mtf_m, ref, fill, start, stop


# ---------------------------------------------------------------------------
# Observation time from header  (returned as JD days; unit cancels in the ratio)
# ---------------------------------------------------------------------------

def get_frame_time(header):
    """Return a comparable observation time (JD days), or None if absent."""
    if "JD" in header:
        try:
            return float(header["JD"])
        except (TypeError, ValueError):
            pass
    if "MJD-OBS" in header:
        try:
            return float(header["MJD-OBS"]) + 2400000.5
        except (TypeError, ValueError):
            pass
    for key in ("DATE-OBS", "DATE_OBS", "DATE"):
        raw = header.get(key)
        if not raw:
            continue
        s = str(raw).strip()
        if "T" not in s:
            tobs = header.get("TIME-OBS")
            if tobs:
                s = s + "T" + str(tobs).strip()
            elif " " in s:
                s = s.replace(" ", "T", 1)
        try:
            from astropy.time import Time
            return float(Time(s, format="isot", scale="utc").jd)
        except Exception:
            continue
    return None


def read_times(file_list):
    """Map each file to its JD time; raise if any frame lacks a usable time."""
    times = {}
    for f in file_list:
        with fits.open(f, memmap=False) as hdul:
            t = get_frame_time(hdul[0].header)
        if t is None:
            raise RuntimeError(
                f"No usable time (JD / MJD-OBS / DATE-OBS) in header of '{f}'")
        times[f] = t
    return times


# ---------------------------------------------------------------------------
# Display stretch: percentile clip + MTF  (shared shape with align.py)
# ---------------------------------------------------------------------------

def mtf_stretch(x, m):
    """Midtones Transfer Function; x in [0..1], m in (0..1), m<0.5 brightens."""
    denom = np.maximum(m + x * (1.0 - 2.0 * m), 1e-10)
    return (1.0 - m) * x / denom


def load_raw(filepath):
    """Load raw FITS data. Returns (raw_array, img_h, img_w)."""
    with fits.open(filepath, memmap=False) as hdul:
        raw = hdul[0].data
        if raw is None:
            raise ValueError(f"No primary image data in {filepath}")
        raw = raw.copy()
    if raw.ndim == 3:
        return raw, raw.shape[1], raw.shape[2]
    if raw.ndim == 2:
        return raw, raw.shape[0], raw.shape[1]
    raise ValueError(f"Unsupported FITS shape: {raw.shape}")


def _percentile_bounds(data):
    """Display black/white points ([0.01%, 99.99%]). Computed once per frame and
    cached, so changing the MTF does not re-run the (slow) percentile."""
    flat = data.ravel().astype(np.float64)
    vmin = float(np.percentile(flat, 0.01))
    vmax = float(np.percentile(flat, 99.99))
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def compute_bounds(raw):
    """Per-channel percentile bounds for a frame (list of (vmin, vmax))."""
    if raw.ndim == 3:
        return [_percentile_bounds(raw[ch]) for ch in range(min(3, raw.shape[0]))]
    return [_percentile_bounds(raw)]


def _stretch_channel(data, mtf_m, bounds):
    """Normalize with cached bounds + MTF -> uint8. float32 for speed (the MTF
    re-stretch runs on every Home/End press)."""
    vmin, vmax = bounds
    work = (data.astype(np.float32) - np.float32(vmin)) / np.float32(vmax - vmin)
    np.clip(work, 0.0, 1.0, out=work)
    work = mtf_stretch(work, np.float32(mtf_m))
    np.multiply(work, np.float32(255.0), out=work)
    return work.astype(np.uint8)


def stretch_to_rgb(raw, mtf_m, bounds):
    """Percentile clip (cached bounds) + MTF -> rgb_uint8 (H,W,3)."""
    if raw.ndim == 3:
        img_h, img_w = raw.shape[1], raw.shape[2]
        rgb = np.zeros((img_h, img_w, 3), dtype=np.uint8)
        for ch in range(min(3, raw.shape[0])):
            rgb[:, :, ch] = _stretch_channel(raw[ch], mtf_m, bounds[ch])
        return rgb
    gray = _stretch_channel(raw, mtf_m, bounds[0])
    return np.stack([gray, gray, gray], axis=-1)


# ---------------------------------------------------------------------------
# PyGame GUI: mark the comet nucleus on the first and last frames
# ---------------------------------------------------------------------------

def run_gui(first_path, last_path, mtf_m, start_pt=None, stop_pt=None):
    """Interactive marking. Returns ((x0,y0),(x1,y1)) or None if cancelled."""
    import pygame

    raw0, h0, w0 = load_raw(first_path)
    raw1, h1, w1 = load_raw(last_path)
    if (h0, w0) != (h1, w1):
        raise RuntimeError(
            f"First/last frame dimensions differ ({w0}x{h0} vs {w1}x{h1}); "
            f"frames must be star-aligned to a common grid.")
    img_w, img_h = w0, h0
    raws = [raw0, raw1]
    paths = [first_path, last_path]
    labels = ["FIRST", "LAST"]
    points = [list(start_pt) if start_pt else None,
              list(stop_pt) if stop_pt else None]
    cur = 0

    pygame.init()
    info = pygame.display.Info()
    win_w = min(1600, info.current_w - 100)
    win_h = min(1000, info.current_h - 100)
    screen = pygame.display.set_mode((win_w, win_h), pygame.RESIZABLE)
    pygame.display.set_caption("cometalign")
    pygame.key.set_repeat(150, 30)
    font = pygame.font.SysFont("consolas,courier,monospace", 16)
    hint_font = pygame.font.SysFont("consolas,courier,monospace", 22)

    status_h = 52  # two text lines: status + hotkey legend
    surfaces = [None, None]
    surf_mtf = [None, None]
    bounds_cache = [None, None]  # per-frame percentile bounds (computed once)

    # Transient hint overlay: shown 5 s solid, then a 10 s fade-out; never blocks.
    hint_text = None
    hint_t0 = 0

    def set_hint(text):
        nonlocal hint_text, hint_t0
        hint_text = text
        hint_t0 = pygame.time.get_ticks()

    # View state (shared across both frames; they share the grid)
    zoom = max(0.02, min(min(win_w / img_w, (win_h - status_h) / img_h), 1.0))
    view_cx = img_w / 2.0
    view_cy = img_h / 2.0

    MTF_STEP = 1.4
    ZOOM_STEP = 1.25
    MIN_ZOOM, MAX_ZOOM = 0.02, 40.0

    nav_keys = {pygame.K_TAB, pygame.K_RETURN, pygame.K_KP_ENTER,
                pygame.K_HOME, pygame.K_END}
    nav_held = set()

    def ensure_surface(i):
        if bounds_cache[i] is None:
            bounds_cache[i] = compute_bounds(raws[i])
        if surfaces[i] is None or surf_mtf[i] != mtf_m:
            rgb = stretch_to_rgb(raws[i], mtf_m, bounds_cache[i])
            surfaces[i] = pygame.surfarray.make_surface(
                np.ascontiguousarray(np.transpose(rgb, (1, 0, 2))))
            surf_mtf[i] = mtf_m

    def img_to_screen(ix, iy):
        uh = win_h - status_h
        return (win_w / 2.0 + (ix - view_cx) * zoom,
                uh / 2.0 + (iy - view_cy) * zoom)

    def screen_to_img(sx, sy):
        uh = win_h - status_h
        return (view_cx + (sx - win_w / 2.0) / zoom,
                view_cy + (sy - uh / 2.0) / zoom)

    def draw_mark(p, color, r):
        sx, sy = img_to_screen(p[0], p[1])
        sx, sy = int(round(sx)), int(round(sy))
        pygame.draw.circle(screen, color, (sx, sy), r, 1)
        pygame.draw.line(screen, color, (sx - r - 5, sy), (sx + r + 5, sy), 1)
        pygame.draw.line(screen, color, (sx, sy - r - 5), (sx, sy + r + 5), 1)

    def render():
        screen.fill((25, 25, 25))
        uh = win_h - status_h
        ensure_surface(cur)
        surf = surfaces[cur]

        half_w = (win_w / 2.0) / zoom
        half_h = (uh / 2.0) / zoom
        cx0 = max(0, int(np.floor(view_cx - half_w)))
        cy0 = max(0, int(np.floor(view_cy - half_h)))
        cx1 = min(img_w, int(np.ceil(view_cx + half_w)))
        cy1 = min(img_h, int(np.ceil(view_cy + half_h)))
        if cx1 > cx0 and cy1 > cy0:
            sub = surf.subsurface(pygame.Rect(cx0, cy0, cx1 - cx0, cy1 - cy0))
            dw = max(1, int(round((cx1 - cx0) * zoom)))
            dh = max(1, int(round((cy1 - cy0) * zoom)))
            scaled = pygame.transform.scale(sub, (dw, dh))
            dsx, dsy = img_to_screen(cx0, cy0)
            screen.blit(scaled, (int(round(dsx)), int(round(dsy))))

        other = 1 - cur
        if points[0] is not None and points[1] is not None:
            a = img_to_screen(*points[0])
            b = img_to_screen(*points[1])
            pygame.draw.line(screen, (0, 140, 200),
                             (int(a[0]), int(a[1])), (int(b[0]), int(b[1])), 1)
        if points[other] is not None:
            draw_mark(points[other], (0, 120, 180), 6)
        if points[cur] is not None:
            draw_mark(points[cur], (255, 60, 60), 9)

        pygame.draw.rect(screen, (0, 0, 0), (0, win_h - status_h, win_w, status_h))
        p = points[cur]
        pstr = f"({p[0]:.1f},{p[1]:.1f})" if p is not None else "unset"
        both = ("BOTH SET - press Enter to confirm"
                if (points[0] is not None and points[1] is not None)
                else "click the nucleus")
        info = (f"{labels[cur]}  {os.path.basename(paths[cur])}  comet={pstr}  "
                f"zoom={zoom:.2f}x  mtf={mtf_m:.4f}  [{both}]")
        keys = ("LMB: nucleus | Tab: frame | Arrows: pan (Shift=faster) | "
                "+/- or wheel: zoom | Home/End: brighter/darker | "
                "Enter: confirm | Q/Esc: cancel")
        screen.blit(font.render(info, True, (215, 215, 215)),
                    (6, win_h - status_h + 4))
        screen.blit(font.render(keys, True, (140, 175, 140)),
                    (6, win_h - status_h + 27))

        # transient hint overlay (non-blocking; 5 s solid, then 10 s fade-out)
        if hint_text is not None:
            elapsed = (pygame.time.get_ticks() - hint_t0) / 1000.0
            if elapsed < 5.0:
                alpha = 255
            elif elapsed < 15.0:
                alpha = int(255 * (1.0 - (elapsed - 5.0) / 10.0))
            else:
                alpha = 0
            if alpha > 0:
                hsurf = hint_font.render(hint_text, True, (255, 245, 180))
                hsurf.set_alpha(alpha)
                pad = 14
                bw = hsurf.get_width() + 2 * pad
                bh = hsurf.get_height() + pad
                box = pygame.Surface((bw, bh), pygame.SRCALPHA)
                box.fill((0, 0, 0, int(alpha * 0.55)))
                bx = max(0, (win_w - bw) // 2)
                screen.blit(box, (bx, 28))
                screen.blit(hsurf, (bx + pad, 28 + pad // 2))

        pygame.display.flip()

    set_hint("Mark the comet nucleus on the FIRST frame")
    render()

    result = None
    running = True
    clock = pygame.time.Clock()
    need_render = False

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False  # cancel
            elif event.type == pygame.VIDEORESIZE:
                win_w, win_h = event.size
                screen = pygame.display.set_mode((win_w, win_h), pygame.RESIZABLE)
                need_render = True
            elif event.type == pygame.MOUSEBUTTONDOWN:
                mx, my = event.pos
                if event.button == 1:
                    if my < win_h - status_h:
                        ix, iy = screen_to_img(mx, my)
                        if 0 <= ix < img_w and 0 <= iy < img_h:
                            points[cur] = [ix, iy]
                            if points[0] is not None and points[1] is not None:
                                set_hint("Both nuclei set - press Enter to confirm")
                            elif cur == 0:
                                set_hint("Tab - switch to the LAST frame")
                            else:
                                set_hint("Tab - go back and mark the FIRST frame")
                            need_render = True
                elif event.button in (4, 5):  # wheel: zoom at cursor
                    bx, by = screen_to_img(mx, my)
                    if event.button == 4:
                        zoom = min(MAX_ZOOM, zoom * ZOOM_STEP)
                    else:
                        zoom = max(MIN_ZOOM, zoom / ZOOM_STEP)
                    uh = win_h - status_h
                    view_cx = bx - (mx - win_w / 2.0) / zoom
                    view_cy = by - (my - uh / 2.0) / zoom
                    need_render = True
            elif event.type == pygame.KEYUP:
                nav_held.discard(event.key)
            elif event.type == pygame.KEYDOWN:
                if event.key in nav_keys:
                    if event.key in nav_held:
                        continue
                    nav_held.add(event.key)
                shift = bool(event.mod & pygame.KMOD_SHIFT)
                pan = (200.0 if shift else 60.0) / zoom

                if event.key in (pygame.K_q, pygame.K_ESCAPE):
                    running = False  # cancel
                elif event.key == pygame.K_TAB:
                    cur = 1 - cur
                    if cur == 1:
                        set_hint("Mark the comet nucleus on the LAST frame, "
                                 "or Tab to go back to the FIRST")
                    else:
                        set_hint("Mark the comet nucleus on the FIRST frame, "
                                 "or Tab for the LAST")
                    need_render = True
                elif event.key == pygame.K_LEFT:
                    view_cx -= pan; need_render = True
                elif event.key == pygame.K_RIGHT:
                    view_cx += pan; need_render = True
                elif event.key == pygame.K_UP:
                    view_cy -= pan; need_render = True
                elif event.key == pygame.K_DOWN:
                    view_cy += pan; need_render = True
                elif event.key in (pygame.K_PLUS, pygame.K_EQUALS, pygame.K_KP_PLUS):
                    zoom = min(MAX_ZOOM, zoom * ZOOM_STEP); need_render = True
                elif event.key in (pygame.K_MINUS, pygame.K_KP_MINUS):
                    zoom = max(MIN_ZOOM, zoom / ZOOM_STEP); need_render = True
                elif event.key == pygame.K_HOME:
                    mtf_m = max(0.001, mtf_m / MTF_STEP); need_render = True
                elif event.key == pygame.K_END:
                    mtf_m = min(0.95, mtf_m * MTF_STEP); need_render = True
                elif event.key in (pygame.K_RETURN, pygame.K_KP_ENTER):
                    if points[0] is not None and points[1] is not None:
                        result = (tuple(points[0]), tuple(points[1]))
                        running = False
                    else:
                        need_render = True

        # Animate the fade-out (5..15 s) by re-rendering each frame; the solid
        # 0..5 s phase is static (already drawn); clear the hint after 15 s.
        if hint_text is not None:
            el = pygame.time.get_ticks() - hint_t0
            if el >= 15000:
                hint_text = None
                need_render = True
            elif el >= 5000:
                need_render = True

        if need_render:
            render()
            need_render = False
        clock.tick(60)

    pygame.quit()
    return result


# ---------------------------------------------------------------------------
# Shift + output
# ---------------------------------------------------------------------------

def shift_frame(data, dx, dy, fill):
    """Subpixel translate by (dx, dy) px (x=col, y=row); RGB per channel;
    out-of-footprint -> fill. Returns float64."""
    from scipy.ndimage import shift as ndi_shift
    work = np.asarray(data, dtype=np.float64)
    if work.ndim == 3:
        out = np.empty_like(work)
        for c in range(work.shape[0]):
            out[c] = ndi_shift(work[c], (dy, dx), order=3, mode="constant",
                               cval=float(fill), prefilter=True)
        return out
    return ndi_shift(work, (dy, dx), order=3, mode="constant",
                     cval=float(fill), prefilter=True)


def to_output(data, orig_dtype):
    """Sanitize NaN/Inf and cast back to the original dtype (clamp+round int)."""
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    if np.issubdtype(orig_dtype, np.integer):
        info = np.iinfo(orig_dtype)
        return np.clip(np.rint(data), info.min, info.max).astype(orig_dtype)
    return data.astype(orig_dtype)


def process_one(infile, outfile, dx, dy, fill, ref):
    with fits.open(infile, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
        if data is None:
            raise RuntimeError("no image data")
        orig = data.dtype
        data = data.copy()

    if abs(dx) < 1e-4 and abs(dy) < 1e-4:
        out = data  # reference frame: leave untouched (no needless resample)
    else:
        out = to_output(shift_frame(data, dx, dy, fill), orig)

    header["HISTORY"] = (f"cometalign.py: comet-aligned (ref={ref}), "
                         f"shift dx={dx:.3f} dy={dy:.3f} px")
    out_dir = os.path.dirname(outfile)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    fits.PrimaryHDU(out, header=header).writeto(outfile, overwrite=True)


def save_points_sidecar(output_spec, first, last, t0, t1, p0, p1):
    base = output_spec.replace("*", "").replace("?", "")
    if base.endswith((".fit", ".fits")):
        base = os.path.splitext(base)[0]
    if base.endswith(("/", "\\")) or os.path.isdir(base):
        base = os.path.join(base, "comet")
    if not base:
        base = "comet"
    path = base + ".cometpts.csv"
    try:
        d = os.path.dirname(path)
        if d and not os.path.exists(d):
            os.makedirs(d, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["frame", "basename", "jd", "x", "y"])
            w.writerow(["first", os.path.basename(first), f"{t0:.8f}",
                        f"{p0[0]:.3f}", f"{p0[1]:.3f}"])
            w.writerow(["last", os.path.basename(last), f"{t1:.8f}",
                        f"{p1[0]:.3f}", f"{p1[1]:.3f}"])
        sys.stderr.write(f"Saved marks -> {path}\n")
    except Exception as e:
        sys.stderr.write(f"Warning: could not save sidecar CSV: {e}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    input_spec, output_spec, mtf_m, ref, fill, start, stop = parse_args(sys.argv)

    io_pairs = batch_utils.build_io_file_lists(input_spec, output_spec)
    if len(io_pairs) < 2:
        sys.stderr.write("Error: need at least 2 frames for comet alignment.\n")
        sys.exit(1)

    input_files = [inp for inp, _ in io_pairs]
    times = read_times(input_files)

    ordered = sorted(input_files, key=lambda f: times[f])
    first, last = ordered[0], ordered[-1]
    t0, t1 = times[first], times[last]
    if t1 <= t0:
        sys.stderr.write("Error: first and last frames share a timestamp; "
                         "cannot interpolate the comet motion.\n")
        sys.exit(1)

    sys.stderr.write(f"{len(input_files)} frames; session span "
                     f"{(t1 - t0) * 24.0:.3f} h\n")
    sys.stderr.write(f"  first: {os.path.basename(first)}\n")
    sys.stderr.write(f"  last : {os.path.basename(last)}\n")

    if start is not None and stop is not None:
        p0 = np.array(start, dtype=float)
        p1 = np.array(stop, dtype=float)
        sys.stderr.write(f"CLI marks: start={tuple(p0)} stop={tuple(p1)}\n")
    else:
        res = run_gui(first, last, mtf_m)
        if res is None:
            sys.stderr.write("Cancelled; nothing written.\n")
            sys.exit(1)
        p0 = np.array(res[0], dtype=float)
        p1 = np.array(res[1], dtype=float)

    save_points_sidecar(output_spec, first, last, t0, t1, p0, p1)
    sys.stderr.write(
        f"Reproduce without the GUI: "
        f"--start {p0[0]:.2f} {p0[1]:.2f} --stop {p1[0]:.2f} {p1[1]:.2f}\n")

    ref_pos = p0 if ref == "first" else p1
    motion = p1 - p0

    total = len(io_pairs)
    for i, (infile, outfile) in enumerate(io_pairs, start=1):
        try:
            f = (times[infile] - t0) / (t1 - t0)
            comet = p0 + f * motion
            dx = float(ref_pos[0] - comet[0])
            dy = float(ref_pos[1] - comet[1])
            process_one(infile, outfile, dx, dy, fill, ref)
            sys.stderr.write(
                f"\r[{i}/{total}] {os.path.basename(outfile)}  "
                f"shift=({dx:+.2f},{dy:+.2f})      ")
            sys.stderr.flush()
        except Exception as e:
            sys.stderr.write(f"\nError on '{infile}': {e}\n")

    sys.stderr.write("\nDone.\n")


if __name__ == "__main__":
    main()
