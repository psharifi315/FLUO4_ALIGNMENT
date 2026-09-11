"""
===============================================================================
align_v2.py
Per-well TH mask -> Fluo-4 alignment, scored, with an auditable pass/fail per
well and the region growth moved OUT of Harmony and INTO the mask.
===============================================================================
WHY THIS REPLACES run_batch_map2_anchor

  1  MAP2 was the wrong anchor. Its dominant structure is the neurite mesh;
     Fluo-4's is somata. ECC was matching a web against a field of blobs, and
     on KR25.14 P1 it failed on 2 of 4 wells returning tx = +148 and -133 px.
     Those two failures were then excluded and the mean of the surviving two -
     which disagreed by 13.3 px - was applied to all 36 wells.

  2  This script registers the TH mask DIRECTLY onto the Fluo-4 image, by
     maximising the mean Fluo-4 intensity underneath the mask. TH somata are
     bright in Fluo-4, so the correct position is the one where every mask
     object sits on a cell. Validated on KR25.15 P1 R2C3: recovered (-32, -79)
     px, identical to the shift derived independently by mask IoU, with
     enrichment 2.27 at the peak against 0.97 unaligned.

  3  The whole search is one FFT cross-correlation (80 ms/well), so EVERY well
     gets its own transform. No plate mean, no anchor requirement, no ECC.

  4  Label IDs are preserved. read_mask_rgba_binary collapsed uint16 labels to
     0/255, which merged touching somata: 32 labels became 28 components on the
     file I checked, a 12.5% object loss before Harmony ever saw the mask.

  5  Every well gets three QC numbers and a verdict, so the wells you cannot
     rescue are DECLARED rather than silently averaged in.

WHAT YOU NEED THAT YOU MAY NOT HAVE YET
  One Fluo-4 intensity image per well (field 1, timepoint 1, or a projection).
  In Harmony this is a single batch image export per measurement, not 36 manual
  ones. Everything else - the TH label masks - you already have.

SETTING THE HARMONY RESIZE TO 0%
  Set CALIBRATE_DILATION = True and run once. The script reports, per plate,
  what fraction of TH mask objects contain a bright Fluo-4 cell as a function of
  outward growth in pixels. Read the smallest growth at which that curve
  plateaus, put it in DILATE_PX, re-run, and the masks written to for_harmony/
  already carry that growth. Harmony's Outer Border then stays at 0% and the
  growth is an explicit, documented number in your methods instead of a hidden
  -100% in a building block.
===============================================================================
"""

# ============================== USER CONFIG ==============================
# ONE ENTRY PER PLATE. Wells are discovered from the R#C# in the filenames -
# you never list wells. Add plates by adding dicts; comment one out to skip it.
#
#   th_mask : the ORIGINAL Harmony THPOSSOMA export for that plate, BEFORE any
#             alignment. Do NOT point this at for_harmony_full - those masks
#             have already been plate-shifted and you would be correcting a
#             correction. The script warns if it thinks that has happened.
#   fluo    : folder of static Fluo-4 images (make_statics.py output).
#   th_icc  : optional, visual QC panel only. "" to skip.
# ---- the only three things to edit ------------------------------------------
# 1. where the per-plate TH mask folders live (the PARENT of the UUID folders)
MASK_ROOT    = r"/content/drive/MyDrive/Results/Phenix/FLUO4_ALIGNMENT/FULL_TH_MASKS/thpossoma2"
STATICS_ROOT = r"/content/drive/MyDrive/Results/Phenix/FLUO4_ALIGNMENT/RAW_ICC_FLUO4_TIMESERIES/STATICS"
OUTPUT_ROOT  = r"/content/drive/MyDrive/Results/Phenix/FLUO4_ALIGNMENT/ALIGN_FINAL_THPOS"

JOBS = {
    "KR25.14P1":   ("1d973782-5b47-4627-b4b1-5181a5882532", "KR25.14P1_LIVE"),
    "KR25.14P2":   ("4ab9295f-c8da-4953-8562-d3ead70a9740", "KR25.14P2_LIVE"),
    "KR25.15P1":   ("5a5434d2-310b-41f4-baac-8e61eeb52bf9", "KR25.15P1_LIVE"),
    "KR25.15P2":   ("3b0e112c-dba0-4fe9-9418-917dbdf23161", "KR25.15P2_LIVE"),
    "KR25.15P3":   ("3ea34d48-d2df-4312-a563-d7cd4a59cbdf", "KR25.15P3_LIVE"),
    "KR25.15P5":   ("955dc245-ceae-47af-9968-4462eff47c54", "KR25.15P5_LIVE"),
    "KR25.21":     ("ce1f47b9-903d-4dba-965b-33f5b6740f71", "KR25.21_LIVE"),
    "KR25.22P1":   ("c4c11143-0f86-472a-ae38-87380ff453b2", "KR25.22P1_LIVE"),
    "KR25.15_384": ("b532145c-8f3a-453f-9ca6-d1319c14bad8", "KR25.15_384_LIVE"),
}

MASK_REQUIRE        = r"F1T1P0_TH"
FLUO_REQUIRE        = r"_fluo_mean"
MAX_SHIFT_PX        = 300
SUBTRACT_BACKGROUND = True
MIN_PEAK_RATIO      = 1.30
DILATE_PX           = 0
CALIBRATE_DILATION  = False
ICC                 = {}

# Optional per-plate POSTHOC ICC statics folder (make_statics.py output, files
# named <well>_ch1.tif ... <well>_ch4.tif). Adding this turns on a SECOND
# contact sheet that is the stronger scientific check: the Fluo-4 comparison
# shows the mask lands on CELLS; the ICC comparison shows it lands on the cells
# that are actually TH-POSITIVE, which is what the thesis claim needs.
#
# You do not have to know which channel is TH. The mask was DERIVED from the TH
# channel, so the TH channel is the one where the aligned mask sits on the
# brightest pixels. The script scores every channel and picks it, then reports
# the scores so you can see the margin.
ICC_CHANNEL = 0   # 0 = auto-detect from the enrichment scores. Set 1-4 to force.
ICC_DETECT_N = 5  # wells sampled for the auto-detection

# One folder; each plate gets its own subfolder inside.
OUTPUT_ROOT = r"/content/drive/MyDrive/Results/Phenix/FLUO4_ALIGNMENT/ALIGN_V4"


# QC thresholds. Run once, look at the printed distribution, then set these.
MIN_ENRICH     = 1.50    # only used when SUBTRACT_BACKGROUND = False
MIN_PEAK_RATIO = 1.30    # peak / best rival peak outside EXCLUDE_R. THE metric:
                         # size-, density- and background-independent. With
                         # background subtraction on, a real match sits well
                         # above 2 and chance sits at ~1.0.
BRIGHT_PCTILE  = 75      # a mask object is "occupied" if its mean Fluo-4
                         # exceeds this percentile of the field
PIXEL_UM       = 0.975
# =========================================================================

import os, re, json
from pathlib import Path
import numpy as np
import cv2
import tifffile as tiff


# ------------------------------- io --------------------------------------
def read_tif_any(path):
    try:
        a = tiff.imread(str(path))
    except Exception:
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None
        with Image.open(str(path)) as im:
            im.seek(0); a = np.array(im)
    if a.ndim == 3:
        a = a[..., 3] if a.shape[-1] == 4 else (
            a[..., 0] if a.shape[-1] in (3, 4) else a[0])
    return a


def read_labels(path):
    """TH mask as an integer label image. Labels are NOT collapsed to binary."""
    a = read_tif_any(path)
    if a.dtype == bool:
        a = a.astype(np.uint16)
    if len(np.unique(a[a > 0])) <= 1:      # already binary: relabel so we
        from scipy import ndimage as ndi   # at least recover components
        a, _ = ndi.label(a > 0)
    return a.astype(np.int32)


def to_float(a):
    a = np.asarray(a, np.float64)
    a[~np.isfinite(a)] = 0
    return a


# Drive FUSE charges roughly a second per directory listing, so os.walk over a
# 36-well mask tree (each well holding MAP/DAPI) is ~108 listings = minutes of
# silence. Instead: ONE listing to find the well entries, then resolve the tif
# only for the wells actually processed.
def index_masks(root, require=None):
    """{well: file-or-directory}. Handles both a tree of <well>/MAP/DAPI/*.tif
    and a flat folder of per-well tifs."""
    root = Path(root)
    if not root.is_dir():
        return {}
    wells = {}
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        if require and not re.search(require, d.name, re.I):
            continue
        m = re.search(r"(R\d+C\d+)", d.name, re.I)
        if m:
            wells[m.group(1).upper()] = d
    if wells:
        return wells
    for p in sorted(root.glob("*.tif*")):          # flat layout
        if require and not re.search(require, p.name, re.I):
            continue
        m = re.search(r"(R\d+C\d+)", p.name, re.I)
        if m:
            wells[m.group(1).upper()] = p
    if wells:
        return wells

    # Nothing at this level. The commonest cause is pointing th_mask at a PARENT
    # of per-plate folders (e.g. TH_masks/THPOSSOMA, which holds one UUID folder
    # per plate, or a <plate>_THPOSSOMA folder holding for_harmony_full etc).
    # Every plate uses the same R#C# names, so mixing them silently pairs each
    # well with an arbitrary plate's mask - which registers at chance and fails
    # every well. Name the candidates instead of returning nothing.
    kids = []
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        n = len(index_masks(d, require))
        if n:
            kids.append((d.name, n))
    if kids:
        print(f"    !! '{root.name}' contains no wells itself, but these "
              f"subfolders do:")
        for n, c in kids[:25]:
            print(f"         {n}   ({c} wells)")
        print("    !! th_mask must point at ONE of those, not at the parent.")
        print("    !! Every plate reuses the same R#C# names, so a parent path")
        print("    !! pairs each well with an arbitrary plate and fails all of them.")
    return {}


def resolve_mask(entry):
    """Turn a well entry into the actual tif. Cached by the caller's dict."""
    entry = Path(entry)
    if entry.is_file():
        return entry
    hits = sorted(entry.glob("*/*/*.tif*")) or sorted(entry.glob("*.tif*"))
    return hits[0] if hits else None


def index_flat(root, require=None):
    """One listing of a flat folder of statics."""
    out = {}
    root = Path(root)
    if not root.is_dir():
        return out
    for p in sorted(root.glob("*.tif*")):
        if require and not re.search(require, p.name, re.I):
            continue
        m = re.search(r"(R\d+C\d+)", p.name, re.I)
        if m:
            out[m.group(1).upper()] = p
    return out


# --------------------------- the objective -------------------------------
def align(F, B, max_shift=MAX_SHIFT_PX, exclude_r=EXCLUDE_R):
    """
    Find the (dx, dy) that maximises mean Fluo-4 intensity under the mask.

    sum(F * roll(B, dy, dx)) is the cross-correlation of F and B, so the whole
    search is one FFT rather than a grid of trials. Verified identical to brute
    force on KR25.15 P1 R2C3.

    Returns dict with dx, dy, enrich (at peak), enrich0 (unshifted),
    peak_ratio (peak / best rival outside exclude_r).
    """
    H, W = F.shape
    n = int(B.sum())
    if n == 0:
        return None
    Fm = float(F.mean())

    # BACKGROUND SUBTRACTION BEFORE THE CORRELATION.
    # The raw objective is sum(F under mask) / (n * mean F). If F carries a
    # large constant background - which a DENSE, brightly stained culture does -
    # then that sum is roughly n * mean(F) wherever the mask is put, enrichment
    # collapses towards 1, and the peak flattens. KR25.22 P1 showed exactly
    # this: 80 mask objects per well, enrichment only 1.06 -> 1.42, peak/rival
    # 1.02. Not a sparse mask and not a bad plate - a swamped objective.
    # Subtracting the mean makes the correlation respond to the cell PATTERN
    # rather than the background LEVEL. Measured on synthetic fields:
    #     sparse, low  background   pk/riv  2.17 -> 12.99
    #     dense,  high background   pk/riv  1.24 ->  5.06
    # and the correct shift is recovered in both cases either way - it is the
    # SHARPNESS of the peak that changes, which is what the verdict keys on.
    Fc = F - Fm if SUBTRACT_BACKGROUND else F

    C = np.fft.irfft2(np.fft.rfft2(Fc) * np.conj(np.fft.rfft2(B.astype(np.float64))),
                      s=(H, W))
    E = C / (n * abs(Fm) if Fm else 1.0)          # enrichment surface, [dy, dx]

    ys, xs = np.mgrid[0:H, 0:W]
    DY = np.where(ys > H // 2, ys - H, ys)
    DX = np.where(xs > W // 2, xs - W, xs)
    win = (np.abs(DY) <= max_shift) & (np.abs(DX) <= max_shift)

    iy, ix = np.unravel_index(np.argmax(np.where(win, E, -np.inf)), E.shape)
    dx, dy = int(DX[iy, ix]), int(DY[iy, ix])
    peak = float(E[iy, ix])
    r = np.hypot(DY - dy, DX - dx)
    rival = float(np.where(win & (r > exclude_r), E, -np.inf).max())
    return dict(dx=dx, dy=dy, enrich=round(peak, 3),
                enrich0=round(float(E[0, 0]), 3),
                rival=round(rival, 3),
                peak_ratio=round(peak / rival, 3) if rival > 0 else np.inf)


def shift_labels(lab, dx, dy):
    """Translate a label image with zero fill. Sign matches align()'s roll."""
    M = np.float32([[1, 0, dx], [0, 1, dy]])
    return cv2.warpAffine(lab.astype(np.int32), M, (lab.shape[1], lab.shape[0]),
                          flags=cv2.INTER_NEAREST, borderMode=cv2.BORDER_CONSTANT,
                          borderValue=0).astype(np.int32)


def dilate_labels(lab, r):
    """Grow every label outward by r px without merging neighbours."""
    if r <= 0:
        return lab
    from scipy import ndimage as ndi
    d, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    out = lab.copy()
    grow = (lab == 0) & (d <= r)
    out[grow] = lab[iy[grow], ix[grow]]
    return out


def occupancy(F, lab, pctile=BRIGHT_PCTILE):
    """Fraction of label objects whose mean Fluo-4 exceeds the field percentile."""
    from scipy import ndimage as ndi
    ids = np.unique(lab[lab > 0])
    if len(ids) == 0:
        return np.nan, 0
    thr = np.percentile(F, pctile)
    means = np.array(ndi.mean(F, lab, ids))
    return float((means > thr).mean()), len(ids)


# ------------------------------ QC panel ---------------------------------
def outline(m):
    m = (m > 0).astype(np.uint8)
    return (m - cv2.erode(m, np.ones((3, 3), np.uint8))) > 0


def n01(a, lo=1, hi=99.5):
    p1, p2 = np.percentile(a, (lo, hi))
    return np.clip((a - p1) / max(p2 - p1, 1e-6), 0, 1)


def index_icc(root):
    """{well: {channel_number: path}} from <well>_ch<N>.tif statics."""
    out = {}
    root = Path(root)
    if not root.is_dir():
        return out
    for p in sorted(root.glob("*.tif*")):
        w = re.search(r"(R\d+C\d+)", p.name, re.I)
        c = re.search(r"_ch(\d+)", p.name, re.I)
        if w and c:
            out.setdefault(w.group(1).upper(), {})[int(c.group(1))] = p
    return out


def to_shape(a, hw):
    return a if a.shape[:2] == tuple(hw) else cv2.resize(
        a.astype(np.float32), (hw[1], hw[0]), interpolation=cv2.INTER_AREA)


def mask_enrich(img, lab):
    """Mean intensity under the mask over the field mean."""
    m = lab > 0
    if not m.any():
        return np.nan
    mu = float(np.asarray(img).mean())
    return float(np.asarray(img)[m].mean() / mu) if mu else np.nan


def icc_panel(F, T, lab, title):
    """Green = live Fluo-4, red = TH ICC, white outline = aligned mask.
    A mask sitting on a yellow cell is on a cell that is both dye-loaded and
    TH-positive - the whole claim, visible."""
    g, r = n01(F), n01(T)
    rgb = np.dstack([r, g, np.zeros_like(g)])
    rgb[outline(lab)] = [1.0, 1.0, 1.0]
    return rgb, title


def panel(F, lab, title, ok):
    g = F - np.percentile(F, 1)
    g = np.clip(g / max(np.percentile(F, 99.5) - np.percentile(F, 1), 1e-6), 0, 1)
    rgb = np.dstack([g, g, g])
    rgb[outline(lab)] = [1.0, 0.25, 0.2] if ok else [1.0, 0.85, 0.1]
    return rgb, title


def contact_sheet(panels, out_png, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(panels)
    c = int(np.ceil(np.sqrt(n))); r = int(np.ceil(n / c))
    fig, ax = plt.subplots(r, c, figsize=(2.0 * c, 2.15 * r), dpi=150)
    ax = np.atleast_1d(ax).ravel()
    for a in ax:
        a.axis("off")
    for a, (img, t) in zip(ax, panels):
        a.imshow(img, interpolation="nearest")
        a.set_title(t, fontsize=5.2, pad=1.5)
    fig.suptitle(f"{label}  -  TH mask on Fluo-4 after per-well alignment "
                 f"(red = pass, yellow = fail)", fontsize=8, y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out_png, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  contact sheet -> {out_png}")


# -------------------------------- main -----------------------------------
def do_plate(cfg):
    PLATE_LABEL = cfg["label"]
    out = Path(OUTPUT_ROOT) / PLATE_LABEL
    (out / "aligned_masks").mkdir(parents=True, exist_ok=True)
    (out / "for_harmony").mkdir(parents=True, exist_ok=True)
    (out / "qa").mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}\n{PLATE_LABEL}\n{'='*70}", flush=True)
    print(f"  listing masks   {cfg['th_mask']}", flush=True)
    masks = index_masks(cfg["th_mask"], MASK_REQUIRE)
    print(f"    -> {len(masks)} well entries", flush=True)
    print(f"  listing statics {cfg['fluo']}", flush=True)
    fluos = index_flat(cfg["fluo"], FLUO_REQUIRE)
    print(f"    -> {len(fluos)} statics", flush=True)
    ths = index_flat(cfg["th_icc"]) if cfg.get("th_icc") else {}
    wells = sorted(set(masks) & set(fluos),
                   key=lambda k: [int(x) for x in re.findall(r"\d+", k)])
    print(f"  {len(wells)} matched well(s)\n", flush=True)
    if not wells:
        print("  !! no wells matched.")
        print(f"     mask keys : {sorted(masks)[:8]}")
        print(f"     fluo keys : {sorted(fluos)[:8]}")
        print(f"     check th_mask / fluo paths and MASK_REQUIRE / FLUO_REQUIRE")
        return []

    # A mask folder holding many more wells than the plate has is the signature
    # of a parent path: the union of several plates' well names.
    if len(masks) > 1.3 * len(fluos):
        print(f"  !! {len(masks)} mask wells but only {len(fluos)} statics.")
        print(f"  !! That usually means th_mask is a PARENT of several plates'")
        print(f"  !! folders, so each well gets an arbitrary plate's mask. Expect")
        print(f"  !! peak/rival near 1.00 and every well to FAIL. Point th_mask")
        print(f"  !! at the single folder for THIS plate and re-run.\n")

    recs, panels, calib = [], [], {r: [] for r in DILATE_TEST}
    iccs = index_icc(cfg["th_icc"]) if cfg.get("th_icc") else {}
    keep = {}          # well -> (Fluo-4, aligned labels) for the ICC stage
    for k in wells:
        F = to_float(read_tif_any(fluos[k]))
        mp = resolve_mask(masks[k])
        if mp is None:
            print(f"  {k}: no tif inside {Path(masks[k]).name}, skipped")
            continue
        L = read_labels(mp)
        if L.shape != F.shape:
            L = cv2.resize(L, (F.shape[1], F.shape[0]),
                           interpolation=cv2.INTER_NEAREST)
        res = align(F, L > 0)
        if res is None:
            print(f"  {k}: empty mask, skipped")
            continue
        La = shift_labels(L, res["dx"], res["dy"])

        # independent check that the warp sign matches the search
        chk = float(F[La > 0].mean() / F.mean()) if (La > 0).any() else 0.0
        if abs(chk - res["enrich"]) > 0.15 * max(res["enrich"], 1e-6):
            print(f"  {k}: WARNING warp/search mismatch "
                  f"({chk:.3f} vs {res['enrich']:.3f}) - check sign convention")

        occ, nobj = occupancy(F, La)
        ok = (res["peak_ratio"] >= MIN_PEAK_RATIO and
              (True if SUBTRACT_BACKGROUND else res["enrich"] >= MIN_ENRICH))
        recs.append(dict(plate=PLATE_LABEL, well=k,
                         dx=res["dx"], dy=res["dy"],
                         shift_px=round(float(np.hypot(res["dx"], res["dy"])), 1),
                         shift_um=round(float(np.hypot(res["dx"], res["dy"])) * PIXEL_UM, 1),
                         enrich=res["enrich"], enrich_unaligned=res["enrich0"],
                         peak_ratio=res["peak_ratio"],
                         n_objects=nobj, occupancy=round(occ, 3) if occ == occ else np.nan,
                         verdict="PASS" if ok else "FAIL"))
        print(f"  {k}: shift ({res['dx']:+4d},{res['dy']:+4d}) = "
              f"{recs[-1]['shift_px']:5.1f} px | enrich {res['enrich']:.2f} "
              f"(was {res['enrich0']:.2f}) | peak/rival {res['peak_ratio']:.2f} | "
              f"occ {occ:.2f} | {recs[-1]['verdict']}")

        if CALIBRATE_DILATION:
            for r in DILATE_TEST:
                o, _ = occupancy(F, dilate_labels(La, r))
                calib[r].append(o)

        Lw = dilate_labels(La, DILATE_PX)
        if WRITE_FLAT_COPY:
            tiff.imwrite(str(out / "aligned_masks" / f"{k}_TH_labels_aligned.tif"),
                         Lw.astype(np.uint16), compression=None)
        tree = out / "for_harmony" / f"{k}_F1T1P0_TH" / "MAP" / "DAPI"
        tree.mkdir(parents=True, exist_ok=True)
        tiff.imwrite(str(tree / "SOX+ Selected_Cell.tiff"),
                     Lw.astype(np.uint16), compression=None)

        img, t = panel(F, Lw, f"{k}  {recs[-1]['shift_px']:.0f}px  "
                              f"e{res['enrich']:.2f}  r{res['peak_ratio']:.2f}", ok)
        panels.append((img, t))
        if iccs:
            keep[k] = (F.astype(np.float32), Lw.astype(np.uint16))

    import pandas as pd

    # ---------------- ICC stage: which channel is TH, and does the mask sit
    # ---------------- on TH-positive cells?
    th_ch, ch_scores = None, {}
    if iccs and keep:
        print("\n  scoring ICC channels (the mask came from TH, so TH is the "
              "channel it lights up) ...", flush=True)
        probe = [k for k in wells if k in iccs and k in keep][:ICC_DETECT_N]
        for k in probe:
            F, L = keep[k]
            for ch, p in sorted(iccs[k].items()):
                try:
                    T = to_shape(to_float(read_tif_any(p)), F.shape)
                except Exception:
                    continue
                ch_scores.setdefault(ch, []).append(mask_enrich(T, L))
        if ch_scores:
            med = {c: float(np.nanmedian(v)) for c, v in ch_scores.items()}
            th_ch = ICC_CHANNEL or max(med, key=med.get)
            print("    channel   mask enrichment in that channel")
            for c in sorted(med):
                tag = "   <-- TH" if c == th_ch else ""
                print(f"      ch{c}        {med[c]:.2f}{tag}")
            if ICC_CHANNEL:
                print(f"    (ch{th_ch} forced by ICC_CHANNEL)")

    if th_ch:
        ip = []
        for k in wells:
            if k not in keep or k not in iccs or th_ch not in iccs[k]:
                continue
            F, L = keep[k]
            try:
                T = to_shape(to_float(read_tif_any(iccs[k][th_ch])), F.shape)
            except Exception:
                continue
            e = mask_enrich(T, L)
            for r in recs:
                if r["well"] == k:
                    r[f"icc_ch{th_ch}_enrich"] = round(e, 3)
            ip.append(icc_panel(F, T, L, f"{k}   TH enrich {e:.2f}"))
        if ip:
            contact_sheet(ip, str(out / "qa" / f"contact_sheet_TH_ch{th_ch}.png"),
                          f"{PLATE_LABEL}  ch{th_ch}=TH (red) + Fluo-4 (green), "
                          f"white = aligned mask")

    D = pd.DataFrame(recs)
    D.to_csv(out / "alignment_qc.csv", index=False)
    contact_sheet(panels, str(out / "qa" / "contact_sheet.png"), PLATE_LABEL)
    if th_ch and f"icc_ch{th_ch}_enrich" in D:
        v = D[f"icc_ch{th_ch}_enrich"].dropna()
        if len(v):
            print(f"\n  TH ICC enrichment under the aligned mask: median "
                  f"{v.median():.2f} (n={len(v)} wells)")
            print("  Above ~1.5 means the mask really is sitting on TH-positive")
            print("  cells, not merely on cells. That is the claim your gate needs.")

    print("\n" + "=" * 70)
    print(f"{PLATE_LABEL}   {len(D)} wells")
    print("=" * 70)
    print(f"  shift        median {D.shift_px.median():.1f} px "
          f"({D.shift_um.median():.1f} um), range {D.shift_px.min():.0f}-{D.shift_px.max():.0f}")
    print(f"  enrichment   median {D.enrich.median():.2f} "
          f"(unaligned {D.enrich_unaligned.median():.2f})")
    print(f"  peak/rival   median {D.peak_ratio.median():.2f}")
    print(f"  occupancy    median {D.occupancy.median():.2f}")
    print(f"  PASS {int((D.verdict=='PASS').sum())} / {len(D)}   "
          f"FAIL: {', '.join(D.loc[D.verdict=='FAIL','well']) or 'none'}")
    print("\n  Per-well spread of shifts is the number that justifies abandoning a")
    print("  plate-level translation. If min and max differ by more than a soma")
    print("  width (~14 px), one shift for the whole plate cannot be correct.")

    if CALIBRATE_DILATION:
        print("\n" + "=" * 70)
        print("DILATION CALIBRATION   -> pick DILATE_PX, keep Harmony at 0%")
        print("=" * 70)
        print("  growth   mean occupancy   gain over previous")
        prev = None
        for r in DILATE_TEST:
            v = float(np.nanmean(calib[r]))
            g = "" if prev is None else f"{v-prev:+.3f}"
            print(f"  {r:3d} px      {v:.3f}            {g}")
            prev = v
        print("\n  Take the smallest growth where the gain drops below ~0.02, put it")
        print("  in DILATE_PX, re-run. The exported mask then carries the growth and")
        print("  Harmony's Outer Border stays at 0%.")

    if D.shift_px.median() < 3 and D.shift_px.max() < 8:
        print("\n  !! every shift is tiny. Either these masks were already aligned")
        print("     (are you pointing th_mask at for_harmony_full?) or the plate")
        print("     genuinely did not move. Check one well by eye before trusting it.")

    with open(out / "summary_v2.json", "w") as f:
        json.dump(dict(plate=PLATE_LABEL, dilate_px=DILATE_PX,
                       min_enrich=MIN_ENRICH, min_peak_ratio=MIN_PEAK_RATIO,
                       max_shift_px=MAX_SHIFT_PX, fluo_require=FLUO_REQUIRE,
                       wells=recs), f, indent=2)
    print(f"\n  alignment_qc.csv, summary_v2.json, for_harmony/  -> {out}")
    return recs


def build_plates():
    """Turn JOBS into the list of dicts do_plate expects, checking every path
    BEFORE any work. A wrong base path or a stray '...' from a pasted message
    otherwise costs a full run per plate."""
    if not JOBS:
        raise SystemExit("JOBS is empty. Add at least one line.")
    out, bad = [], []
    for label, (mask, fluo) in JOBS.items():
        cands = list(mask) if isinstance(mask, (list, tuple)) else [mask]
        fp = Path(STATICS_ROOT) / fluo
        ip = ICC.get(label, "")
        for tag, p in (("statics", fp),) + ((("icc", Path(ip)),) if ip else ()):
            sp = str(p)
            if "..." in sp:
                bad.append(f"{label} {tag}: contains '...' - that is an ellipsis "
                           f"from a chat message, not a path")
            elif not Path(sp).is_dir():
                bad.append(f"{label} {tag}: does not exist\n        {sp}")
        for c in cands:
            mp = Path(c) if (str(c).startswith(("/", "\\")) or ":" in str(c)[:3]) \
                 else Path(MASK_ROOT) / c
            lbl = label if len(cands) == 1 else f"{label}__{Path(c).name[:8]}"
            sp = str(mp)
            if "..." in sp:
                bad.append(f"{lbl} mask: contains '...' - that is an ellipsis "
                           f"from a chat message, not a path")
            elif not mp.is_dir():
                bad.append(f"{lbl} mask: does not exist\n        {sp}")
            out.append(dict(label=lbl, th_mask=str(mp), fluo=str(fp), th_icc=str(ip)))
    print("checking paths ...")
    if bad:
        print("\n  !! " + "\n  !! ".join(bad))
        raise SystemExit(f"\n{len(bad)} bad path(s). Nothing was run.\n"
                         f"MASK_ROOT    = {MASK_ROOT}\n"
                         f"STATICS_ROOT = {STATICS_ROOT}")
    for j in out:
        print(f"  {j['label']:<11} mask {Path(j['th_mask']).name[:40]:<42}"
              f"statics {Path(j['fluo']).name}")
    print(f"  all {len(out)} plate(s) OK\n")
    return out


def main():
    import pandas as pd
    PLATES = build_plates()
    allrecs = []
    for cfg in PLATES:
        try:
            allrecs += do_plate(cfg)
        except Exception as e:
            print(f"  !! {cfg.get('label','?')} failed: {type(e).__name__}: {str(e)[:120]}")
    if not allrecs:
        raise SystemExit("\nNothing aligned.")
    A = pd.DataFrame(allrecs)
    Path(OUTPUT_ROOT).mkdir(parents=True, exist_ok=True)
    A.to_csv(Path(OUTPUT_ROOT) / "alignment_qc_ALL.csv", index=False)

    print("\n" + "=" * 70)
    print("ACROSS PLATES   -> this table is your exclusions table")
    print("=" * 70)
    g = (A.groupby("plate")
          .agg(wells=("well", "size"),
               pass_=("verdict", lambda v: int((v == "PASS").sum())),
               shift_med=("shift_px", "median"), shift_min=("shift_px", "min"),
               shift_max=("shift_px", "max"), enrich_med=("enrich", "median"),
               ratio_med=("peak_ratio", "median"), occ_med=("occupancy", "median"))
          .reset_index())
    g["spread_px"] = (g.shift_max - g.shift_min).round(1)
    print(g.round(2).to_string(index=False))
    print("\n  spread_px is the range of per-well shifts on that plate. Anything")
    print("  above ~14 px (one soma) means a single plate-level translation could")
    print("  never have been right, which is the sentence your methods needs.")

    # where several candidate masks were tried for one plate, declare a winner
    A["base"] = A.plate.str.split("__").str[0]
    multi = [b for b, g in A.groupby("base") if g.plate.nunique() > 1]
    if multi:
        print("\n" + "=" * 70)
        print("CANDIDATE MASKS TRIED  -> winner by median peak/rival")
        print("=" * 70)
        for b in multi:
            g = (A[A.base == b].groupby("plate")
                 .agg(wells=("well", "size"), ratio=("peak_ratio", "median"),
                      enrich=("enrich", "median"), occ=("occupancy", "median"),
                      pass_=("verdict", lambda v: int((v == "PASS").sum())))
                 .sort_values("ratio", ascending=False))
            print(f"\n  {b}")
            for pl, r in g.iterrows():
                cand = pl.split("__")[-1]
                v = ("<-- WINNER" if r.ratio >= 1.30 else
                     "borderline" if r.ratio >= 1.15 else "no")
                print(f"    {cand:<10} {int(r.wells):>3} wells  peak/rival "
                      f"{r.ratio:5.2f}  enrich {r.enrich:5.2f}  "
                      f"occ {r.occ:4.2f}  PASS {int(r.pass_):>2}  {v}")
            if g.ratio.max() < 1.15:
                print("    none convincing - no mask in the candidate set belongs")
                print("    to this plate, or the plate is too sparse to register.")
    fails = A[A.verdict == "FAIL"]
    if len(fails):
        print(f"\n  FAIL wells ({len(fails)} of {len(A)}):")
        for p, sub in fails.groupby("plate"):
            print(f"    {p}: {', '.join(sub.well)}")
    print(f"\n  alignment_qc_ALL.csv -> {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
