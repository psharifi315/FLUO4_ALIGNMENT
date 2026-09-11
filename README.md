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


if __name__ == "__main__":
    main()
