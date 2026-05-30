"""Long-read polyploid phasing via anchor extension (paraphase-inspired).

Used in place of `whatshap polyphase` when ploidy > 2 and reads are long.
Algorithm:
  1. Build the (read × variant) allele matrix the same way as the short-read
     path — for long reads, pileup-truncate is still efficient because each
     read deposits its allele at every variant column it spans.
  2. Pick anchor variants — biallelic phaseable SNVs whose GT multiset is the
     most balanced (a·b for [a copies of 0, b copies of alt]). We want about
     ⌈log2(K)⌉+1 anchors so K haplotypes get distinct barcodes.
  3. Each read with all anchors observed contributes a barcode = its alleles
     at the anchor columns. The K most-populous distinct barcodes become
     cluster seeds.
  4. Every read is assigned to the closest seed (mismatch rate on jointly
     observed sites). For long reads, one pass is almost always correct.
  5. Single k-modes refinement pass — replaces anchor-only seeds with
     per-site majority centers, then reassigns. Catches reads where the
     anchor was noisy.
  6. Same downstream pipeline as polyphase_short: cluster→GT-multiset
     assignment, phase blocks (with a relaxed connecting-reads threshold —
     a single long read connects many variants), emit VCF with phased GT + PS.

Falls back to full iterative k-modes if anchor seeding can't produce K
distinct complete-barcode clusters (low coverage, no good biallelic
anchors, etc.).
"""

import os
from collections import Counter
from math import ceil, log2
from typing import Any, Dict, List, Optional

import numpy as np

from ..logging import get_logger
from .polyphase_short import (
    _assign_clusters_to_gt,
    _build_matrix,
    _connectivity_within_clusters,
    _emit_vcf,
    _kmodes,
    _load_variants,
    _make_phase_blocks,
    _pairwise_to_center,
    _passthrough,
    _per_site_cluster_counts,
    _symlink_through,
    _update_center,
)

# Long-read tunables
MIN_CONNECTING_READS_LONG = 2     # one long read can connect many variants
N_REFINEMENT_ITERS = 1            # one k-modes pass after anchor seeding

logger = get_logger("polyphase_long")


# -- anchor selection -------------------------------------------------------
def _balanced_score(gt: List[int]) -> int:
    """For biallelic SNV: a·b where a, b are counts of allele 0 and non-0.

    Max at ⌊K/2⌋·⌈K/2⌉ (balanced multiset). Zero for homozygous.
    """
    if not gt:
        return 0
    a = sum(1 for x in gt if x == 0)
    b = len(gt) - a
    return a * b


def _pick_anchors(
    variants: List[Dict[str, Any]],
    phaseable_idx: List[int],
    K: int,
) -> List[int]:
    """Return matrix column indices of top biallelic anchors by balance score.

    Target ⌈log2(K)⌉+1 anchors so K haplotypes can be distinguished by barcode.
    """
    target = max(2, ceil(log2(max(K, 2))) + 1)
    candidates: List = []
    for c, vi in enumerate(phaseable_idx):
        v = variants[vi]
        if len(v["alt"]) != 1:
            continue
        score = _balanced_score(v["gt"])
        if score > 0:
            candidates.append((-score, c))   # negative so sort ascending = best first
    candidates.sort()
    return [c for _, c in candidates[:target]]


# -- anchor-based seeding ---------------------------------------------------
def _anchor_seeds(
    matrix: np.ndarray, anchor_cols: List[int], K: int
) -> Optional[np.ndarray]:
    """Pick K seed read indices from the K most-populous distinct anchor barcodes.

    Returns None when there aren't K distinct complete barcodes available —
    the caller falls back to iterative k-modes in that case.
    """
    n_reads = matrix.shape[0]
    if not anchor_cols or n_reads < K:
        return None

    bc = matrix[:, anchor_cols]
    complete_mask = (bc >= 0).all(axis=1)
    if int(complete_mask.sum()) < K:
        return None
    complete_rows = np.where(complete_mask)[0]
    complete = bc[complete_mask]

    counter = Counter(tuple(int(x) for x in row) for row in complete)
    if len(counter) < K:
        return None
    top_barcodes = [bc_t for bc_t, _ in counter.most_common(K)]

    seed_rows: List[int] = []
    used = set()
    for bc_t in top_barcodes:
        for ri, row in zip(complete_rows, complete):
            ri_int = int(ri)
            if ri_int in used:
                continue
            if tuple(int(x) for x in row) == bc_t:
                seed_rows.append(ri_int)
                used.add(ri_int)
                break
    if len(seed_rows) != K:
        return None
    return np.array(seed_rows, dtype=np.int64)


def _assign_to_centers(matrix: np.ndarray, centers: np.ndarray) -> np.ndarray:
    """Assign each read to the closest center by jointly-observed mismatch rate."""
    K = centers.shape[0]
    dists = np.stack(
        [_pairwise_to_center(matrix, centers[k]) for k in range(K)], axis=0
    )
    return np.argmin(dists, axis=0).astype(np.int32)


def _refine_one_pass(
    matrix: np.ndarray, labels: np.ndarray, K: int, max_allele: int
) -> np.ndarray:
    """Recompute per-cluster majority centers, then reassign reads."""
    n_cols = matrix.shape[1]
    centers = np.full((K, n_cols), -1, dtype=np.int8)
    for k in range(K):
        members = matrix[labels == k]
        centers[k] = _update_center(members, n_cols, max_allele)
    return _assign_to_centers(matrix, centers)


# -- main driver ------------------------------------------------------------
def polyphase_long_reads(
    bam: Any,
    in_vcf: str,
    out_vcf: str,
    ploidy: int,
) -> None:
    """Polyploid phaser for long reads. Drop-in for `whatshap polyphase`."""
    if bam.use_existing and os.path.exists(out_vcf) and os.path.getsize(out_vcf):
        return
    if ploidy < 2:
        _symlink_through(in_vcf, out_vcf, bam.sentieon)
        return

    variants = _load_variants(in_vcf)
    if not any(v["phaseable"] for v in variants):
        _passthrough(in_vcf, out_vcf, bam.sentieon)
        return

    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    matrix, _read_ids, phaseable_idx = _build_matrix(bam_path, bam.ref, variants)
    logger.debug(
        f"polyphase_long: ploidy={ploidy}, phaseable_sites={len(phaseable_idx)}, "
        f"reads_in_matrix={matrix.shape[0]}"
    )

    n_phaseable = len(phaseable_idx)
    if matrix.shape[0] < ploidy or n_phaseable < 2:
        _passthrough(in_vcf, out_vcf, bam.sentieon)
        return

    max_allele = int(max(len(variants[i]["alt"]) for i in phaseable_idx))

    # Anchor-based seeding ----------------------------------------------------
    anchor_cols = _pick_anchors(variants, phaseable_idx, ploidy)
    seed_rows = _anchor_seeds(matrix, anchor_cols, ploidy) if anchor_cols else None

    if seed_rows is None:
        logger.debug("anchor seeding insufficient — falling back to k-modes")
        labels, _centers = _kmodes(matrix, ploidy, max_allele)
    else:
        seed_centers = matrix[seed_rows]
        labels = _assign_to_centers(matrix, seed_centers)
        for _ in range(N_REFINEMENT_ITERS):
            new_labels = _refine_one_pass(matrix, labels, ploidy, max_allele)
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels

    # Downstream is shared with polyphase_short -------------------------------
    counts = _per_site_cluster_counts(matrix, labels, ploidy, max_allele)

    assigned_per_site: List[List[int]] = []
    confident_per_site = np.zeros(n_phaseable, dtype=bool)
    for c, vi in enumerate(phaseable_idx):
        gt_ms = variants[vi]["gt"]
        if len(gt_ms) != ploidy:
            assigned_per_site.append([])
            confident_per_site[c] = False
            continue
        assigned, conf = _assign_clusters_to_gt(counts[:, c, :], list(gt_ms), ploidy)
        assigned_per_site.append(assigned)
        confident_per_site[c] = conf

    # Block boundaries from PER-CLUSTER connectivity (min over clusters of reads
    # spanning each adjacent pair), same as the short-read path. Breaking where a
    # haplotype is weakly linked isolates phase-switch points so the downstream
    # decomposition (assign(), per-PS-block) can't mis-assign an over-merged
    # switched run. Long reads keep blocks long anyway (per-cluster connectivity
    # stays high), so this stays well above whatshap's fragmentation while staying
    # decomposition-safe.
    connectivity = _connectivity_within_clusters(matrix, labels, ploidy)
    variant_positions = [variants[vi]["pos"] for vi in phaseable_idx]
    ps_per_col = _make_phase_blocks(
        confident_per_site,
        connectivity,
        variant_positions,
        min_conn=MIN_CONNECTING_READS_LONG,
    )

    _emit_vcf(in_vcf, out_vcf, variants, phaseable_idx, assigned_per_site, ps_per_col, bam.sentieon)
