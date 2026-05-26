"""Short-read polyploid phasing via read clustering (k-modes).

Used in place of `whatshap polyphase` when ploidy > 2 and reads are short.
Algorithm:
  1. Pull every read overlapping multiple phaseable (heterozygous SNV) sites,
     build a (read × variant) allele matrix with -1 for unobserved.
  2. Run k-modes with K = ploidy: distance is mismatch rate over jointly
     observed sites; cluster centers are per-site majority alleles. Use
     farthest-first init + a few random restarts; pick the lowest-cost run.
  3. Per variant, solve a small assignment problem to map clusters onto the
     original GT multiset so the phased GT is a permutation of the unphased one.
  4. Break phase sets where adjacent variants are not connected by enough
     within-cluster reads, or where a cluster has insufficient coverage.
  5. Stream the input VCF, attach phased GT and PS to each phaseable site.
"""

import os
import re
import subprocess
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pysam
import vcflib

from ..logging import get_logger

# -- tunables ---------------------------------------------------------------
MIN_BASE_QUAL = 13            # min base quality at a variant site
MIN_MAP_QUAL = 20             # min mapping quality
MAX_KMODES_ITER = 25
N_RESTARTS = 3
MIN_VARIANTS_PER_READ = 2     # reads observing fewer sites carry no phase info
MIN_READS_PER_CLUSTER_AT_SITE = 2
MIN_CONNECTING_READS = 2      # min reads spanning a variant pair within a cluster
RNG_SEED = 0x5EDD0

logger = get_logger("polyphase_short")


# -- VCF I/O ----------------------------------------------------------------
def _parse_gt(gt_str: str) -> Optional[List[int]]:
    """Parse '0/1/0/1' or '0|1|0|1' into a list of ints. None if any allele is '.'."""
    if not gt_str:
        return None
    parts = re.split(r"[|/]", gt_str)
    out = []
    for p in parts:
        if p == "." or p == "":
            return None
        try:
            out.append(int(p))
        except ValueError:
            return None
    return out


def _is_snv(ref: str, alts: List[str]) -> bool:
    """True iff ref and every alt is a single base (no indel/MNV)."""
    if len(ref) != 1:
        return False
    return all(len(a) == 1 for a in alts)


def _load_variants(in_vcf: str) -> List[Dict[str, Any]]:
    """Return one dict per record. Phaseable iff SNV and GT has >1 distinct allele.

    Records are returned in file order. Each dict:
      chrom, pos (0-based), ref, alt (list), gt (list[int]), phaseable (bool)
    """
    out: List[Dict[str, Any]] = []
    vcf = vcflib.VCF(in_vcf)
    try:
        for v in vcf:
            gt = _parse_gt(v.samples[0].get("GT", "")) if v.samples else None
            alts = list(v.alt) if v.alt else []
            phaseable = (
                gt is not None
                and len(set(gt)) > 1
                and _is_snv(v.ref, alts)
            )
            out.append({
                "chrom": v.chrom,
                "pos": v.pos,           # 0-based
                "ref": v.ref,
                "alt": alts,
                "gt": gt,
                "phaseable": phaseable,
            })
    finally:
        vcf.close()
    return out


# -- matrix building --------------------------------------------------------
def _open_bam(bam_path: str, ref: str) -> pysam.AlignmentFile:
    if bam_path.endswith(".cram"):
        return pysam.AlignmentFile(bam_path, "rc", reference_filename=ref)
    return pysam.AlignmentFile(bam_path, "rb")


def _match_snv_allele(base: str, ref_base: str, alt_bases: List[str]) -> Optional[int]:
    base = base.upper()
    if base == ref_base.upper():
        return 0
    for i, a in enumerate(alt_bases, start=1):
        if base == a.upper():
            return i
    return None


def _build_matrix(
    bam_path: str,
    ref: str,
    variants: List[Dict[str, Any]],
) -> Tuple[np.ndarray, List[str], List[int]]:
    """Pileup the union region once; record per-(read, phaseable-variant) allele.

    Returns:
      matrix: (n_reads, n_phaseable) int8, with -1 for unobserved
      read_ids: query_name for each row (for debugging)
      phaseable_idx: index into the full `variants` list for each matrix column
    """
    phaseable_idx = [i for i, v in enumerate(variants) if v["phaseable"]]
    if len(phaseable_idx) < 2:
        return np.zeros((0, len(phaseable_idx)), dtype=np.int8), [], phaseable_idx

    # phaseable variants grouped by chromosome (preserving order)
    chroms: Dict[str, List[int]] = {}
    for col, vi in enumerate(phaseable_idx):
        chroms.setdefault(variants[vi]["chrom"], []).append(col)

    pos_to_col: Dict[Tuple[str, int], int] = {
        (variants[phaseable_idx[c]]["chrom"], variants[phaseable_idx[c]]["pos"]): c
        for c in range(len(phaseable_idx))
    }

    read_alleles: Dict[str, Dict[int, int]] = {}

    bamh = _open_bam(bam_path, ref)
    try:
        for chrom, cols in chroms.items():
            positions = [variants[phaseable_idx[c]]["pos"] for c in cols]
            start_1b = min(positions) + 1
            end_1b = max(positions) + 1
            region_str = f"{chrom}:{start_1b}-{end_1b}"
            try:
                pileup = bamh.pileup(
                    region=region_str,
                    truncate=True,
                    min_base_quality=MIN_BASE_QUAL,
                    min_mapping_quality=MIN_MAP_QUAL,
                    stepper="samtools",
                )
            except ValueError:
                # No reads / bad region
                continue
            for col_obj in pileup:
                key = (chrom, col_obj.reference_pos)
                if key not in pos_to_col:
                    continue
                mat_col = pos_to_col[key]
                vi = phaseable_idx[mat_col]
                ref_b = variants[vi]["ref"]
                alts = variants[vi]["alt"]
                for pread in col_obj.pileups:
                    if pread.is_del or pread.is_refskip:
                        continue
                    if pread.query_position is None:
                        continue
                    read = pread.alignment
                    if (
                        read.is_unmapped
                        or read.is_secondary
                        or read.is_supplementary
                        or read.is_duplicate
                    ):
                        continue
                    seq = read.query_sequence
                    if seq is None:
                        continue
                    base = seq[pread.query_position]
                    allele = _match_snv_allele(base, ref_b, alts)
                    if allele is None:
                        continue
                    rid = read.query_name
                    read_alleles.setdefault(rid, {})[mat_col] = allele
    finally:
        bamh.close()

    read_ids = [r for r, sites in read_alleles.items() if len(sites) >= MIN_VARIANTS_PER_READ]
    n_reads = len(read_ids)
    n_cols = len(phaseable_idx)
    matrix = np.full((n_reads, n_cols), -1, dtype=np.int8)
    for r, rid in enumerate(read_ids):
        for c, a in read_alleles[rid].items():
            matrix[r, c] = a
    return matrix, read_ids, phaseable_idx


# -- k-modes clustering -----------------------------------------------------
def _pairwise_to_center(reads: np.ndarray, center: np.ndarray) -> np.ndarray:
    """For each row in `reads`, mismatch rate vs `center` on jointly observed sites.

    Reads with no overlap get distance = 1.0 (worst).
    """
    obs = (reads >= 0) & (center >= 0)  # broadcasts
    n_obs = obs.sum(axis=1)
    n_diff = ((reads != center) & obs).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(n_obs > 0, n_diff / np.maximum(n_obs, 1), 1.0)
    return d


def _update_center(members: np.ndarray, n_cols: int, max_allele: int) -> np.ndarray:
    """Per-site majority allele among `members` (each row is a read).

    Sites with no observation in any member become -1.
    """
    center = np.full(n_cols, -1, dtype=np.int8)
    if members.shape[0] == 0:
        return center
    for c in range(n_cols):
        col = members[:, c]
        obs = col[col >= 0]
        if obs.size == 0:
            continue
        counts = np.bincount(obs.astype(np.int64), minlength=max_allele + 1)
        center[c] = int(np.argmax(counts))
    return center


def _farthest_first_init(matrix: np.ndarray, K: int, rng: np.random.Generator) -> np.ndarray:
    """Pick K rows as initial centers by farthest-first traversal.

    Seed: highest-coverage read (most observed sites). Ties broken randomly.
    """
    coverage = (matrix >= 0).sum(axis=1)
    n_reads = matrix.shape[0]
    # Pick seed among top-coverage reads (randomized tie-break)
    top = np.flatnonzero(coverage == coverage.max())
    seed = int(rng.choice(top))
    chosen = [seed]
    dists = _pairwise_to_center(matrix, matrix[seed])
    for _ in range(1, K):
        # Pick the read with max distance to any chosen center (largest min-distance)
        nxt = int(np.argmax(dists))
        if nxt in chosen:
            # Fall back to a random uncovered pick
            remaining = [i for i in range(n_reads) if i not in chosen]
            if not remaining:
                break
            nxt = int(rng.choice(remaining))
        chosen.append(nxt)
        new_d = _pairwise_to_center(matrix, matrix[nxt])
        dists = np.minimum(dists, new_d)
    while len(chosen) < K:
        # If we ran out (very few reads), duplicate seeds
        chosen.append(chosen[-1])
    return np.array(chosen[:K], dtype=np.int64)


def _kmodes_once(
    matrix: np.ndarray,
    K: int,
    max_allele: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """One k-modes run. Returns (labels, centers, total_cost)."""
    n_reads, n_cols = matrix.shape
    seed_idx = _farthest_first_init(matrix, K, rng)
    centers = matrix[seed_idx].copy()
    labels = np.zeros(n_reads, dtype=np.int32)

    for _ in range(MAX_KMODES_ITER):
        # Assignment: for each read, distance to each center → argmin
        # Stack distances: shape (K, n_reads)
        dists = np.stack([_pairwise_to_center(matrix, centers[k]) for k in range(K)], axis=0)
        new_labels = np.argmin(dists, axis=0).astype(np.int32)
        if np.array_equal(new_labels, labels):
            labels = new_labels
            break
        labels = new_labels
        # Update each center
        for k in range(K):
            members = matrix[labels == k]
            centers[k] = _update_center(members, n_cols, max_allele)

    # Final cost: sum of read→assigned-center distances
    final_d = np.stack([_pairwise_to_center(matrix, centers[k]) for k in range(K)], axis=0)
    cost = float(final_d[labels, np.arange(n_reads)].sum())
    return labels, centers, cost


def _kmodes(matrix: np.ndarray, K: int, max_allele: int) -> Tuple[np.ndarray, np.ndarray]:
    """Run k-modes with N_RESTARTS starts; keep the lowest-cost result.

    If n_reads < K we collapse to as many clusters as we have reads.
    """
    n_reads = matrix.shape[0]
    if n_reads == 0:
        return np.zeros(0, dtype=np.int32), np.full((K, matrix.shape[1]), -1, dtype=np.int8)
    if n_reads < K:
        labels = np.arange(n_reads, dtype=np.int32)
        centers = np.full((K, matrix.shape[1]), -1, dtype=np.int8)
        for k in range(n_reads):
            centers[k] = matrix[k]
        return labels, centers

    best = None
    for r in range(N_RESTARTS):
        rng = np.random.default_rng(RNG_SEED + r)
        labels, centers, cost = _kmodes_once(matrix, K, max_allele, rng)
        if best is None or cost < best[2]:
            best = (labels, centers, cost)
    assert best is not None
    return best[0], best[1]


# -- cluster → GT-multiset assignment ---------------------------------------
def _per_site_cluster_counts(
    matrix: np.ndarray, labels: np.ndarray, K: int, max_allele: int
) -> np.ndarray:
    """counts[k, c, a] = # reads in cluster k observing allele a at column c."""
    n_cols = matrix.shape[1]
    counts = np.zeros((K, n_cols, max_allele + 1), dtype=np.int32)
    for k in range(K):
        members = matrix[labels == k]
        for c in range(n_cols):
            col = members[:, c]
            obs = col[col >= 0]
            if obs.size == 0:
                continue
            bc = np.bincount(obs.astype(np.int64), minlength=max_allele + 1)
            counts[k, c] = bc
    return counts


def _assign_clusters_to_gt(
    counts_at_site: np.ndarray, gt_multiset: List[int], K: int
) -> Tuple[List[int], bool]:
    """Map K clusters to the K alleles in the multiset to maximize total support.

    counts_at_site: shape (K, max_allele+1) — vote count per cluster per allele.
    gt_multiset: list of K alleles (with repeats), e.g. [0,0,1,1].
    Returns (assigned_alleles_per_cluster, confident_bool).
    `confident` is False if any cluster has zero observations at this site.
    """
    # Build a profit matrix M[k, slot] = count of cluster k for gt_multiset[slot]
    M = np.array(
        [[counts_at_site[k, gt_multiset[s]] for s in range(K)] for k in range(K)],
        dtype=np.int64,
    )
    # Maximize total assignment value via scipy if available, else brute force.
    # For K ≤ 10, brute force over permutations of slots is fine (10! = 3.6M).
    # We negate for "minimize" semantics.
    try:
        from scipy.optimize import linear_sum_assignment
        row_ind, col_ind = linear_sum_assignment(-M)
        slot_for_cluster = list(col_ind)
    except Exception:  # pragma: no cover
        from itertools import permutations
        best_val = -1
        slot_for_cluster = list(range(K))
        for perm in permutations(range(K)):
            v = sum(M[k, perm[k]] for k in range(K))
            if v > best_val:
                best_val = v
                slot_for_cluster = list(perm)
    assigned = [gt_multiset[slot] for slot in slot_for_cluster]
    # Confidence: every cluster had at least one observation
    cluster_total = counts_at_site.sum(axis=1)
    confident = bool(np.all(cluster_total >= MIN_READS_PER_CLUSTER_AT_SITE))
    return assigned, confident


# -- phase block boundaries -------------------------------------------------
def _connectivity_within_clusters(
    matrix: np.ndarray, labels: np.ndarray, K: int
) -> np.ndarray:
    """For each adjacent column pair, min over clusters of #reads observing both.

    Returns array shape (n_cols - 1,). Used to decide PS breaks.
    """
    n_cols = matrix.shape[1]
    if n_cols < 2:
        return np.zeros(0, dtype=np.int32)
    conn = np.full(n_cols - 1, 1 << 30, dtype=np.int32)
    for k in range(K):
        members = matrix[labels == k]
        if members.shape[0] == 0:
            conn[:] = 0
            continue
        obs = members >= 0
        # spans[c] = reads with obs at column c AND column c+1
        spans = (obs[:, :-1] & obs[:, 1:]).sum(axis=0)
        conn = np.minimum(conn, spans.astype(np.int32))
    return conn


def _make_phase_blocks(
    confident_per_site: np.ndarray,
    connectivity: np.ndarray,
    variant_positions: List[int],
    min_conn: int = MIN_CONNECTING_READS,
) -> List[int]:
    """Assign a phase-set ID to each phaseable column, or 0 if unphased.

    PS ID = position-1-based of first phaseable variant in the block (whatshap convention).
    A block is a maximal run of confident sites with sufficient connectivity to the next.
    """
    n = len(confident_per_site)
    ps = [0] * n
    i = 0
    while i < n:
        if not confident_per_site[i]:
            i += 1
            continue
        j = i
        while (
            j + 1 < n
            and confident_per_site[j + 1]
            and connectivity[j] >= min_conn
        ):
            j += 1
        if j > i:  # block must contain ≥2 sites to be useful
            block_id = variant_positions[i] + 1
            for k in range(i, j + 1):
                ps[k] = block_id
        i = j + 1
    return ps


# -- main driver ------------------------------------------------------------
def polyphase_short_reads(
    bam: Any,
    in_vcf: str,
    out_vcf: str,
    ploidy: int,
) -> None:
    """Polyploid phaser for short reads. Drop-in for `whatshap polyphase`.

    Output: VCF with phased GT (pipe-separated) and PS tag on phased records.
    Records that cannot be phased (indels, isolated het sites, low coverage)
    are emitted with their original unphased GT.
    """
    if bam.use_existing and os.path.exists(out_vcf) and os.path.getsize(out_vcf):
        return
    if ploidy < 2:
        # Caller should have handled this; emit a symlink for safety.
        _symlink_through(in_vcf, out_vcf, bam.sentieon)
        return

    variants = _load_variants(in_vcf)
    if not any(v["phaseable"] for v in variants):
        _passthrough(in_vcf, out_vcf, bam.sentieon)
        return

    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    matrix, read_ids, phaseable_idx = _build_matrix(bam_path, bam.ref, variants)
    logger.debug(
        f"polyphase_short: ploidy={ploidy}, phaseable_sites={len(phaseable_idx)}, "
        f"reads_in_matrix={matrix.shape[0]}"
    )

    n_phaseable = len(phaseable_idx)
    if matrix.shape[0] < ploidy or n_phaseable < 2:
        _passthrough(in_vcf, out_vcf, bam.sentieon)
        return

    max_allele = int(max(len(variants[i]["alt"]) for i in phaseable_idx))
    labels, _centers = _kmodes(matrix, ploidy, max_allele)
    counts = _per_site_cluster_counts(matrix, labels, ploidy, max_allele)

    # Per-site assignment + confidence
    assigned_per_site: List[List[int]] = []
    confident_per_site = np.zeros(n_phaseable, dtype=bool)
    for c, vi in enumerate(phaseable_idx):
        gt_ms = variants[vi]["gt"]
        if len(gt_ms) != ploidy:
            # Input GT ploidy mismatch — pass through unphased.
            assigned_per_site.append([])
            confident_per_site[c] = False
            continue
        assigned, conf = _assign_clusters_to_gt(counts[:, c, :], list(gt_ms), ploidy)
        assigned_per_site.append(assigned)
        confident_per_site[c] = conf

    connectivity = _connectivity_within_clusters(matrix, labels, ploidy)
    variant_positions = [variants[vi]["pos"] for vi in phaseable_idx]
    ps_per_col = _make_phase_blocks(confident_per_site, connectivity, variant_positions)

    # Emit
    _emit_vcf(in_vcf, out_vcf, variants, phaseable_idx, assigned_per_site, ps_per_col, bam.sentieon)


# -- emit helpers -----------------------------------------------------------
def _emit_vcf(
    in_vcf: str,
    out_vcf: str,
    variants: List[Dict[str, Any]],
    phaseable_idx: List[int],
    assigned_per_site: List[List[int]],
    ps_per_col: List[int],
    sentieon: str,
) -> None:
    # Map (chrom, pos) → (matrix column index, ps, assigned alleles)
    site_info: Dict[Tuple[str, int], Tuple[int, int, List[int]]] = {}
    for c, vi in enumerate(phaseable_idx):
        v = variants[vi]
        site_info[(v["chrom"], v["pos"])] = (c, ps_per_col[c], assigned_per_site[c])

    in_vh = vcflib.VCF(in_vcf)
    out_vh = vcflib.VCF(out_vcf, "w")
    try:
        update_lines = [
            '##FORMAT=<ID=PS,Number=1,Type=Integer,'
            'Description="Phase set identifier (segdup-caller polyphase)">'
        ]
        out_vh.copy_header(in_vh, update=update_lines)
        out_vh.emit_header()
        for v in in_vh:
            info = site_info.get((v.chrom, v.pos))
            if info is not None:
                _, ps, assigned = info
                if ps > 0 and assigned:
                    v.samples[0]["GT"] = "|".join(str(a) for a in assigned)
                    v.samples[0]["PS"] = ps
                    v.line = None
            out_vh.emit(v)
    finally:
        out_vh.close()
        in_vh.close()
    _index_vcf(out_vcf, sentieon)


def _passthrough(in_vcf: str, out_vcf: str, sentieon: str) -> None:
    """Copy input VCF unchanged when nothing can be phased (still need PS header)."""
    in_vh = vcflib.VCF(in_vcf)
    out_vh = vcflib.VCF(out_vcf, "w")
    try:
        out_vh.copy_header(in_vh)
        out_vh.emit_header()
        for v in in_vh:
            out_vh.emit(v)
    finally:
        out_vh.close()
        in_vh.close()
    _index_vcf(out_vcf, sentieon)


def _symlink_through(in_vcf: str, out_vcf: str, sentieon: str) -> None:
    if os.path.exists(out_vcf):
        os.remove(out_vcf)
    os.symlink(os.path.abspath(in_vcf), out_vcf)
    in_idx = in_vcf + ".tbi"
    out_idx = out_vcf + ".tbi"
    if os.path.exists(in_idx):
        if os.path.exists(out_idx):
            os.remove(out_idx)
        os.symlink(os.path.abspath(in_idx), out_idx)
    else:
        _index_vcf(out_vcf, sentieon)


def _index_vcf(vcf_path: str, sentieon: str) -> None:
    cmd = f"{sentieon} util vcfindex {vcf_path}"
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()
