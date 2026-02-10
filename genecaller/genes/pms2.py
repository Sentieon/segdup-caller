"""PMS2 gene with block-based gene conversion interpretation."""

import json
import re
from dataclasses import dataclass
from typing import Dict, Any, List, Optional, Set

from genecaller.gene import Gene
from genecaller.gene_model import GeneConversionSegment
from genecaller.logging import get_logger
from genecaller.util import get_data_file


@dataclass
class ConversionEvent:
    """Interpreted gene conversion event."""

    region: str
    start: int
    end: int
    converted_alleles: int
    direction: str
    event_type: str
    clinical_significance: str
    is_fragmented: bool
    is_recombinant: bool
    block_coverage: Dict[int, float]
    contributing_segments: List[GeneConversionSegment]
    inferred_mechanism: Optional[str] = None
    inferred_swap_region: Optional[str] = None
    is_inferred: bool = False


class PMS2(Gene):
    """PMS2 gene with block-based gene conversion interpretation.

    See docs/PMS2.md for detailed documentation on block structure and classification.
    """

    # Block IDs (B3 is stable/rarely converts; B4 is active/frequently converts)
    HOTSPOT_BLOCK = 2  # Exon 13 (nested in B1)
    LIMIT_BLOCK = 1  # Exon 14/15, 3' terminus
    UPSTREAM_BLOCK = 4  # Active upstream block (B3 stable, B4 active - interlaced)

    # Penetration detection parameters (block coordinates loaded from manifest)
    MIN_PENETRATION_DEPTH = 50
    UPSTREAM_EDGE_ZONE = 2000  # Segment must start within 2kb of upstream block start

    # Coverage thresholds
    HOTSPOT_THRESHOLD = 0.80
    LIMIT_THRESHOLD = 0.80
    DONUT_THRESHOLD = 0.50

    MAX_GAP_TOLERANCE = 500
    TOTAL_CN = 4

    def __init__(self, cfg: dict, ref_file: str) -> None:
        super().__init__(cfg, ref_file)
        self.logger = get_logger(__name__)
        self.manifest_blocks: List[Dict] = []
        self.site_to_block: Dict[int, int] = {}

        # Block coordinates (loaded from manifest)
        self.block2_start: Optional[int] = None
        self.block2_end: Optional[int] = None
        self.upstream_block_start: Optional[int] = None
        self.upstream_block_end: Optional[int] = None

        self.block1_boundary: Optional[int] = None
        if self.conversion_regions:
            self.block1_boundary = self.conversion_regions[0].get(
                "three_prime_terminus"
            )
            if self.block1_boundary:
                self.logger.debug(f"3' terminus boundary: {self.block1_boundary}")

        self._load_block_manifest()

    def _load_block_manifest(self) -> None:
        """Load pms2_psv_manifest.json and build lookup tables."""
        manifest_path = get_data_file("data/pms2_psv_manifest.json")
        if not manifest_path:
            self.logger.warning("PMS2 block manifest not found")
            return

        try:
            with open(manifest_path, "r") as f:
                self.manifest_blocks = json.load(f)

            for block in self.manifest_blocks:
                block_id = block["block_id"]
                for site in block["psv_sites"]:
                    if site not in self.site_to_block:
                        self.site_to_block[site] = block_id
                    elif block_id == self.HOTSPOT_BLOCK:
                        self.site_to_block[site] = block_id

                # Extract block coordinates from manifest
                if block_id == self.HOTSPOT_BLOCK:
                    self.block2_start = block["start_pos"]
                    self.block2_end = block["end_pos"]
                elif block_id == self.UPSTREAM_BLOCK:
                    self.upstream_block_start = block["start_pos"]
                    self.upstream_block_end = block["end_pos"]

            self.logger.info(
                f"Loaded PMS2 block manifest: {len(self.manifest_blocks)} blocks, "
                f"{len(self.site_to_block)} PSV sites"
            )
            if self.block2_start:
                self.logger.debug(
                    f"Hotspot block (B2): {self.block2_start}-{self.block2_end}"
                )
            if self.upstream_block_start:
                self.logger.debug(
                    f"Upstream block (B4): {self.upstream_block_start}-{self.upstream_block_end}"
                )
        except Exception as e:
            self.logger.error(f"Error loading PMS2 block manifest: {e}")
            self.manifest_blocks = []
            self.site_to_block = {}

    @property
    def block2_region(self) -> str:
        """Generate B2 region string from loaded coordinates."""
        if self.block2_start and self.block2_end:
            return f"{self.chr}:{self.block2_start}-{self.block2_end}"
        return ""

    def _get_block_by_id(self, block_id: int) -> Optional[Dict]:
        for block in self.manifest_blocks:
            if block["block_id"] == block_id:
                return block
        return None

    def _convert_to_pseudogene_coords(
        self, start: int, end: int
    ) -> tuple[int, int, str]:
        """Convert PMS2 coordinates to PMS2CL coordinates.

        When conversion direction is PMS2 -> PMS2CL (negative alleles),
        the event is occurring in PMS2CL, so coordinates should be reported
        in PMS2CL space.

        Returns (converted_start, converted_end, region_string).
        """
        try:
            # Use forward mapping (g1tog2): PMS2 -> PMS2CL
            conv_start = self.mapping.get(start, direction="forward")
            conv_end = self.mapping.get(end, direction="forward")
            # PMS2CL coordinates may be in reverse order due to strand
            if conv_start > conv_end:
                conv_start, conv_end = conv_end, conv_start
            region = f"{self.chr}:{conv_start}-{conv_end}"
            return conv_start, conv_end, region
        except (ValueError, KeyError) as e:
            self.logger.warning(f"Could not convert coordinates {start}-{end}: {e}")
            # Fall back to original coordinates
            return start, end, f"{self.chr}:{start}-{end}"

    def _get_b1_subregion_sites(self) -> tuple[List[int], List[int]]:
        """Get B1 sites split by upstream/downstream of B2."""
        b1_block = self._get_block_by_id(self.LIMIT_BLOCK)
        if not b1_block:
            return [], []
        if self.block2_start is None or self.block2_end is None:
            return [], []

        b1_sites = b1_block["psv_sites"]
        upstream = [s for s in b1_sites if s < self.block2_start]
        downstream = [s for s in b1_sites if s > self.block2_end]
        return upstream, downstream

    def _calculate_site_coverage(
        self, sites: List[int], segments: List[GeneConversionSegment]
    ) -> float:
        """Calculate fraction of PSV sites covered by segments."""
        if not sites:
            return 0.0

        covered = 0
        for site in sites:
            for seg in segments:
                if seg.start <= site <= seg.end:
                    covered += 1
                    break
        return covered / len(sites)

    def _detect_donut_hole(
        self,
        segments: List[GeneConversionSegment],
    ) -> Optional[int]:
        """Detect donut-hole pattern using allele mismatch between B2 and B1 borders.

        A donut hole is detected when B1 borders (both upstream and downstream of B2)
        have more converted alleles than B2. This handles cases where there IS coverage
        in B2 but with fewer conversions than expected based on the surrounding borders.

        Prerequisites: Both B1 border regions must have sufficient coverage (DONUT_THRESHOLD)
        to ensure reliable allele count comparison.

        Returns the allele mismatch value (negative) if donut hole detected,
        None otherwise. Only returns when mismatch is negative (fewer conversions
        in B2 than borders). Positive mismatch (more in B2) is expected and ignored.

        The mismatch is computed conservatively as the maximum of (B2 - upstream)
        and (B2 - downstream) to use the least negative value, reducing false positives.
        """
        upstream_sites, downstream_sites = self._get_b1_subregion_sites()
        if not upstream_sites or not downstream_sites:
            return None

        # Require sufficient coverage in both border regions for reliable comparison
        upstream_cov = self._calculate_site_coverage(upstream_sites, segments)
        downstream_cov = self._calculate_site_coverage(downstream_sites, segments)

        if upstream_cov < self.DONUT_THRESHOLD or downstream_cov < self.DONUT_THRESHOLD:
            return None

        if self.block2_start is None or self.block2_end is None:
            return None

        # Get allele counts for B2 and each border region separately
        b2_alleles = []
        upstream_border_alleles = []
        downstream_border_alleles = []

        for seg in segments:
            if seg.conversion_allele_count is None:
                continue

            covers_b2 = seg.start < self.block2_end and seg.end > self.block2_start
            covers_upstream = seg.start <= max(upstream_sites) and seg.end >= min(
                upstream_sites
            )
            covers_downstream = seg.start <= max(downstream_sites) and seg.end >= min(
                downstream_sites
            )

            if covers_b2:
                b2_alleles.append(seg.conversion_allele_count)
            if covers_upstream:
                upstream_border_alleles.append(seg.conversion_allele_count)
            if covers_downstream:
                downstream_border_alleles.append(seg.conversion_allele_count)

        # Need both border regions to have allele data for comparison
        if not upstream_border_alleles or not downstream_border_alleles:
            return None

        b2_count = max(b2_alleles, key=abs) if b2_alleles else 0
        upstream_count = max(upstream_border_alleles, key=abs)
        downstream_count = max(downstream_border_alleles, key=abs)

        # Calculate mismatch for each border, take maximum (least negative = most conservative)
        upstream_mismatch = b2_count - upstream_count
        downstream_mismatch = b2_count - downstream_count
        allele_mismatch = max(upstream_mismatch, downstream_mismatch)

        if allele_mismatch < 0:
            self.logger.info(
                f"Donut-hole pattern (allele mismatch): B2={b2_count}, "
                f"upstream={upstream_count}, downstream={downstream_count} "
                f"-> mismatch={allele_mismatch}"
            )
            return allele_mismatch

        if allele_mismatch > 0:
            self.logger.debug(
                f"No donut-hole (positive mismatch): B2={b2_count}, "
                f"upstream={upstream_count}, downstream={downstream_count} "
                f"-> more conversions in B2 than borders"
            )
        return None

    def _check_upstream_penetration(
        self, segments: List[GeneConversionSegment]
    ) -> bool:
        """Check if any segment penetrates into upstream block (B4) edge zone."""
        if self.upstream_block_start is None:
            return False

        min_penetration_point = self.upstream_block_start + self.MIN_PENETRATION_DEPTH
        max_start_point = self.upstream_block_start + self.UPSTREAM_EDGE_ZONE

        for seg in segments:
            if seg.start <= max_start_point and seg.end >= min_penetration_point:
                self.logger.debug(
                    f"Upstream penetration: segment {seg.start}-{seg.end} "
                    f"(start <= {max_start_point}, end >= {min_penetration_point})"
                )
                return True
        return False

    def _detect_allele_mismatch(
        self, segments: List[GeneConversionSegment]
    ) -> Optional[int]:
        """Detect allele count mismatch between B2-covering and upstream-covering segments.

        When B2 and upstream both have conversion but with different allele counts,
        this may indicate a reciprocal swap in the B2 region.

        Returns the inferred swap allele count if mismatch detected, None otherwise.
        Only returns negative values (fewer conversions in core than fringe),
        as positive values (more conversions in core) are expected and do not
        require inferred reciprocal events.
        """
        if self.upstream_block_start is None:
            return None
        if self.block2_start is None or self.block2_end is None:
            return None

        b2_alleles = [0]
        upstream_alleles = []

        for seg in segments:
            if seg.conversion_allele_count is None:
                continue

            covers_b2 = seg.start < self.block2_end and seg.end > self.block2_start
            in_upstream = (
                seg.start > self.block2_end and seg.end > self.upstream_block_start
            )

            if covers_b2:
                b2_alleles.append(seg.conversion_allele_count)
            elif in_upstream:
                upstream_alleles.append(seg.conversion_allele_count)

        if not upstream_alleles:
            return None

        b2_count = max(b2_alleles, key=abs)
        upstream_count = max(upstream_alleles, key=abs)

        if abs(b2_count) > abs(upstream_count) > 0:
            # More conversions in core (upstream) than fringe (B2) - this is expected
            self.logger.debug(
                f"Allele mismatch (positive): B2={b2_count}, upstream={upstream_count} "
                f"-> more conversions in core than fringe, no inferred event needed"
            )
            return None
        swap_allele = upstream_count - b2_count
        return swap_allele if swap_allele else None

    def _calculate_block_coverage(
        self, segments: List[GeneConversionSegment]
    ) -> Dict[int, float]:
        """Calculate fraction of each block covered by segments."""
        coverage = {}

        for block in self.manifest_blocks:
            block_id = block["block_id"]
            block_sites = set(block["psv_sites"])
            n_sites = block["n_sites"]

            if n_sites == 0:
                coverage[block_id] = 0.0
                continue

            covered_sites: Set[int] = set()
            for seg in segments:
                for site in block_sites:
                    if seg.start <= site <= seg.end:
                        covered_sites.add(site)

            coverage[block_id] = len(covered_sites) / n_sites

        return coverage

    def _classify_event(
        self,
        coverage: Dict[int, float],
        segments: List[GeneConversionSegment],
        direction: str,
    ) -> List[ConversionEvent]:
        """Classify event type from block coverage. See docs/PMS2.md for details."""
        if not segments:
            return []

        # Ensure block2 coordinates are loaded from manifest
        if self.block2_start is None or self.block2_end is None:
            self.logger.warning("Block2 coordinates not loaded - skipping classification")
            return []

        hotspot_cov = coverage.get(self.HOTSPOT_BLOCK, 0.0)
        limit_cov = coverage.get(self.LIMIT_BLOCK, 0.0)
        upstream_cov = coverage.get(self.UPSTREAM_BLOCK, 0.0)

        sorted_segs = sorted(segments, key=lambda s: s.start)
        start = sorted_segs[0].start
        end = max(seg.end for seg in sorted_segs)
        region = f"{self.chr}:{start}-{end}"

        is_contiguous = True
        for i in range(1, len(sorted_segs)):
            gap = sorted_segs[i].start - sorted_segs[i - 1].end
            if gap > self.MAX_GAP_TOLERANCE:
                is_contiguous = False
                self.logger.debug(
                    f"Gap of {gap}bp between {sorted_segs[i - 1].end} and {sorted_segs[i].start}"
                )

        is_recombinant = (
            self.block1_boundary is not None and start <= self.block1_boundary
        )

        allele_counts = [
            seg.conversion_allele_count
            for seg in segments
            if seg.conversion_allele_count is not None
        ]
        if not allele_counts:
            return []
        converted_alleles = max(allele_counts, key=abs)
        is_fragmented = len(segments) > 1

        if direction == "pathogenic":
            direction_label = "PMS2CL->PMS2"
            clinical = "Potential loss of MMR function - relevant for Lynch syndrome"
        else:
            direction_label = "PMS2->PMS2CL"
            clinical = "Functional gain in pseudogene - likely benign"

        events = []

        def add_fragment_annotation(event_type: str) -> str:
            if is_fragmented:
                return event_type + (
                    " [Disconnected Fragments]"
                    if not is_contiguous
                    else " [Fragmented]"
                )
            return event_type

        upstream_penetrated = self._check_upstream_penetration(sorted_segs)

        # A. EXTENDED CONVERSION: B2 high + upstream (B4) penetrated
        if hotspot_cov >= self.HOTSPOT_THRESHOLD and upstream_penetrated:
            # Check for allele mismatch indicating reciprocal swap
            swap_allele = self._detect_allele_mismatch(sorted_segs)

            if swap_allele is not None:
                # Allele mismatch detected - infer reciprocal swap
                # For reciprocal conversion: primary event extends back to B2_START,
                # inferred event spans from B2_START to primary event's original start
                primary_start = self.block2_start
                primary_region = f"{self.chr}:{primary_start}-{end}"
                inferred_end = start  # Original start of the observed segments
                inferred_region = f"{self.chr}:{self.block2_start}-{inferred_end}"

                base_type = "Extended Conversion [Allele Mismatch]"
                if is_recombinant:
                    event_type = f"Recombinant {base_type} ({direction_label})"
                else:
                    event_type = f"{base_type} ({direction_label})"

                primary_event = ConversionEvent(
                    region=primary_region,
                    start=primary_start,
                    end=end,
                    converted_alleles=converted_alleles,
                    direction=direction,
                    event_type=add_fragment_annotation(event_type),
                    clinical_significance=clinical,
                    is_fragmented=is_fragmented,
                    is_recombinant=is_recombinant,
                    block_coverage=coverage,
                    contributing_segments=sorted_segs,
                    inferred_mechanism="Reciprocal Swap (Allele Mismatch)",
                    inferred_swap_region=inferred_region,
                )
                events.append(primary_event)

                # Inferred swap event (opposite sign: if B2 is missing -1, it received +1)
                inferred_allele = -swap_allele
                swap_dir = "pathogenic" if inferred_allele > 0 else "gain"
                swap_dir_label = (
                    "PMS2CL->PMS2" if inferred_allele > 0 else "PMS2->PMS2CL"
                )
                swap_clin = (
                    "Pathogenic (Reciprocal Crossover) - PMS2 promoter drives pseudogene"
                    if inferred_allele > 0
                    else "Functional gain - likely benign"
                )

                inferred_event = ConversionEvent(
                    region=inferred_region,
                    start=self.block2_start,
                    end=inferred_end,
                    converted_alleles=inferred_allele,
                    direction=swap_dir,
                    event_type=f"Inferred Reciprocal Swap (Exon 13) ({swap_dir_label})",
                    clinical_significance=swap_clin,
                    is_fragmented=False,
                    is_recombinant=False,
                    block_coverage={self.HOTSPOT_BLOCK: 1.0},
                    contributing_segments=[],
                    inferred_mechanism="Reciprocal Swap (Allele Mismatch)",
                    inferred_swap_region=inferred_region,
                    is_inferred=True,
                )
                events.append(inferred_event)
                self.logger.info(
                    f"Inferred reciprocal swap from allele mismatch: {inferred_allele:+d}"
                )
                return events

            # No allele mismatch - regular extended conversion
            if is_recombinant:
                event_type = f"Recombinant Extended Conversion (Exon 11-15 + 3' Ext) ({direction_label})"
            else:
                event_type = f"Extended Conversion (Exon 11-15) ({direction_label})"

            events.append(
                ConversionEvent(
                    region=region,
                    start=start,
                    end=end,
                    converted_alleles=converted_alleles,
                    direction=direction,
                    event_type=add_fragment_annotation(event_type),
                    clinical_significance=clinical,
                    is_fragmented=is_fragmented,
                    is_recombinant=is_recombinant,
                    block_coverage=coverage,
                    contributing_segments=sorted_segs,
                )
            )
            return events

        # B. CLASSIC/RECOMBINANT: B2 high only
        if hotspot_cov >= self.HOTSPOT_THRESHOLD:
            # Determine exon range based on actual block coverage
            if limit_cov >= self.LIMIT_THRESHOLD:
                exon_range = "Exon 13-15"  # Both B2 and B1 covered
            else:
                exon_range = "Exon 13"  # Only B2 covered

            if is_recombinant:
                event_type = (
                    f"Recombinant Conversion (3' Extension) ({direction_label})"
                )
            else:
                event_type = f"Classic Conversion ({exon_range}) ({direction_label})"

            events.append(
                ConversionEvent(
                    region=region,
                    start=start,
                    end=end,
                    converted_alleles=converted_alleles,
                    direction=direction,
                    event_type=add_fragment_annotation(event_type),
                    clinical_significance=clinical,
                    is_fragmented=is_fragmented,
                    is_recombinant=is_recombinant,
                    block_coverage=coverage,
                    contributing_segments=sorted_segs,
                )
            )
            return events

        # C. INFERRED: upstream penetrated without B2 -> infer reciprocal swap
        if upstream_penetrated:
            donut_mismatch = self._detect_donut_hole(sorted_segs)
            mismatch_allele = self._detect_allele_mismatch(sorted_segs)

            # Determine mismatch source: donut hole takes precedence, then allele mismatch
            if donut_mismatch is not None:
                subtype = "[Donut-Hole Pattern]"
                inferred_note = "Definitive reciprocal swap"
                effective_mismatch = donut_mismatch
            elif mismatch_allele is not None:
                subtype = "[Allele Mismatch]"
                inferred_note = "Reciprocal swap from allele mismatch"
                effective_mismatch = mismatch_allele
            else:
                # No mismatch detected - upstream penetrated but no evidence of swap
                # This is expected pattern, no inferred event needed
                base_type = "Extended Conversion"
                if limit_cov >= self.LIMIT_THRESHOLD and is_recombinant:
                    event_type = f"Recombinant {base_type} ({direction_label})"
                else:
                    event_type = f"{base_type} ({direction_label})"

                events.append(
                    ConversionEvent(
                        region=region,
                        start=start,
                        end=end,
                        converted_alleles=converted_alleles,
                        direction=direction,
                        event_type=add_fragment_annotation(event_type),
                        clinical_significance=clinical,
                        is_fragmented=is_fragmented,
                        is_recombinant=is_recombinant,
                        block_coverage=coverage,
                        contributing_segments=sorted_segs,
                    )
                )
                return events

            # For reciprocal conversion: primary event extends back to B2_START,
            # inferred event spans from B2_START to primary event's original start
            primary_start = self.block2_start
            primary_region = f"{self.chr}:{primary_start}-{end}"
            inferred_end = start  # Original start of the observed segments
            inferred_region = f"{self.chr}:{self.block2_start}-{inferred_end}"

            base_type = "Extended Conversion"
            if limit_cov >= self.LIMIT_THRESHOLD and is_recombinant:
                event_type = f"Recombinant {base_type} {subtype} ({direction_label})"
            else:
                event_type = f"{base_type} {subtype} ({direction_label})"

            primary_event = ConversionEvent(
                region=primary_region,
                start=primary_start,
                end=end,
                converted_alleles=converted_alleles,
                direction=direction,
                event_type=add_fragment_annotation(event_type),
                clinical_significance=clinical,
                is_fragmented=is_fragmented,
                is_recombinant=is_recombinant,
                block_coverage=coverage,
                contributing_segments=sorted_segs,
                inferred_mechanism="Reciprocal Swap",
                inferred_swap_region=inferred_region,
            )
            events.append(primary_event)

            # Inferred swap - use negated mismatch allele
            # (opposite sign: if B2 is missing -1, it received +1)
            inferred_allele = -effective_mismatch
            swap_dir = "pathogenic" if inferred_allele > 0 else "gain"
            swap_dir_label = "PMS2CL->PMS2" if inferred_allele > 0 else "PMS2->PMS2CL"
            swap_clin = (
                "Pathogenic (Reciprocal Crossover) - PMS2 promoter drives pseudogene"
                if inferred_allele > 0
                else "Functional gain - likely benign"
            )

            inferred_event = ConversionEvent(
                region=inferred_region,
                start=self.block2_start,
                end=inferred_end,
                converted_alleles=inferred_allele,
                direction=swap_dir,
                event_type=f"Inferred Reciprocal Swap (Exon 13) ({swap_dir_label})",
                clinical_significance=swap_clin,
                is_fragmented=False,
                is_recombinant=False,
                block_coverage={self.HOTSPOT_BLOCK: 1.0},
                contributing_segments=[],
                inferred_mechanism="Reciprocal Swap",
                inferred_swap_region=inferred_region,
                is_inferred=True,
            )
            events.append(inferred_event)
            self.logger.info(f"Inferred reciprocal swap in B2: {inferred_note}")
            return events

        # D. B1 high with donut-hole pattern
        if limit_cov >= self.LIMIT_THRESHOLD:
            donut_mismatch = self._detect_donut_hole(sorted_segs)

            if donut_mismatch is not None:
                # Donut-hole pattern: B1 borders converted, B2 is the hole
                # Primary event: observed B1 conversion (original region)
                # Inferred event: reciprocal swap in B2 (the hole)
                event_type = (
                    f"Partial Conversion [Donut-Hole Pattern] ({direction_label})"
                )
                primary_event = ConversionEvent(
                    region=region,
                    start=start,
                    end=end,
                    converted_alleles=converted_alleles,
                    direction=direction,
                    event_type=add_fragment_annotation(event_type),
                    clinical_significance="Partial event with definitive reciprocal swap signature",
                    is_fragmented=is_fragmented,
                    is_recombinant=is_recombinant,
                    block_coverage=coverage,
                    contributing_segments=sorted_segs,
                    inferred_mechanism="Reciprocal Swap",
                    inferred_swap_region=self.block2_region,
                )
                events.append(primary_event)

                # Inferred swap in B2 (the hole) - use negated donut mismatch
                # (opposite sign: if B2 is missing -1, it received +1)
                inferred_allele = -donut_mismatch
                swap_dir = "pathogenic" if inferred_allele > 0 else "gain"
                swap_dir_label = (
                    "PMS2CL->PMS2" if inferred_allele > 0 else "PMS2->PMS2CL"
                )
                swap_clin = (
                    "Pathogenic (Reciprocal Crossover) - PMS2 promoter drives pseudogene"
                    if inferred_allele > 0
                    else "Functional gain - likely benign"
                )

                inferred_event = ConversionEvent(
                    region=self.block2_region,
                    start=self.block2_start,
                    end=self.block2_end,
                    converted_alleles=inferred_allele,
                    direction=swap_dir,
                    event_type=f"Inferred Reciprocal Swap (Exon 13) ({swap_dir_label})",
                    clinical_significance=swap_clin,
                    is_fragmented=False,
                    is_recombinant=False,
                    block_coverage={self.HOTSPOT_BLOCK: 1.0},
                    contributing_segments=[],
                    inferred_mechanism="Reciprocal Swap",
                    inferred_swap_region=self.block2_region,
                    is_inferred=True,
                )
                events.append(inferred_event)
                return events

            # E. Partial tail without donut-hole
            if is_recombinant:
                event_type = (
                    f"Recombinant Conversion (3' Extension) ({direction_label})"
                )
                clinical_sig = clinical
            else:
                event_type = f"Partial Conversion (3' Tail) ({direction_label})"
                clinical_sig = "Partial event in 3' region - may need confirmation"

            events.append(
                ConversionEvent(
                    region=region,
                    start=start,
                    end=end,
                    converted_alleles=converted_alleles,
                    direction=direction,
                    event_type=add_fragment_annotation(event_type),
                    clinical_significance=clinical_sig,
                    is_fragmented=is_fragmented,
                    is_recombinant=is_recombinant,
                    block_coverage=coverage,
                    contributing_segments=sorted_segs,
                )
            )
            return events

        # No significant event
        self.logger.debug(
            f"No significant event: Block coverage B1={limit_cov:.2f}, "
            f"B2={hotspot_cov:.2f}, B4={upstream_cov:.2f}"
        )
        return []

    def _interpret_conversion_segments(
        self, raw_segments: List[GeneConversionSegment]
    ) -> List[ConversionEvent]:
        """Group segments by direction and apply block-based classification."""
        if not self.manifest_blocks:
            self.logger.warning("No block manifest - skipping interpretation")
            return []

        if not raw_segments:
            return []

        pathogenic_segs = [
            s
            for s in raw_segments
            if s.conversion_allele_count is not None and s.conversion_allele_count > 0
        ]
        gain_segs = [
            s
            for s in raw_segments
            if s.conversion_allele_count is not None and s.conversion_allele_count < 0
        ]

        self.logger.info(
            f"PMS2: {len(pathogenic_segs)} pathogenic, {len(gain_segs)} gain segments"
        )

        events = []

        if pathogenic_segs:
            coverage = self._calculate_block_coverage(pathogenic_segs)
            self.logger.debug(f"Pathogenic coverage: {coverage}")
            for event in self._classify_event(coverage, pathogenic_segs, "pathogenic"):
                events.append(event)
                self.logger.info(
                    f"  {'Inferred' if event.is_inferred else 'Pathogenic'}: {event.event_type}"
                )

        if gain_segs:
            coverage = self._calculate_block_coverage(gain_segs)
            self.logger.debug(f"Gain coverage: {coverage}")
            for event in self._classify_event(coverage, gain_segs, "gain"):
                events.append(event)
                self.logger.info(
                    f"  {'Inferred' if event.is_inferred else 'Gain'}: {event.event_type}"
                )

        return events

    def prepare_output(self) -> Dict[str, Any]:
        """Prepare output with PMS2-specific block-based interpretation."""
        result_data = super().prepare_output()
        interpretation_lines = []

        pms2_cn = result_data["Copy numbers"].get("PMS2", 2)
        pms2cl_cn = result_data["Copy numbers"].get("PMS2CL", 2)
        interpretation_lines.append(f"PMS2: CN={pms2_cn}")
        interpretation_lines.append(f"PMS2CL: CN={pms2cl_cn}")

        interpreted_events = self._interpret_conversion_segments(
            self.converted_segments
        )

        if interpreted_events:
            result_data["gene_conversions"] = []
            interpretation_lines.append(
                f"Gene Conversion Events: {len(interpreted_events)}"
            )

            for i, event in enumerate(interpreted_events, 1):
                # Convert coordinates to PMS2CL for gain events (negative alleles)
                # The event is occurring in PMS2CL when PMS2 sequence replaces PMS2CL
                if event.converted_alleles < 0:
                    _, _, display_region = self._convert_to_pseudogene_coords(
                        event.start, event.end
                    )
                    pms2_ref_region = event.region  # Keep original PMS2 coordinates
                    if event.inferred_swap_region:
                        # Parse and convert inferred_swap_region too
                        match = re.match(r"chr\d+:(\d+)-(\d+)", event.inferred_swap_region)
                        if match:
                            swap_start, swap_end = int(match.group(1)), int(match.group(2))
                            _, _, display_swap_region = self._convert_to_pseudogene_coords(
                                swap_start, swap_end
                            )
                            pms2_ref_swap_region = event.inferred_swap_region
                        else:
                            display_swap_region = event.inferred_swap_region
                            pms2_ref_swap_region = None
                    else:
                        display_swap_region = None
                        pms2_ref_swap_region = None
                else:
                    display_region = event.region
                    pms2_ref_region = None
                    display_swap_region = event.inferred_swap_region
                    pms2_ref_swap_region = None

                conversion_data = {
                    "region": display_region,
                    "converted_alleles": int(event.converted_alleles),
                    "is_fusion_candidate": event.is_recombinant,
                    "interpretation": event.event_type,
                    "clinical_significance": event.clinical_significance,
                    "is_inferred": event.is_inferred,
                }

                # Include original PMS2 coordinates as reference for gain events
                if pms2_ref_region:
                    conversion_data["pms2_reference_region"] = pms2_ref_region

                if event.inferred_mechanism:
                    conversion_data["inferred_mechanism"] = event.inferred_mechanism
                if display_swap_region:
                    conversion_data["inferred_swap_region"] = display_swap_region
                if pms2_ref_swap_region:
                    conversion_data["pms2_reference_swap_region"] = pms2_ref_swap_region

                conversion_sites = []
                for seg in event.contributing_segments:
                    conversion_sites.extend(
                        [s.id for s in seg.signals if s.is_conversion_site()]
                    )
                conversion_data["conversion_sites"] = list(set(conversion_sites))

                result_data["gene_conversions"].append(conversion_data)

                if event.is_inferred:
                    interpretation_lines.append(
                        f"  Event {i} (Inferred): {event.event_type}"
                    )
                else:
                    interpretation_lines.append(f"  Event {i}: {event.event_type}")

                interpretation_lines.append(f"    Region: {display_region}")
                if pms2_ref_region:
                    interpretation_lines.append(f"    PMS2 Reference: {pms2_ref_region}")
                interpretation_lines.append(
                    f"    Allele change: {event.converted_alleles:+d}"
                )
                interpretation_lines.append(
                    f"    Clinical: {event.clinical_significance}"
                )

                if not event.is_inferred:
                    block_names = {1: "B1", 2: "B2", 4: "B4"}
                    blocks_str = ", ".join(
                        f"{block_names.get(b, f'B{b}')}={c:.2f}"
                        for b, c in sorted(event.block_coverage.items())
                        if c > 0.01 and b <= 4
                    )
                    if blocks_str:
                        interpretation_lines.append(f"    Evidence: {blocks_str}")
                    if event.inferred_mechanism:
                        inferred_region_display = display_swap_region or event.inferred_swap_region
                        interpretation_lines.append(
                            f"    Inferred: {event.inferred_mechanism} in {inferred_region_display}"
                        )
                        if pms2_ref_swap_region:
                            interpretation_lines.append(
                                f"    PMS2 Reference (swap): {pms2_ref_swap_region}"
                            )

                if event.is_inferred:
                    interpretation_lines.append(
                        "    Note: Derived from artifactual gap pattern"
                    )

                if event.is_recombinant and not event.is_inferred:
                    interpretation_lines.append(
                        "    Recombinant: Extends to 3' terminus (Exon 14/15)"
                    )

                if event.is_fragmented:
                    interpretation_lines.append(
                        f"    Note: Event reconstructed from {len(event.contributing_segments)} fragments"
                    )

        elif self.converted_segments:
            interpretation_lines.append(
                "Gene Conversion: Raw signals below block thresholds"
            )
            interpretation_lines.append(
                f"  ({len(self.converted_segments)} raw segment(s) - may be noise)"
            )
        else:
            interpretation_lines.append("No gene conversion events detected")

        result_data["cn_interpretation"] = "\n".join(interpretation_lines)

        return result_data
