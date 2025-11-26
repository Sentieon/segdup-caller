from genecaller.gene_model import (
    CBSGeneSegmenter,
    HMMSEGGeneSegmenter,
    PELTGeneSegmenter,
    calc_conversion_signals,
)
from genecaller.logging import get_logger
import vcflib
from typing import Optional
from dataclasses import dataclass


class GeneConversionDetector:
    def __init__(self, gene, conversion_config: dict):
        self.gene = gene
        self.config = conversion_config
        self.variants = []
        self.signals = []
        self.breakpoints = []
        self.segmenter = None
        self.logger = get_logger(__name__)
        self._init_from_config()

    def _init_from_config(self):
        self.segment_method = self.config.get("segment_method", "HMMSEG")
        self.segmenter_params = self.config.get("segmenter_params", {})
        self.conversion_region = self.config.get("conversion_region", None)
        self.copy_numbers = self.config.get("copy_numbers", [])
        self.total_copy_number = sum(self.copy_numbers)
        self.min_num_signals = self.config.get("min_num_signals", 2)

        # Fusion detection parameters (default: disabled)
        self.fusion_detection = self.config.get("fusion_detection", False)
        self.gene_strand = self.config.get("gene_strand", None)  # "+" or "-"
        self.three_prime_terminus = self.config.get("three_prime_terminus", None)
        self.five_prime_terminus = self.config.get("five_prime_terminus", None)
        self.recombination_region = self.config.get("recombination_region", None)

        if not self.conversion_region or not self.copy_numbers:
            raise ValueError(
                "Conversion region and copy numbers must be specified in the config."
            )
        self.chrom, self.start, self.end = (
            self.conversion_region.split(":")[0],
            int(self.conversion_region.split(":")[1].split("-")[0]),
            int(self.conversion_region.split(":")[1].split("-")[1]),
        )
        # initial self.segments to be entire region
        self.segments = []

        self.logger.info(
            f"Initialized GeneConversionDetector for region {self.conversion_region}"
        )
        self.logger.debug(f"  Segment method: {self.segment_method}")
        self.logger.debug(
            f"  Copy numbers: {self.copy_numbers} (total: {self.total_copy_number})"
        )
        self.logger.debug(f"  Min num signals: {self.min_num_signals}")
        self.logger.debug(f"  Segmenter params: {self.segmenter_params}")

        # Log fusion detection parameters if enabled
        if self.fusion_detection:
            if self.gene_strand and (
                self.three_prime_terminus or self.five_prime_terminus
            ):
                self.logger.info("  Fusion detection enabled")
                self.logger.debug(f"    Gene strand: {self.gene_strand}")
                if self.three_prime_terminus:
                    self.logger.debug(f"    3' terminus: {self.three_prime_terminus}")
                if self.five_prime_terminus:
                    self.logger.debug(f"    5' terminus: {self.five_prime_terminus}")
                if self.recombination_region:
                    self.logger.debug(
                        f"    Recombination region: {self.recombination_region}"
                    )
            else:
                self.logger.warning(
                    "  Fusion detection enabled but missing required parameters "
                    "(gene_strand and/or terminus coordinates). Fusion detection will be skipped."
                )

    def init_segmenter(self, method, params):
        # Ensure total_cn is in params
        if "total_cn" not in params:
            params["total_cn"] = self.total_copy_number

        # Pass min_num_signals to segmenters
        # For CBS and PELT, this becomes min_width (number of signals)
        # For HMMSEG, used for post-processing merging
        if "min_width" not in params:
            params["min_width"] = self.min_num_signals

        # Add recombination_regions from gene-level config if available
        # Priority: gene.recombination_regions (top-level) > params (from config)
        if self.gene.recombination_regions and "recombination_regions" not in params:
            params["recombination_regions"] = self.gene.recombination_regions

        # Add transition rates from gene-level config if available
        if hasattr(self.gene, "config"):
            if (
                "transition_rate_background" in self.gene.config
                and "transition_rate_background" not in params
            ):
                params["transition_rate_background"] = self.gene.config[
                    "transition_rate_background"
                ]
            if (
                "transition_rate_hotspot" in self.gene.config
                and "transition_rate_hotspot" not in params
            ):
                params["transition_rate_hotspot"] = self.gene.config[
                    "transition_rate_hotspot"
                ]

        if method == "CBS":
            self.segmenter = CBSGeneSegmenter(self.signals, params)
        elif method == "HMMSEG":
            # HMMSEG doesn't use min_width during segmentation, but for post-processing
            params["min_num_signals"] = self.min_num_signals
            self.segmenter = HMMSEGGeneSegmenter(self.signals, params)
        elif method == "PELT":
            self.segmenter = PELTGeneSegmenter(self.signals, params)
        else:
            raise ValueError(f"Unknown segmentation method: {method}")

    def load_variants(self, vcf_files):
        self.logger.debug(f"Loading variants from {len(vcf_files)} VCF files")
        for vcf_file in vcf_files:
            vcf = vcflib.VCF(vcf_file)
            vars_in_region = [
                var for var in vcf.range(self.chrom, self.start, self.end)
            ]
            self.logger.debug(f"  {vcf_file}: {len(vars_in_region)} variants")
            self.variants.extend(vars_in_region)
        self.logger.info(
            f"Loaded {len(self.variants)} total variants in region {self.chrom}:{self.start}-{self.end}"
        )
        liftover_bam = self.gene.read_data["short_read"]["liftover"]
        records = []
        for v in self.variants:
            if not v.id or v.id == ".":
                continue
            if len(v.ref) > 1 or any(len(alt) > 1 for alt in v.alt):
                ads = v.samples[0].get("AD", [])
                if not ads or len(ads) < 2:
                    continue
            else:
                pileuprec = liftover_bam.get_pileup(v.chrom, v.pos)
                if pileuprec is None:
                    continue
                ref_cnt = pileuprec.AD.get(v.ref, 0)
                alt_cnt = pileuprec.AD.get(v.alt[0], 0)
                if sum(pileuprec.AD.values()) - (ref_cnt + alt_cnt) > 2:
                    continue
                ads = [ref_cnt, alt_cnt]
            records.append([v.id, v.pos, ads])
        self.signals = calc_conversion_signals(records)
        self.logger.info(f"Calculated {len(self.signals)} conversion signals")

    def _parse_region(self, region_str: str) -> tuple:
        """
        Parse a genomic region string into (chrom, start, end).

        Args:
            region_str: Region string in format "chr:start-end"

        Returns:
            Tuple of (chrom, start, end)
        """
        chrom = region_str.split(":")[0]
        start = int(region_str.split(":")[1].split("-")[0])
        end = int(region_str.split(":")[1].split("-")[1])
        return chrom, start, end

    def _is_in_recombination_region(self, breakpoint: int) -> bool:
        """Check if breakpoint falls within any recombination region (supports comma-separated regions)."""
        if not self.recombination_region:
            return True  # No constraint, accept all breakpoints

        # Split by comma to handle multiple disjoint regions
        for region_str in self.recombination_region.split(","):
            chrom, coord_str = region_str.split(":")[:2]
            parts = coord_str.strip().split("-")
            rec_start = int(parts[0])
            rec_end = int(parts[1])

            if rec_start <= breakpoint <= rec_end:
                self.logger.debug(
                    f"    Breakpoint {breakpoint} is in recombination region {chrom}:{rec_start}-{rec_end}"
                )
                return True

        self.logger.debug(
            f"    Breakpoint {breakpoint} not in any recombination region: {self.recombination_region}"
        )
        return False

    def _is_fusion_by_position(self, segment, all_conversion_segments):
        """Check if segment is a fusion based on 3'/5' terminus position. Sets segment.fusion_type."""
        if not self.fusion_detection:
            return False

        if not self.gene_strand:
            self.logger.debug("    Fusion check skipped: gene_strand not configured")
            return False

        if not self.three_prime_terminus and not self.five_prime_terminus:
            self.logger.debug(
                "    Fusion check skipped: no terminus coordinates configured"
            )
            return False

        if not all_conversion_segments:
            return False

        # Initialize fusion type to None
        segment.fusion_type = None

        # Check 3' terminus fusion
        if self.three_prime_terminus:
            contains_3prime = segment.start <= self.three_prime_terminus <= segment.end

            if contains_3prime:
                # Determine breakpoint based on strand
                if self.gene_strand == "-":
                    breakpoint = segment.start
                elif self.gene_strand == "+":
                    breakpoint = segment.end
                else:
                    self.logger.warning(
                        f"Unknown gene_strand value: {self.gene_strand}"
                    )
                    return False

                self.logger.debug(
                    f"    3' fusion check ({self.gene_strand} strand): segment=[{segment.start}-{segment.end}], "
                    f"three_prime_terminus={self.three_prime_terminus}, "
                    f"contains_terminus={contains_3prime}, breakpoint={breakpoint}"
                )

                # Check if breakpoint is in any recombination region (supports multiple disjoint regions)
                if self._is_in_recombination_region(breakpoint):
                    segment.fusion_type = "3_PRIME"
                    return True

        # Check 5' terminus fusion
        if self.five_prime_terminus:
            contains_5prime = segment.start <= self.five_prime_terminus <= segment.end

            if contains_5prime:
                # Determine breakpoint based on strand
                if self.gene_strand == "-":
                    breakpoint = segment.end
                elif self.gene_strand == "+":
                    breakpoint = segment.start
                else:
                    self.logger.warning(
                        f"Unknown gene_strand value: {self.gene_strand}"
                    )
                    return False

                self.logger.debug(
                    f"    5' fusion check ({self.gene_strand} strand): segment=[{segment.start}-{segment.end}], "
                    f"five_prime_terminus={self.five_prime_terminus}, "
                    f"contains_terminus={contains_5prime}, breakpoint={breakpoint}"
                )

                # Check if breakpoint is in any recombination region (supports multiple disjoint regions)
                if self._is_in_recombination_region(breakpoint):
                    segment.fusion_type = "5_PRIME"
                    return True

        return False

    def detect_conversions(self):
        if not self.total_copy_number or not self.signals:
            self.logger.warning(
                "No conversion detection: total_copy_number=0 or no signals"
            )
            return []

        self.logger.info(
            f"Starting conversion detection with {self.segment_method} method"
        )

        if self.segmenter is None:
            self.init_segmenter(self.segment_method, self.segmenter_params)
            assert self.segmenter is not None

        self.segmenter.segment()
        self.breakpoints = self.segmenter.breakpoints
        self.logger.info(
            f"Segmentation complete: {len(self.segmenter.segments)} segments, {len(self.breakpoints)} breakpoints"
        )
        self.logger.debug(f"  Breakpoints: {self.breakpoints}")

        # First pass: Collect all conversion segments (with sufficient signals)
        conversion_segments = []
        filtered_count = 0
        outside_gene_count = 0

        for s in self.segmenter.segments:
            num_signals = s.num_sites
            converted_cnt = s.allele_counts[1] - self.copy_numbers[1]

            self.logger.debug(
                f"  Segment {s.start}-{s.end}: allele_counts={s.allele_counts}, "
                f"converted_cnt={converted_cnt}, num_signals={num_signals}"
            )

            # Apply min_num_signals filter
            if num_signals < self.min_num_signals:
                self.logger.debug(
                    f"    Filtered out: {num_signals} signals < min_num_signals={self.min_num_signals}"
                )
                filtered_count += 1
                continue

            # Filter segments completely outside gene boundary (strand-aware)
            if self.three_prime_terminus and self.gene_strand:
                is_outside_gene = False
                if self.gene_strand == "-":
                    # Negative strand: gene extends from three_prime_terminus toward higher coordinates
                    # Filter segments completely below (less than) three_prime_terminus
                    if s.end <= self.three_prime_terminus:
                        is_outside_gene = True
                        self.logger.debug(
                            f"    Filtered out (- strand): segment end {s.end} <= three_prime_terminus {self.three_prime_terminus}"
                        )
                elif self.gene_strand == "+":
                    # Positive strand: gene extends from lower coordinates toward three_prime_terminus
                    # Filter segments completely above (greater than) three_prime_terminus
                    if s.start >= self.three_prime_terminus:
                        is_outside_gene = True
                        self.logger.debug(
                            f"    Filtered out (+ strand): segment start {s.start} >= three_prime_terminus {self.three_prime_terminus}"
                        )

                if is_outside_gene:
                    outside_gene_count += 1
                    continue

            # Keep segments with non-zero conversion count (both positive and negative)
            # Gene-specific interpretation will determine biological meaning
            if converted_cnt != 0:
                s.conversion_allele_count = converted_cnt
                conversion_segments.append(s)

        if filtered_count > 0:
            self.logger.info(
                f"Filtered out {filtered_count} segments with < {self.min_num_signals} signals"
            )
        if outside_gene_count > 0:
            self.logger.info(
                f"Filtered out {outside_gene_count} segments completely outside gene boundary"
            )

        # Second pass: Mark fusion candidates among conversion segments
        for s in conversion_segments:
            is_fusion_candidate = self._is_fusion_by_position(s, conversion_segments)
            s.is_fusion_candidate = is_fusion_candidate

            # Determine event type and include fusion type if applicable
            if is_fusion_candidate:
                fusion_type = getattr(s, "fusion_type", None)
                event_type = f"FUSION ({fusion_type})" if fusion_type else "FUSION"
            else:
                event_type = "CONVERSION"

            self.logger.info(
                f"    {event_type} detected: {s.start}-{s.end}, "
                f"converted {s.conversion_allele_count} alleles, {s.num_sites} signals"
            )

        # Note: Adjacent segment merging is handled in the HMM segmenter
        # (_merge_adjacent_identical_segments in gene_model.py) to prevent
        # the _merge_small_segments post-processing from creating adjacent
        # segments with identical allele counts in the first place.
        self.segments = conversion_segments
        self.logger.info(f"Detected {len(self.segments)} conversion segments")
        return self.segments

    @staticmethod
    def load_known_conversions(db_file: str) -> list["KnownConversion"]:
        """
        Load known gene conversion events from a database file.

        Expected format:
        ##file_format=KnownGeneConversionDB_v1.0
        ##reference_genome=GRCh38
        ##description=A database of known gene conversion events (complex alleles) and their signature variants.
        ##column_definitions=<ID=event_name,Description="Common name or ID for the conversion event">
        ##column_definitions=<ID=recipient_gene,Description="Gene symbol of the recipient locus">
        ##column_definitions=<ID=recipient_locus,Description="Genomic region (chr:start-end) containing all signature variants">
        ##column_definitions=<ID=donor_gene,Description="Gene symbol of the donor locus">
        ##column_definitions=<ID=donor_locus,Description="Genomic region (chr:start-end) of the corresponding donor sequence">
        ##column_definitions=<ID=signature_variants,Description="Pipe-separated list of signature variant IDs">
        ##column_definitions=<ID=event_type,Description="CONVERSION or FUSION (optional, default: CONVERSION)">
        ##column_definitions=<ID=min_match_threshold,Description="Minimum match fraction 0.0-1.0 (optional, default: 0.8)">
        #event_name	recipient_gene	recipient_locus	donor_gene	donor_locus	signature_variants	event_type	min_match_threshold

        Note: event_type and min_match_threshold columns are optional for backward compatibility.

        Args:
            db_file: Path to the known conversion database file

        Returns:
            List of KnownConversion objects
        """
        logger = get_logger(__name__)
        known_conversions = []

        try:
            with open(db_file, "r") as f:
                for line in f:
                    line = line.strip()

                    # Skip empty lines and header lines
                    if (
                        not line
                        or line.startswith("##")
                        or line.startswith("#event_name")
                    ):
                        continue

                    fields = line.split("\t")
                    if len(fields) < 6:
                        logger.warning(
                            f"Skipping malformed line (expected at least 6 fields): {line}"
                        )
                        continue

                    # Required fields
                    event_name = fields[0]
                    recipient_gene = fields[1]
                    recipient_locus = fields[2]
                    donor_gene = fields[3]
                    donor_locus = fields[4]
                    signature_variants = fields[5]

                    # Optional fields (backward compatible)
                    event_type = fields[6] if len(fields) > 6 else "CONVERSION"
                    min_match_threshold = float(fields[7]) if len(fields) > 7 else 0.8

                    # Parse signature variants (pipe-separated list of variant IDs)
                    # Handle "NA" or empty string as no signature variants
                    if signature_variants and signature_variants.upper() != "NA":
                        variant_ids = [
                            v.strip()
                            for v in signature_variants.split("|")
                            if v.strip()
                        ]
                    else:
                        variant_ids = []

                    known_conversion = KnownConversion(
                        event_name=event_name,
                        recipient_gene=recipient_gene,
                        recipient_locus=recipient_locus,
                        donor_gene=donor_gene,
                        donor_locus=donor_locus,
                        signature_variant_ids=variant_ids,
                        event_type=event_type,
                        min_match_threshold=min_match_threshold,
                    )
                    known_conversions.append(known_conversion)

            logger.info(
                f"Loaded {len(known_conversions)} known conversion events from {db_file}"
            )
            return known_conversions

        except FileNotFoundError:
            logger.error(f"Known conversion database file not found: {db_file}")
            return []
        except Exception as e:
            logger.error(f"Error loading known conversion database: {e}")
            return []

    @staticmethod
    def match_known_conversion(
        detected_segment,
        known_conversions: list["KnownConversion"],
        gene_names: list[str],
        match_threshold: float = 0.8,
    ) -> Optional[list[tuple[str, float, set]]]:
        """
        Match detected conversion segment against known conversion events.

        This method handles:
        1. CONVERSION events: Require signature variant matching above threshold
        2. FUSION events: Position-based detection, signature variants optional (for disambiguation)
        3. Nested/overlapping events: Prioritize most comprehensive matches

        Args:
            detected_segment: GeneConversionSegment with signals and is_fusion_candidate flag
            known_conversions: List of KnownConversion objects from the database
            gene_names: List of gene names for current gene (e.g., ["GBA1", "GBAP1"])
            match_threshold: Default minimum match fraction (can be overridden by event-specific threshold)

        Returns:
            List of tuples (event_name, match_score, matched_variant_ids) for all non-overlapping matches
            that meet the threshold, sorted by number of matched variants (descending).
            Returns None if no matches meet the threshold.
        """
        logger = get_logger(__name__)

        if not detected_segment or not known_conversions:
            return None

        # Filter known conversions to match current gene names
        filtered_conversions = [
            kc for kc in known_conversions if kc.recipient_gene in gene_names
        ]

        if not filtered_conversions:
            return None

        # Extract variant IDs from detected segment's signals
        detected_signals = getattr(detected_segment, "signals", [])
        detected_variant_ids = set(
            signal.id for signal in detected_signals if signal.id
        )

        logger.debug(
            f"Matching {len(detected_variant_ids)} detected variants against {len(filtered_conversions)} known events (filtered from {len(known_conversions)} by gene names: {gene_names})"
        )
        logger.debug(
            f"  Segment is fusion candidate: {getattr(detected_segment, 'is_fusion_candidate', False)}"
        )

        # Calculate matches for all events
        all_matches = []
        segment_is_fusion = getattr(detected_segment, "is_fusion_candidate", False)

        for known_conv in filtered_conversions:
            # Use event-specific threshold if available, otherwise use default
            event_threshold = known_conv.min_match_threshold

            # Event type must match segment type
            if known_conv.event_type == "FUSION":
                # FUSION events can only match fusion segments
                if not segment_is_fusion:
                    logger.debug(
                        f"  Skipping {known_conv.event_name} (FUSION): "
                        f"segment not marked as fusion candidate"
                    )
                    continue

                # If no signature variants specified, any fusion matches (for disambiguation only)
                if not known_conv.signature_variant_ids:
                    logger.debug(
                        f"  {known_conv.event_name} (FUSION): matched by position (no signature variants)"
                    )
                    all_matches.append(
                        {
                            "event_name": known_conv.event_name,
                            "event_type": "FUSION",
                            "score": 1.0,
                            "matched_ids": set(),
                            "num_matched": 0,
                            "total_variants": 0,
                        }
                    )
                    continue
            else:
                # CONVERSION events can only match non-fusion segments
                if segment_is_fusion:
                    logger.debug(
                        f"  Skipping {known_conv.event_name} (CONVERSION): "
                        f"segment is marked as fusion candidate"
                    )
                    continue

            # Variant-based matching (for both CONVERSION and FUSION with signature variants)
            signature_ids = set(known_conv.signature_variant_ids)
            matched_ids = signature_ids.intersection(detected_variant_ids)

            if not signature_ids:
                continue

            match_score = len(matched_ids) / len(signature_ids)

            logger.debug(
                f"  {known_conv.event_name} ({known_conv.event_type}): "
                f"{len(matched_ids)}/{len(signature_ids)} variants matched "
                f"(score={match_score:.2f}, threshold={event_threshold:.2f})"
            )

            if match_score >= event_threshold:
                all_matches.append(
                    {
                        "event_name": known_conv.event_name,
                        "event_type": known_conv.event_type,
                        "score": match_score,
                        "matched_ids": matched_ids,
                        "num_matched": len(matched_ids),
                        "total_variants": len(signature_ids),
                    }
                )

        if not all_matches:
            logger.info(f"No matches meet threshold (default={match_threshold})")
            return None

        # Sort by number of matched variants (descending), then by score
        all_matches.sort(key=lambda x: (x["num_matched"], x["score"]), reverse=True)

        # Greedy selection of non-overlapping events
        # Start with the most comprehensive match and exclude any events that share variants
        selected_matches = []
        used_variant_ids = set()

        for match in all_matches:
            # Check if this event's variants overlap with already selected events
            if not match["matched_ids"].intersection(used_variant_ids):
                selected_matches.append(
                    (match["event_name"], match["score"], match["matched_ids"])
                )
                used_variant_ids.update(match["matched_ids"])

                if match["total_variants"] > 0:
                    logger.info(
                        f"Selected match: {match['event_name']} ({match['event_type']}) "
                        f"({match['num_matched']}/{match['total_variants']} variants, score={match['score']:.2f})"
                    )
                else:
                    logger.info(
                        f"Selected match: {match['event_name']} ({match['event_type']}) "
                        f"(position-based, no signature variants)"
                    )
            else:
                overlap_count = len(match["matched_ids"].intersection(used_variant_ids))
                logger.debug(
                    f"Skipping {match['event_name']}: {overlap_count} variants overlap with selected events"
                )

        logger.info(f"Final selection: {len(selected_matches)} non-overlapping events")
        return selected_matches if selected_matches else None


@dataclass
class KnownConversion:
    """Represents a known gene conversion event from the database."""

    event_name: str
    recipient_gene: str
    recipient_locus: str
    donor_gene: str
    donor_locus: str
    signature_variant_ids: list[str]
    event_type: str = "CONVERSION"  # "CONVERSION" or "FUSION"
    min_match_threshold: float = 0.8  # Can be lower for fusions (e.g., 0.6)
