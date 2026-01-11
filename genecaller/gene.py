from .util import IntervalList, GeneMapping, get_data_file, Reference
from .logging import get_logger
from .gene_model import create_gene_priors, GenePriors, GeneConversionSegment
from .conversion_detector import GeneConversionDetector
import copy
from genecaller.bam_process import CopyNumberModel, Phased_vcf
from concurrent.futures import ThreadPoolExecutor
import heapq
import vcflib
from typing import List, Dict, Tuple, Any, Optional
import numpy as np


class GeneRegion:
    cn = 2
    cn_diff = 0
    total_cn = 2
    total_cn_diff = 0
    var_factor = 1.0

    def __init__(
        self,
        region: str,
        region_id: int,
        name: str = "",
        dup: int = 1,
        longdel_id: int = -1,
    ) -> None:
        self.region = region
        self.region_id = 1 << region_id if region_id >= 0 else 0
        self.name = name
        self.dup = dup
        self.longdel_id = 1 << longdel_id if longdel_id >= 0 else 0  # 0: no deletion
        self.called_vcf = {}
        self.phased_vcf = {}
        self.matched_region = []
        self.orig_region_id = region_id
        self.orig_longdel_id = longdel_id

    # merge with another region if they are the same region_id, dup, longdel_id and contiguous
    # return True if merged, False otherwise
    def merge(self, region: "GeneRegion") -> bool:
        if (
            self.region_id != region.region_id
            or self.dup != region.dup
            or self.longdel_id != region.longdel_id
        ):
            return False
        merged = (
            IntervalList(region=self.region)
            .union(region.region)
            .to_region_str()
            .split(",")
        )
        if len(merged) != 1:
            return False
        self.region = ",".join(merged)
        return True

    # test if two regions overlap
    # return True if overlap, False otherwise
    def is_overlap(self, region: str) -> bool:
        overlap = (
            IntervalList(region=self.region)
            .intersect(IntervalList(region=region))
            .to_region_str()
            .split(",")
        )
        return len(overlap) > 0 and overlap[0] != ""

    # Perform liftover for the gene region
    # Modify the region in place
    def liftover(self, mapping: "GeneMapping") -> None:
        try:
            r = mapping.interval(self.region, direction="auto")
        except Exception:
            raise RuntimeError("Failed to liftover")
        self.region = r

    # return size of the region
    def size(self) -> int:
        _, s, e = self.get_region()
        return e - s

    # return chrom, start, end of the region
    def get_region(self) -> Tuple[str, int, int]:
        flds = self.region.split(":")
        c = flds[0]
        locs = flds[1].split("-")
        return c, int(locs[0]), int(locs[1])

    # set region to c:s-e
    def set_region(self, c: str, s: int, e: int) -> None:
        self.region = f"{c}:{s}-{e}"

    # Insert region2 into region1 (region2 should be a subset of region1)
    # Return the resulting up to three sub-regions of region1
    @staticmethod
    def split_by_region(
        region1: "GeneRegion", region2: "GeneRegion", min_size: int = 50
    ) -> List["GeneRegion"]:
        if region2.size() < min_size:
            return [region1]
        diff_region = (
            IntervalList(region=region1.region)
            .subtract(region2.region)
            .to_region_str()
            .split(",")
        )
        if region2.dup < 0:
            region2.dup = region1.dup
        result_regions = [region2]
        for dr in diff_region:
            if not dr:
                continue
            r = copy.deepcopy(region1)
            r.region = dr
            result_regions.append(r)
        regions = [
            r
            for r in sorted(result_regions, key=lambda r: r.region)
            if r.size() > min_size
        ]
        c, s, e = region1.get_region()
        _, s1, e1 = regions[0].get_region()
        _, s2, e2 = regions[-1].get_region()
        if s1 != s:
            regions[0].set_region(c, s, e1)
        if e2 != e:
            regions[-1].set_region(c, s2, e)
        return regions

    @staticmethod
    def remove_region_by_longdel(
        regions: List["GeneRegion"], longdel_id: int
    ) -> List["GeneRegion"]:
        del_mask = ~(1 << longdel_id)
        i = 0
        while i < len(regions):
            r = regions[i]
            if r.longdel_id >> longdel_id & 1:
                # Clear the longdel_id bit
                r.longdel_id &= del_mask
                prev_merged = next_merged = False
                if i > 0 and regions[i - 1].merge(r):
                    r = regions[i - 1]
                    prev_merged = True
                if i < len(regions) - 1 and r.merge(regions[i + 1]):
                    next_merged = True

                if prev_merged and next_merged:
                    # Merge all three regions: prev + current + next.
                    regions.pop(i + 1)  # Remove next
                    regions.pop(i)  # remove current
                    continue
                if prev_merged:
                    # Merged with previous, remove current.
                    regions.pop(i)
                    continue
                if next_merged:
                    # Merged with next, remove next
                    regions.pop(i + 1)
                    # Don't increment i, check the same position again
                    continue
                # If no merge possible, leave the region and move to next
            i += 1
        return regions

    @staticmethod
    def insert(
        regions: List["GeneRegion"],
        region_to_insert: "GeneRegion",
        ignore_region_id: bool = False,
        min_overlap: int = 100,
        keep_orig_id: bool = False,
        use_orig: bool = False,
    ) -> List["GeneRegion"]:
        result_regions = []
        overlaps = []
        for i, r in enumerate(regions):
            ovl = IntervalList(region=r.region).intersect(region_to_insert.region)
            if ovl.size() > min_overlap:
                overlaps.append((i, ovl.to_region_str()))
        if not overlaps:
            return regions
        last_i = -1
        for ovl in overlaps:
            i = ovl[0]
            if use_orig:
                other_region = copy.deepcopy(region_to_insert)
                region = copy.deepcopy(regions[i])
            else:
                other_region = copy.deepcopy(regions[i])
                region = copy.deepcopy(region_to_insert)
            region.region = ovl[1]
            if keep_orig_id:
                region.orig_region_id = (
                    [other_region.region_id, region.region_id]
                    if not use_orig
                    else [region.region_id, other_region.region_id]
                )
                region.orig_longdel_id = (
                    [other_region.longdel_id, region.longdel_id]
                    if not use_orig
                    else [region.longdel_id, other_region.longdel_id]
                )
            if not use_orig:
                if ignore_region_id:
                    region.region_id = other_region.region_id
                else:
                    region.region_id |= other_region.region_id
                region.longdel_id |= other_region.longdel_id
            region.dup = other_region.dup
            if region.called_vcf or other_region.called_vcf:
                for seq_key, vcf_dict in other_region.called_vcf.items():
                    if seq_key in region.called_vcf:
                        region.called_vcf[seq_key].update(vcf_dict)
                    else:
                        region.called_vcf[seq_key] = copy.deepcopy(vcf_dict)
            if region.phased_vcf or other_region.phased_vcf:
                for seq_key, vcf_dict in other_region.phased_vcf.items():
                    if seq_key in region.phased_vcf:
                        region.phased_vcf[seq_key].update(vcf_dict)
                    else:
                        region.phased_vcf[seq_key] = copy.deepcopy(vcf_dict)
            rs = GeneRegion.split_by_region(
                regions[i], region, min_size=int(min_overlap / 2)
            )
            result_regions += regions[last_i + 1 : i] + rs
            last_i = i
        result_regions += regions[last_i + 1 :]
        return result_regions

    @staticmethod
    def merge_by_cn(regions: list["GeneRegion"]) -> list["GeneRegion"]:
        if not regions:
            return []

        final_regions = []
        i = 0

        while i < len(regions):
            current_region = copy.deepcopy(regions[i])
            c, s, e = current_region.get_region()

            # Find all consecutive regions with the same CN and merge them
            j = i + 1
            while j < len(regions):
                next_region = regions[j]
                _, s1, e1 = next_region.get_region()

                # Check if regions are adjacent and have same CN
                if next_region.cn == current_region.cn and abs(s1 - e) < 10:
                    # Merge the regions
                    current_region.region_id |= next_region.region_id
                    current_region.longdel_id |= next_region.longdel_id
                    e = e1  # Extend the end coordinate
                    j += 1
                else:
                    break

            # Set the final merged region coordinates
            current_region.region = f"{c}:{s}-{e}"
            final_regions.append(current_region)

            # Move to the next unprocessed region
            i = j

        return final_regions

    @staticmethod
    def get_regions_by_coord(
        regions: List["GeneRegion"], reg: str
    ) -> List["GeneRegion"]:
        region = IntervalList(reg)
        return [r for r in regions if region.intersect(r.region).size()]


class Gene:
    high_mq_regions = {}
    cn_prior_std = 2
    conversion_db = None  # Path to known gene conversion database
    _log_priors = {}  # DEPRECATED: kept for backward compatibility
    _longdel_log_priors = {}  # Pre-calculated log priors for long deletions
    _priors: Optional[GenePriors] = None  # GenePriors instance for CN prior calculation
    converted_segments: List[GeneConversionSegment] = []

    def __init__(self, cfg: dict, ref_file: str) -> None:
        self.config = cfg.get("config", {})
        self.ref = Reference(ref_file)
        self.gene_strand = cfg["gene_strand"]
        self.gene_names = cfg["gene_names"]
        self.gene_regions = cfg["gene_regions"]
        self.conversion_regions = cfg.get("conversion_regions", [])
        self.recombination_regions = cfg.get("recombination_regions", None)
        self.result_vcfs = []
        self.diff_vcf = [get_data_file(p) for p in cfg["diff_vcf"]]
        assert all(self.diff_vcf)
        mapping_file = get_data_file(cfg["map"])
        assert mapping_file is not None, f"Mapping file not found: {cfg['map']}"
        self.mapping = GeneMapping(mapping_file, ref_genome=self.ref)
        if "liftover_regions" in cfg:
            self.liftover_region_names = []
            self.liftover_target_regions = []
            for n, r in cfg["liftover_regions"]:
                self.liftover_region_names.append(n)
                self.liftover_target_regions.append(r)
        else:
            self.liftover_region_names = self.gene_names
            self.liftover_target_regions = self.gene_regions

        # Parse and validate exclusion regions (optional)
        self.exclusion_regions = cfg.get("exclusion_regions", [])
        if self.exclusion_regions:
            # Validate: exclusions must be outside liftover regions
            liftover_itv = IntervalList(",".join(self.liftover_target_regions))
            exclusion_itv = IntervalList(",".join(self.exclusion_regions))
            overlap = liftover_itv.intersect(exclusion_itv)
            if overlap.size() > 0:
                raise ValueError(
                    f"Exclusion regions overlap with liftover regions at {overlap.to_region_str()}. "
                    f"Exclusions are only allowed outside liftover_target_regions."
                )

        total_dup_region = ",".join(self.liftover_target_regions)
        self.nonliftover_region = [
            IntervalList(region=g).subtract(total_dup_region).to_region_str()
            for _, g in enumerate(self.gene_regions)
        ]

        # Subtract exclusion regions from nonliftover regions
        if self.exclusion_regions:
            exclusion_str = ",".join(self.exclusion_regions)
            self.nonliftover_region = [
                IntervalList(region=r).subtract(exclusion_str).to_region_str()
                if r
                else ""
                for r in self.nonliftover_region
            ]
        self.liftover_regions = [
            ",".join([g for j, g in enumerate(self.liftover_target_regions) if j != i])
            for i in range(len(self.liftover_target_regions))
        ]
        self.realign_regions = [
            IntervalList(region=g).pad(cfg.get("padding", 200)).to_region_str()
            for g in self.liftover_target_regions
        ]
        self.setup_regions(cfg)
        self.chr = cfg["gene_regions"][0].split(",")[0].split(":")[0]
        self.variant_call_threads = 4  # Default to 4 threads for vara
        self.logger = get_logger(self.__class__.__name__)
        self._initialize_priors()
        self._initialize_longdel_priors(cfg)

    def _initialize_priors(self) -> None:
        """Initialize CN prior model using factory pattern."""
        # Create GenePriors instance using factory
        self._priors = create_gene_priors(self.config, self.all_vars["cns"])

    def _initialize_longdel_priors(self, cfg: dict) -> None:
        """Pre-calculate log priors for long deletions"""
        import numpy as np

        # Check if gene-specific longdel priors are defined
        if "longdel_priors" in self.config:
            longdel_priors = self.config["longdel_priors"]

            # Validate that we have longdel definitions
            if "longdel" not in cfg:
                self.logger.warning(
                    "longdel_priors specified but no longdel definitions found"
                )
                return

            # Get longdel names from config
            longdel_names = [del_info[0] for del_info in cfg["longdel"]]

            # Check for missing or extra priors
            missing_priors = set(longdel_names) - set(longdel_priors.keys())
            extra_priors = set(longdel_priors.keys()) - set(longdel_names)

            if missing_priors:
                self.logger.warning(f"Missing longdel priors for: {missing_priors}")
            if extra_priors:
                self.logger.warning(
                    f"Extra longdel priors (not in longdel list): {extra_priors}"
                )

            # Calculate log probabilities for existing deletions
            for longdel_name in longdel_names:
                if longdel_name in longdel_priors:
                    prob = longdel_priors[longdel_name]
                    if prob <= 0:
                        self.logger.warning(
                            f"Longdel prior for {longdel_name} is {prob} <= 0. Setting to very small value."
                        )
                        self._longdel_log_priors[longdel_name] = -1e9
                    else:
                        self._longdel_log_priors[longdel_name] = np.log(prob)
                else:
                    # Use default small probability for missing priors
                    self._longdel_log_priors[longdel_name] = np.log(0.001)

            self.logger.debug(f"Using longdel priors: {longdel_priors}")
        else:
            # No longdel priors specified - use uniform or default
            if "longdel" in cfg:
                longdel_names = [del_info[0] for del_info in cfg["longdel"]]
                default_prob = 0.01  # Default 1% probability for each deletion
                for longdel_name in longdel_names:
                    self._longdel_log_priors[longdel_name] = np.log(default_prob)
                self.logger.debug(
                    f"Using default longdel priors: {default_prob} for each deletion"
                )

    def setup_regions(self, config: dict) -> None:
        self.all_regions = []  # full list of GeneRegion objects
        self.segdup_regions = []  # list of GeneRegion objects for all segdup regions
        self.liftover_segdup_regions = []  # list of GeneRegion objects for all segdup regions after liftover
        ignore_nondup = config.get("ignore_nondup", True)
        no_custom = "cn_regions" not in config
        for i, r in enumerate(self.gene_regions):
            region_id = i + 1 if no_custom else 0
            # if "liftover_regions" in config and ignore_nondup:
            #     region_id = -1
            self.all_regions.append(GeneRegion(r, region_id, self.gene_names[i], i + 1))
        self.all_vars = {"cns": {}, "dels": {}}
        if not ignore_nondup:
            self.all_vars["cns"][0] = {"name": "non-duplication", "cn_diff": 0}
        if no_custom:
            for i in range(len(self.liftover_target_regions)):
                self.all_vars["cns"][i + 1] = {
                    "name": self.liftover_region_names[i],
                    "cn_diff": 0,
                }
        num_nondup = 0 if ignore_nondup else len(self.gene_regions)
        if "liftover_regions" in config:
            # change all gene regions to non-dup by default
            for r in self.all_regions:
                r.dup = 0
            liftover_regions = [r.split(",") for r in self.liftover_target_regions]
            for j in range(len(liftover_regions[0])):
                regions = []
                sd_regions = []
                for i, name in enumerate(self.liftover_region_names):
                    n = name if not j else f"{name}{j}"
                    regions.append(
                        GeneRegion(
                            liftover_regions[i][j],
                            (num_nondup + i + 1) if no_custom else 0,
                            n,
                            i + 1,
                        )
                    )
                    sd_regions.append(
                        GeneRegion(
                            liftover_regions[i][j],
                            (num_nondup + i + 1) if no_custom else 0,
                            n,
                            i + 1,
                        )
                    )
                for r in regions:
                    self.all_regions = GeneRegion.insert(self.all_regions, r)
                if not self.segdup_regions:
                    self.segdup_regions = [[r] for r in sd_regions]
                else:
                    for k, r in enumerate(sd_regions):
                        self.segdup_regions[k].append(r)
        else:
            for r in self.all_regions:
                self.segdup_regions.append([copy.deepcopy(r)])
        self.orig_segdup_regions = copy.deepcopy(self.segdup_regions)
        # long deletions
        self.longdel_names = []
        if "longdel" in config:
            longdel = config["longdel"]
            del_id = 0
            for del_name, longdel in longdel:
                del_region = GeneRegion(longdel, 0, del_name, longdel_id=del_id)
                self.all_regions = GeneRegion.insert(
                    self.all_regions, del_region, ignore_region_id=True
                )
                self.segdup_regions = [
                    GeneRegion.insert(r, del_region, ignore_region_id=True)
                    for r in self.segdup_regions
                ]
                self.all_vars["dels"][del_id] = {"name": del_name, "cn_diff": 0}
                self.longdel_names.append(del_name)
                del_id += 1
        # CN-specific regions
        if "cn_regions" in config:
            self.cn_regions = config["cn_regions"]
            region_id = 1
            for name, r in config["cn_regions"]:
                region = GeneRegion(r, region_id, name, -1)
                self.all_regions = GeneRegion.insert(self.all_regions, region)
                for i, r in enumerate(self.segdup_regions):
                    self.segdup_regions[i] = GeneRegion.insert(r, region)
                self.all_vars["cns"][region_id] = {"name": name, "cn_diff": 0}
                region_id += 1
        # merge liftover segdup regions
        for i, regions in enumerate(self.segdup_regions):
            lo_regs = [copy.deepcopy(r) for r in regions]
            for j, regs in enumerate(self.segdup_regions):
                if i == j:
                    continue
                for r in regs:
                    r = copy.deepcopy(r)
                    r.liftover(self.mapping)
                    lo_regs = GeneRegion.insert(lo_regs, r, keep_orig_id=True)
            self.liftover_segdup_regions.append(lo_regs)

    def prepare_output(self) -> Dict[str, Any]:
        result_data = {"Copy numbers": {}, "Variants": self.result_vcfs}
        cns = result_data["Copy numbers"]
        for _, v in self.all_vars["cns"].items():
            cns[v["name"]] = 2 + v["cn_diff"]
        for _, v in self.all_vars["dels"].items():
            cns[v["name"]] = v["cn_diff"]

        # Load known conversions if database is configured
        known_conversions = []
        if self.conversion_db:
            known_conversions = GeneConversionDetector.load_known_conversions(
                self.conversion_db
            )

        if self.conversion_regions:
            result_data["gene_conversions"] = []
            for seg in self.converted_segments:
                is_fusion = getattr(seg, "is_fusion_candidate", False)
                conversion_data = {
                    "region": f"{self.chr}:{seg.start}-{seg.end}",
                    "converted_alleles": int(seg.conversion_allele_count),
                    "is_fusion_candidate": bool(is_fusion),
                }
                if is_fusion:
                    fusion_type = getattr(seg, "fusion_type", None)
                    if fusion_type:
                        conversion_data["fusion_type"] = fusion_type
                else:
                    conversion_data["conversion_sites"] = [
                        s.id for s in seg.signals if s.is_conversion_site()
                    ]

                # Match against known conversions if available
                if known_conversions:
                    matches = GeneConversionDetector.match_known_conversion(
                        seg, known_conversions, self.gene_names, match_threshold=0.8
                    )
                    if matches:
                        if not is_fusion:
                            conversion_data["interpretation"] = [
                                {
                                    "event_name": event_name,
                                    "match_score": score,
                                    "matched_variants": list(matched_ids),
                                }
                                for event_name, score, matched_ids in matches
                            ]
                        else:
                            conversion_data["interpretation"] = [
                                {
                                    "event_name": event_name,
                                }
                                for event_name, _, _ in matches
                            ]
                    else:
                        event_type = "fusion" if is_fusion else "conversion"
                        # conversion_data["interpretation"] = (
                        #     f"Novel {event_type} (no known match)"
                        # )
                else:
                    conversion_data["interpretation"] = (
                        "No conversion database configured"
                    )

                result_data["gene_conversions"].append(conversion_data)

        return result_data

    def solve_longdel(self) -> None:
        self.cn_model = CopyNumberModel(self.read_data, self, self.params)
        result = {}  # longdel_id -> cn_diff
        if not self.longdel_names:
            return
        for regions in self.segdup_regions:
            del_region = {}
            solved = []
            non_del = []
            for r in regions:
                if (
                    r.longdel_id > 0 and r.longdel_id & (r.longdel_id - 1) == 0
                ):  # single deletion id bit set
                    if r.longdel_id not in del_region:
                        del_region[r.longdel_id] = []
                    del_region[r.longdel_id].append(r.region)
                    if r.longdel_id not in solved:
                        solved.append(r.longdel_id)
                elif r.longdel_id == 0:
                    non_del.append(r.region)
            if not del_region:
                continue
            del_region_list = [",".join(del_region[id]) for id in solved]
            # Get corresponding longdel names for the solved deletion IDs
            solved_longdel_names = [self.longdel_names[id] for id in solved]
            longdels = self.cn_model.detect_longdels(
                ",".join(non_del),
                del_region_list,
                self.read_data["short_read"]["liftover"],
                solved_longdel_names,
            )
            for i, diff in enumerate(longdels):
                result[solved[i]] = diff
        if len(solved) != len(self.longdel_names):
            del_region = {}
            non_del = []
            for r in self.all_regions:
                if r.dup:
                    continue
                if (
                    r.longdel_id > 0 and r.longdel_id & (r.longdel_id - 1) == 0
                ):  # single deletion id bit set
                    if r.longdel_id not in del_region:
                        del_region[r.longdel_id] = []
                    del_region[r.longdel_id].append(r.region)
                    if r.longdel_id not in solved:
                        solved.append(r.longdel_id)
                elif r.longdel_id == 0:
                    non_del.append(r.region)
            if not del_region:
                return
            del_region_list = [",".join(del_region[id]) for id in solved]
            # Get corresponding longdel names for the solved deletion IDs
            solved_longdel_names = [self.longdel_names[id] for id in solved]
            longdels = self.cn_model.detect_longdels(
                ",".join(non_del),
                del_region_list,
                self.read_data["short_read"]["liftover"],
                solved_longdel_names,
            )
            for i, diff in enumerate(longdels):
                result[solved[i]] = diff
        self.logger.info("Result:")
        for id, cn_diff in result.items():
            self.logger.info(f" * {self.longdel_names[id]}: delta_CN = {cn_diff}")
            if not cn_diff:
                continue
            # remove del from del_regions
            self.all_regions = GeneRegion.remove_region_by_longdel(self.all_regions, id)
            for i, r in enumerate(self.segdup_regions):
                self.segdup_regions[i] = GeneRegion.remove_region_by_longdel(r, id)
            self.longdel_names.pop(id)

    def segmentation(self) -> None:
        # process each region. Get high mapq regions
        self.logger.info("Analyzing read stats of each gene")
        for name, rd in self.read_data.items():
            bam, params = rd["bam"], rd["params"]
            self.logger.info(f" * Running {name}...")
            segments = bam.get_segments(",".join(self.gene_regions))
            high_qual_segments = [s for s in segments if s.mq >= params["min_map_qual"]]
            high_mq_regions = ",".join(
                [
                    IntervalList.interval_str(s.chrom, s.start, s.end)
                    for s in high_qual_segments
                ]
            )
            high_mq_regions = (
                IntervalList(high_mq_regions)
                .intersect(",".join([r for r in self.liftover_target_regions if r]))
                .union(",".join([r for r in self.nonliftover_region if r]))
                .to_region_str()
            )
            self.high_mq_regions[name] = high_mq_regions

    def solve_cn(self) -> None:
        region_ids = list(self.all_vars["cns"].keys())
        del_ids = list(self.all_vars["dels"].keys())

        self.cn_model = CopyNumberModel(self.read_data, self, self.params)
        liftover_ads = self.cn_model.get_segdup_ads("short_read")
        for regions in self.liftover_segdup_regions:
            for r in regions:
                itv = IntervalList(r.region)
                r.ads = {
                    k: v for k, v in liftover_ads.items() if itv.get(self.chr, k, k + 1)
                }

        region_states = range(-2, 5)  # [-2, -1, 0, 1, 2, 3] -> 6 possibilities
        del_states = range(-2, 1)  # [-2, -1, 0] -> 3 possibilities

        # Initialize all variables to neutral state (0)
        for rid in region_ids:
            self.all_vars["cns"][rid]["cn_diff"] = 0
        for did in del_ids:
            self.all_vars["dels"][did]["cn_diff"] = 0

        # Cache for previously calculated states to avoid recomputation
        state_cache = {}

        def get_current_state() -> tuple[int, ...]:
            """Get current state as a hashable tuple"""
            cn_states = tuple(
                self.all_vars["cns"][rid]["cn_diff"] for rid in region_ids
            )
            del_states = tuple(self.all_vars["dels"][did]["cn_diff"] for did in del_ids)
            return cn_states + del_states

        def set_state(state: tuple[int, ...]) -> None:
            """Set the current state from a tuple"""
            for i, rid in enumerate(region_ids):
                self.all_vars["cns"][rid]["cn_diff"] = state[i]
            for j, did in enumerate(del_ids):
                self.all_vars["dels"][did]["cn_diff"] = state[len(region_ids) + j]

        def get_cached_prob(state: tuple[int, ...]) -> float:
            """Get probability from cache or calculate and cache it"""
            if state not in state_cache:
                set_state(state)
                prob = self.cal_logprob_cn()
                state_cache[state] = prob

                # Log the state being evaluated
                cn = ", ".join(
                    [
                        v["name"] + ": " + str(v["cn_diff"])
                        for v in self.all_vars["cns"].values()
                    ]
                )
                msg = f"config: CN: {cn}"
                if len(self.all_vars["dels"]) > 0:
                    d = ",".join(
                        [
                            v["name"] + ": " + str(v["cn_diff"])
                            for v in self.all_vars["dels"].values()
                        ]
                    )
                    msg += f", DEL: {d}"
                self.logger.debug(f"{msg}. log_prob: {prob}")
            return state_cache[state]

        # Hill climbing algorithm
        current_state = get_current_state()
        current_prob = get_cached_prob(current_state)

        # Log initial state
        cn = ", ".join(
            [
                v["name"] + ": " + str(v["cn_diff"])
                for v in self.all_vars["cns"].values()
            ]
        )
        msg = f"Initial state - CN: {cn}"
        if len(self.all_vars["dels"]) > 0:
            d = ", ".join(
                [
                    v["name"] + ": " + str(v["cn_diff"])
                    for v in self.all_vars["dels"].values()
                ]
            )
            msg += f", DEL: {d}"
        self.logger.info(f"{msg}. log_prob: {current_prob}")

        step = 0
        improved = True
        while improved:
            improved = False
            best_neighbor_state = None
            best_neighbor_prob = current_prob

            # Test all neighbors (±1 for each variable)
            # Single variable moves
            for i in range(len(region_ids)):
                for delta in [-1, 1]:
                    new_cn_diff = current_state[i] + delta
                    if new_cn_diff in region_states:
                        neighbor_state = list(current_state)
                        neighbor_state[i] = new_cn_diff
                        neighbor_state = tuple(neighbor_state)

                        neighbor_prob = get_cached_prob(neighbor_state)
                        if neighbor_prob > best_neighbor_prob:
                            best_neighbor_state = neighbor_state
                            best_neighbor_prob = neighbor_prob
                            improved = True

            for i in range(len(del_ids)):
                for delta in [-1, 1]:
                    del_idx = len(region_ids) + i
                    new_del_diff = current_state[del_idx] + delta
                    if new_del_diff in del_states:
                        neighbor_state = list(current_state)
                        neighbor_state[del_idx] = new_del_diff
                        neighbor_state = tuple(neighbor_state)

                        neighbor_prob = get_cached_prob(neighbor_state)
                        if neighbor_prob > best_neighbor_prob:
                            best_neighbor_state = neighbor_state
                            best_neighbor_prob = neighbor_prob
                            improved = True

            # Two variable moves - region pairs
            for i in range(len(region_ids)):
                for j in range(i + 1, len(region_ids)):
                    for delta_i in [-1, 1]:
                        for delta_j in [-1, 1]:
                            new_cn_diff_i = current_state[i] + delta_i
                            new_cn_diff_j = current_state[j] + delta_j
                            if (
                                new_cn_diff_i in region_states
                                and new_cn_diff_j in region_states
                            ):
                                neighbor_state = list(current_state)
                                neighbor_state[i] = new_cn_diff_i
                                neighbor_state[j] = new_cn_diff_j
                                neighbor_state = tuple(neighbor_state)

                                neighbor_prob = get_cached_prob(neighbor_state)
                                if neighbor_prob > best_neighbor_prob:
                                    best_neighbor_state = neighbor_state
                                    best_neighbor_prob = neighbor_prob
                                    improved = True

            # Two variable moves - deletion pairs
            for i in range(len(del_ids)):
                for j in range(i + 1, len(del_ids)):
                    for delta_i in [-1, 1]:
                        for delta_j in [-1, 1]:
                            del_idx_i = len(region_ids) + i
                            del_idx_j = len(region_ids) + j
                            new_del_diff_i = current_state[del_idx_i] + delta_i
                            new_del_diff_j = current_state[del_idx_j] + delta_j
                            if (
                                new_del_diff_i in del_states
                                and new_del_diff_j in del_states
                            ):
                                neighbor_state = list(current_state)
                                neighbor_state[del_idx_i] = new_del_diff_i
                                neighbor_state[del_idx_j] = new_del_diff_j
                                neighbor_state = tuple(neighbor_state)

                                neighbor_prob = get_cached_prob(neighbor_state)
                                if neighbor_prob > best_neighbor_prob:
                                    best_neighbor_state = neighbor_state
                                    best_neighbor_prob = neighbor_prob
                                    improved = True

            # Two variable moves - region and deletion pairs
            for i in range(len(region_ids)):
                for j in range(len(del_ids)):
                    for delta_i in [-1, 1]:
                        for delta_j in [-1, 1]:
                            del_idx = len(region_ids) + j
                            new_cn_diff = current_state[i] + delta_i
                            new_del_diff = current_state[del_idx] + delta_j
                            if (
                                new_cn_diff in region_states
                                and new_del_diff in del_states
                            ):
                                neighbor_state = list(current_state)
                                neighbor_state[i] = new_cn_diff
                                neighbor_state[del_idx] = new_del_diff
                                neighbor_state = tuple(neighbor_state)

                                neighbor_prob = get_cached_prob(neighbor_state)
                                if neighbor_prob > best_neighbor_prob:
                                    best_neighbor_state = neighbor_state
                                    best_neighbor_prob = neighbor_prob
                                    improved = True

            # Move to the best neighbor if improvement found
            if improved and best_neighbor_state is not None:
                step += 1
                current_state = best_neighbor_state
                current_prob = best_neighbor_prob
                set_state(current_state)

                # Log the optimization step
                cn = ", ".join(
                    [
                        v["name"] + ": " + str(v["cn_diff"])
                        for v in self.all_vars["cns"].values()
                    ]
                )
                msg = f"Step {step} - CN: {cn}"
                if len(self.all_vars["dels"]) > 0:
                    d = ", ".join(
                        [
                            v["name"] + ": " + str(v["cn_diff"])
                            for v in self.all_vars["dels"].values()
                        ]
                    )
                    msg += f", DEL: {d}"
                self.logger.info(f"{msg}. log_prob: {current_prob}")

        # Ensure final optimal state is set in self.all_vars
        set_state(current_state)
        self.logger.debug(f"Hill climbing converged to log_prob: {current_prob}")
        self.logger.debug(f"States evaluated: {len(state_cache)}")
        # set the called CN values for regions
        for r in self.all_regions:
            r.cn_diff = self.get_cn_diff(r.region_id, r.longdel_id)
            r.cn = 2 + r.cn_diff
        for regions in self.segdup_regions:
            for r in regions:
                r.cn_diff = self.get_cn_diff(r.region_id, r.longdel_id)
                r.cn = 2 + r.cn_diff
        for regions in self.liftover_segdup_regions:
            for r in regions:
                diffs = [
                    self.get_cn_diff(r.orig_region_id[j], r.orig_longdel_id[j])
                    for j in range(len(r.orig_region_id))
                ]
                r.cn_diff = sum(diffs)
                r.cn = 2 * len(self.liftover_segdup_regions) + r.cn_diff
        self.merged_regions = GeneRegion.merge_by_cn(self.all_regions)
        self.merged_segdup_regions = [
            GeneRegion.merge_by_cn(regions) for regions in self.segdup_regions
        ]
        self.merged_liftover_segdup_regions = [
            GeneRegion.merge_by_cn(regions) for regions in self.liftover_segdup_regions
        ]

    def get_cn_diff_by_region_id(self, region_id: int) -> int:
        cn_diff = 0
        if region_id == 0:
            if 0 in self.all_vars["cns"]:
                cn_diff = self.all_vars["cns"][0]["cn_diff"]
        else:
            i = 0
            while (1 << i) <= region_id:
                # Check if the bit at the current position is 'on' (i.e., 1)
                if (region_id >> i) & 1 and i in self.all_vars["cns"]:
                    cn_diff += self.all_vars["cns"][i]["cn_diff"]
                i += 1
        return cn_diff

    def get_prior_by_region_id(self, region_id: int) -> float:
        prior = 0.0
        if region_id == 0:
            if 0 in self.all_vars["cns"]:
                region_name = self.all_vars["cns"][0]["name"]
                cn_diff = self.all_vars["cns"][0]["cn_diff"]
                cn = cn_diff + 2
                prior = (
                    self._priors.get_log_prior(cn, region_name) if self._priors else 0.0
                )
        else:
            i = 0
            while (1 << i) <= region_id:
                # Check if the bit at the current position is 'on' (i.e., 1)
                if (region_id >> i) & 1 and i in self.all_vars["cns"]:
                    region_name = self.all_vars["cns"][i]["name"]
                    cn_diff = self.all_vars["cns"][i]["cn_diff"]
                    cn = cn_diff + 2
                    prior += (
                        self._priors.get_log_prior(cn, region_name)
                        if self._priors
                        else 0.0
                    )
                i += 1
        return prior

    def get_cn_diff_by_longdel_id(self, longdel_id: int) -> int:
        cn_diff = 0
        if longdel_id > 0:  # is a deletion
            i = 0
            while (1 << i) <= longdel_id:
                # Check if the bit at the current position is 'on' (i.e., 1)
                if (longdel_id >> i) & 1 and i in self.all_vars["dels"]:
                    cn_diff += self.all_vars["dels"][i]["cn_diff"]
                i += 1
        return cn_diff

    def get_cn_diff(self, region_id: int, longdel_id: int) -> int:
        cn_diff = 0
        cn_diff += self.get_cn_diff_by_region_id(region_id)
        cn_diff += self.get_cn_diff_by_longdel_id(longdel_id)
        return cn_diff

    def cal_logprob_cn(self) -> float:
        bam = self.read_data["short_read"]["bam"]
        liftover_bam = self.read_data["short_read"]["liftover"]
        high_mq_regions = IntervalList(self.high_mq_regions["short_read"])
        total_prob = 0.0
        min_region_size = 400
        # check validity
        for r in self.all_regions:
            cn_diff = self.get_cn_diff(r.region_id, r.longdel_id)
            if cn_diff < -2:
                return -1e9
        for r in self.all_regions:
            mq_region = high_mq_regions.intersect(r.region)
            if mq_region.size() < min_region_size:
                continue
            cn_diff = self.get_cn_diff(r.region_id, r.longdel_id)
            cn = cn_diff + 2

            # Calculate CN prior (for the final copy number)
            # Use region-specific priors if available (via region name)
            cn_log_prior = self.get_prior_by_region_id(r.region_id)

            # Calculate longdel priors for deletions affecting this region
            longdel_log_prior = 0.0
            if r.longdel_id > 0:
                # Check each deletion bit that's set
                for del_id in range(len(self.longdel_names)):
                    if (
                        r.longdel_id >> del_id
                    ) & 1:  # This deletion affects this region
                        del_name = self.longdel_names[del_id]
                        cn_diff_for_del = (
                            self.all_vars["dels"].get(del_id, {}).get("cn_diff", 0)
                        )

                        if del_name in self._longdel_log_priors:
                            if cn_diff_for_del != 0:  # Deletion is present
                                longdel_log_prior += self._longdel_log_priors[del_name]
                            else:  # Deletion is absent
                                # Use log(1 - p) for no deletion
                                p_del = np.exp(self._longdel_log_priors[del_name])
                                longdel_log_prior += (
                                    np.log(1 - p_del) if p_del < 0.999 else -1e6
                                )

            total_log_prior = cn_log_prior + longdel_log_prior
            region_prob = bam.calc_logprob(
                mq_region.to_region_str(), cn, log_prior=total_log_prior
            )
            total_prob += region_prob
        cn_neutral = 2 * len(self.liftover_segdup_regions)
        for regions in self.liftover_segdup_regions:
            for r in regions:
                diffs = [
                    self.get_cn_diff(r.orig_region_id[j], r.orig_longdel_id[j])
                    for j in range(len(r.orig_region_id))
                ]

                # Calculate CN priors and longdel priors for each paralogous region
                cn_log_priors = []
                longdel_log_priors = []

                for j in range(len(r.orig_region_id)):
                    diff = diffs[j]
                    # CN prior
                    cn = diff + 2
                    cn_log_priors.append(
                        self.get_prior_by_region_id(r.orig_region_id[j])
                    )
                    # Longdel priors for this paralogous region
                    longdel_prior = 0.0
                    longdel_id = r.orig_longdel_id[j]
                    if longdel_id > 0:
                        for del_id in range(len(self.longdel_names)):
                            if (
                                longdel_id >> del_id
                            ) & 1:  # This deletion affects this region
                                del_name = self.longdel_names[del_id]
                                cn_diff_for_del = (
                                    self.all_vars["dels"]
                                    .get(del_id, {})
                                    .get("cn_diff", 0)
                                )

                                if del_name in self._longdel_log_priors:
                                    if cn_diff_for_del != 0:  # Deletion is present
                                        longdel_prior += self._longdel_log_priors[
                                            del_name
                                        ]
                                    else:  # Deletion is absent
                                        p_del = np.exp(
                                            self._longdel_log_priors[del_name]
                                        )
                                        longdel_prior += (
                                            np.log(1 - p_del) if p_del < 0.999 else -1e6
                                        )

                    longdel_log_priors.append(longdel_prior)

                cn_total = cn_neutral + sum(diffs)
                log_prior = sum(cn_log_priors) + sum(longdel_log_priors)
                diff0 = diffs[0]
                ratio = (2 + diff0) / cn_total if cn_total else 0
                liftover_prob = liftover_bam.calc_logprob(
                    r.region, cn_total, r.ads, ratio, log_prior=log_prior
                )
                total_prob += liftover_prob
        return total_prob

    def call_cn(self, read_data: Dict[str, Any], params: Dict[str, Any]) -> None:
        self.read_data = read_data
        self.params = params
        # self.solve_longdel()
        self.segmentation()
        self.solve_cn()

    def _call_variant(self, params: dict) -> bool:
        """Helper method to call variants for a single region."""
        bam = params["bam"]
        region = params["region"]
        out_vcf = params["out_vcf"]
        phased_vcf = params["phased_vcf"]
        seq_key = params["seq_key"]
        bam_type = params["bam_type"]

        try:
            # Call variants and phase
            bam.call_variant(out_vcf, params)
            bam.phase_vcf(out_vcf, phased_vcf, params["ploidy"])

            self.logger.debug(
                f"Completed variant calling for {bam_type} {seq_key} bam in region {region}"
            )

            return True
        except Exception as e:
            self.logger.error(
                f"Error calling variants for {bam_type} {seq_key} bam in region {region}: {e}"
            )
            return False

    def call_small_var(self) -> None:
        seq_keys = ["short_read"]
        if "long_read" in self.read_data:
            seq_keys.append("long_read")
        tasks = []
        for seq_key in seq_keys:
            bam = self.read_data[seq_key]["bam"]
            orig_params = self.read_data[seq_key]["params"]
            orig_params["seq_key"] = seq_key
            orig_params["bam_type"] = "orig"
            for i, r in enumerate(self.merged_regions):
                if not r.cn:
                    continue
                params = copy.deepcopy(orig_params)
                params["bam"] = bam
                params["region"] = r.region
                params["ploidy"] = r.cn
                if len(self.diff_vcf) == 1:
                    params["dbsnp"] = self.diff_vcf[0]
                else:
                    if r.dup == 0:
                        params.pop("dbsnp", None)
                    else:
                        params["dbsnp"] = self.diff_vcf[r.dup - 1]
                params["out_vcf"] = f"{params['prefix']}.{r.name}.{i}.raw.vcf.gz"
                params["phased_vcf"] = f"{params['prefix']}.{r.name}.{i}.phased.vcf.gz"
                r.called_vcf[seq_key] = {"orig": params["out_vcf"]}
                r.phased_vcf[seq_key] = {"orig": params["phased_vcf"]}
                tasks.append(params)
            bam = self.read_data[seq_key]["liftover"]
            orig_params["bam_type"] = "liftover"
            for i, regions in enumerate(self.merged_liftover_segdup_regions):
                for j, r in enumerate(regions):
                    if not r.cn:
                        continue
                    params = copy.deepcopy(orig_params)
                    params["bam"] = bam
                    params["region"] = r.region
                    params["ploidy"] = r.cn
                    if len(self.diff_vcf) == 1:
                        params["dbsnp"] = self.diff_vcf[0]
                    else:
                        params["dbsnp"] = self.diff_vcf[i]
                    params["out_vcf"] = (
                        f"{params['prefix']}.{r.name}.liftover.{i}_{j}.raw.vcf.gz"
                    )
                    params["phased_vcf"] = (
                        f"{params['prefix']}.{r.name}.liftover.{i}_{j}.phased.vcf.gz"
                    )
                    r.called_vcf[seq_key] = {"liftover": params["out_vcf"]}
                    r.phased_vcf[seq_key] = {"liftover": params["phased_vcf"]}
                    tasks.append(params)
        if not tasks:
            self.logger.info(
                "No regions with copy number > 0, skipping variant calling"
            )
            return

        max_threads = getattr(self, "variant_call_threads", 4)  # Default to 4 threads
        num_threads = min(max_threads, len(tasks))

        self.logger.info(
            f"Starting parallel variant calling with {num_threads} threads for {len(tasks)} tasks"
        )

        try:
            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                # Submit region variant calling tasks
                result_future = []
                for task in tasks:
                    future = executor.submit(self._call_variant, task)
                    result_future.append(future)

                # check task results
                for future in result_future:
                    result = future.result()
                    if not result:
                        raise RuntimeError("One of the variant calling tasks failed.")

        except Exception as e:
            self.logger.error(f"Parallel variant calling failed: {e}")
            raise

        self.logger.info("Completed parallel variant calling")

    def resolve_phased(self, params: Dict[str, Any]) -> None:
        # first merge all liftover_segdup regions into original regions
        for i, segdup_regions in enumerate(self.orig_segdup_regions):
            for r in segdup_regions:
                self.merged_liftover_segdup_regions[i] = GeneRegion.insert(
                    self.merged_liftover_segdup_regions[i], r, use_orig=True
                )
        for r in self.merged_regions:
            r.total_cn = r.cn
            r.dup = 0
        for i in range(len(self.merged_liftover_segdup_regions[0])):
            all_matched = []
            for j in range(len(self.merged_liftover_segdup_regions)):
                region = self.merged_liftover_segdup_regions[j][i]
                self.merged_regions = GeneRegion.insert(
                    self.merged_regions, region, use_orig=True
                )
                # matches are references to elements of self.merged_regions (linear).
                # this opereration may break one or more of self.merged_regions
                # to match the boundary of region
                matches = GeneRegion.get_regions_by_coord(
                    self.merged_regions, region.region
                )
                # modification here is done on elements of self.merged_segments (linear)
                for matched in matches:
                    # matched = matches[0]
                    matched.total_cn = region.cn
                    matched.total_cn_diff = region.cn_diff
                    matched.dup = j + 1
                all_matched.append(matches)
            max_matched = max([len(m) for m in all_matched])
            min_matched = min([len(m) for m in all_matched])
            # Validate that all liftover segments have consistent matching
            assert max_matched == min_matched, "Not all liftover segments are matched."
            for j, matches in enumerate(all_matched):
                for matched in matches:
                    matched.matched_region = []
                    for jj, mm in enumerate(all_matched):
                        if jj == j:
                            continue
                        for m in mm:
                            lifted = copy.deepcopy(m)
                            lifted.liftover(self.mapping)
                            if (
                                lifted.region == matched.region
                                or IntervalList(matched.region)
                                .subtract(lifted.region)
                                .size()
                                < 20
                            ):
                                matched.matched_region.append(m)

        # Now all regions are CN-distinct with full info about its CN and its matching segdup
        # For each region, resolve phased small variants
        all_vars = []
        vcf_proc = Phased_vcf(self, self.read_data)
        seq_keys = ["short_read"]
        if "long_read" in self.read_data:
            seq_keys.append("long_read")
        for r in self.merged_regions:
            if not r.cn:
                continue
            region_itv = IntervalList(r.region)
            resolved = {}
            if r.dup and r.total_cn - r.cn > 0:
                cn1 = r.cn
                cn2 = r.total_cn - r.cn
                for seq_key in seq_keys:
                    if "liftover" not in r.phased_vcf[seq_key]:
                        orig_fname = r.phased_vcf[seq_key].get("orig")
                        if orig_fname:
                            orig_vcf = vcflib.VCF(orig_fname)
                            orig_vars = [
                                v for v in orig_vcf if (v.chrom, v.pos) in region_itv and 'MLrejected' not in v.filter
                            ]
                            orig_vcf.close()
                        else:
                            orig_vars = []
                        lift_vars = []  # No liftover variants when no liftover file
                        resolved[seq_key] = lift_vars, orig_vars
                        continue
                    run_params = {
                        "region": r.region,
                        "orig_region": r.region,
                        "gene_id": r.dup - 1,
                        "seq": seq_key,
                        "cns": (cn1, cn2),
                        "cn_diff": self.get_cn_diff_by_longdel_id(r.longdel_id),
                        "mismatch": params.get("max_mismatch", vcf_proc.max_mismatch),
                    }
                    # right now only consider resolving two segup regions
                    try:
                        matched_vcf = r.matched_region[0].phased_vcf[seq_key]["orig"]
                    except (KeyError, IndexError, AttributeError):
                        matched_vcf = None
                    resolved[seq_key] = vcf_proc.assign(
                        r.phased_vcf[seq_key]["liftover"],
                        r.phased_vcf[seq_key].get("orig"),
                        matched_vcf,
                        run_params,
                    )
            else:
                for seq_key in seq_keys:
                    lift_fname = r.phased_vcf[seq_key].get("liftover")
                    if lift_fname:
                        lift_vcf = vcflib.VCF(lift_fname)
                        lift_vars = [
                            v for v in lift_vcf if (v.chrom, v.pos) in region_itv and 'MLrejected' not in v.filter
                        ]
                        lift_vcf.close()
                    else:
                        lift_vars = []
                    orig_fname = r.phased_vcf[seq_key].get("orig")
                    if orig_fname:
                        orig_vcf = vcflib.VCF(orig_fname)
                        orig_vars = [
                            v for v in orig_vcf if (v.chrom, v.pos) in region_itv and 'MLrejected' not in v.filter
                        ]
                        orig_vcf.close()
                    else:
                        orig_vars = []
                    resolved[seq_key] = lift_vars, orig_vars
            # resolve short and long read
            if len(resolved) > 1:
                short_liftover = iter(resolved["short_read"][0])
                short_vcf = iter(resolved["short_read"][1])
                long_vcf = iter(resolved["long_read"][1])
                out_vcf1 = []
                for v1 in vcf_proc.vcf_grouper(short_liftover, short_vcf, long_vcf):
                    out_vcf1.append(v1)
            else:
                out_vcf1 = resolved["short_read"][1]
            for v in out_vcf1:
                heapq.heappush(all_vars, (v.pos, f"{v.ref}_{','.join(v.alt)}", v))
        # output all variants
        out_vcf_fname = f"{params['outdir']}/{params['sample_name']}.{self.gene_names[0]}.result.vcf.gz"
        out_vcf = vcflib.VCF(out_vcf_fname, "w")
        diff_vcf = vcflib.VCF(self.merged_regions[0].phased_vcf["short_read"]["orig"])
        out_vcf.copy_header(
            diff_vcf,
            update='##ALT=<ID=DEL,Description="Deletion relative to the reference">',
        )
        out_vcf.emit_header()
        while all_vars:
            var = heapq.heappop(all_vars)
            out_vcf.emit(var[2])
        out_vcf.close()
        self.result_vcfs = out_vcf_fname

    def detect_gene_conversions(self) -> None:
        """Detect gene conversions across all configured conversion regions."""
        if not self.conversion_regions:
            return

        self.logger.info(
            f"Performing gene conversion detection on {len(self.conversion_regions)} region(s)"
        )

        all_segments = []
        for conversion_config in self.conversion_regions:
            # Get conversion region coordinates
            chrom, region_start, region_end = (
                conversion_config["region"].split(":")[0],
                int(conversion_config["region"].split(":")[1].split("-")[0]),
                int(conversion_config["region"].split(":")[1].split("-")[1]),
            )

            # Find all overlapping duplicated regions and group by copy number
            # Split conversion detection if CN changes within the conversion region
            segments_by_cn = []  # List of (cns, start, end, vcf_files)
            current_cns = None
            current_start = region_start
            current_vcfs = []

            for r in self.merged_regions:
                if not r.dup:
                    continue

                # Check if region overlaps with conversion region
                r_chrom, r_start, r_end = r.get_region()
                if r_chrom != chrom:
                    continue
                if r_end <= region_start or r_start >= region_end:
                    continue

                # Get copy numbers for this region: [ref_count, alt_count]
                region_cns = [r.cn, r.total_cn - r.cn]

                # Check if CN changed - need to split detection
                if current_cns is None:
                    # First overlapping region
                    current_cns = region_cns
                    current_vcfs = [r.phased_vcf["short_read"]["liftover"]]
                elif region_cns == current_cns:
                    # Same CN, continue accumulating
                    current_vcfs.append(r.phased_vcf["short_read"]["liftover"])
                else:
                    # CN changed - save current segment and start new one
                    segments_by_cn.append(
                        (current_cns, current_start, r_start, current_vcfs)
                    )
                    current_cns = region_cns
                    current_start = r_start
                    current_vcfs = [r.called_vcf["short_read"]["liftover"]]

            # Add final segment
            if current_cns is not None and current_vcfs:
                segments_by_cn.append(
                    (current_cns, current_start, region_end, current_vcfs)
                )

            # Run conversion detection for each CN-uniform segment
            for cns, seg_start, seg_end, vcf_files in segments_by_cn:
                cfg = copy.deepcopy(conversion_config)
                cfg["conversion_region"] = f"{chrom}:{seg_start}-{seg_end}"
                cfg["copy_numbers"] = cns
                cfg["gene_strand"] = self.gene_strand

                detector = GeneConversionDetector(self, cfg)
                detector.load_variants(set(vcf_files))
                segments = detector.detect_conversions()
                all_segments.extend(segments)

        self.logger.info(
            f"Gene conversion detection complete: found {len(all_segments)} conversion segment(s)"
        )
        for i, seg in enumerate(all_segments, 1):
            self.logger.info(
                f"  Conversion {i}: {seg.start}-{seg.end} ({seg.end - seg.start} bp), "
                f"{seg.num_sites} sites, converted {seg.conversion_allele_count} alleles"
            )

        self.converted_segments = all_segments
