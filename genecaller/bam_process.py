import vcflib
import pysam
import os
import subprocess
import sys
from .util import IntervalList, get_data_file
from .logging import get_logger
import heapq
from collections import namedtuple, Counter
from scipy.stats import binom, gamma
import tempfile
import numpy as np
import pandas as pd
import re
import copy
import bisect
import multiprocessing
from typing import List, Dict, Tuple, Optional, Any, Union, Callable

Segment = namedtuple("Segment", "chrom, start, end, cnt, mq, mean, std, mean0, std0")
DepthSegment = namedtuple("DepthSegment", "chrom, start, end, mq, mean, mean0")


class Bam:
    def __init__(self, bam: str, ref: str, param: Dict[str, Any]) -> None:
        self.bam = bam
        self.bamh = None
        self.header = None
        self.long_read = False
        self.ref = ref
        self.depths = {}
        self.segments = {}
        self._depth_coords_cache = {}  # Cache for sorted coordinate lists
        self.cnv_model = ""
        self.lr_model = ""
        self.sr_model = ""
        self.prefix = None
        self.sample_name = "SAMPLE"
        self.use_existing = True
        self.clipped_bam = None
        self.outdir = None
        self.tmpdir = param["tmpdir"]
        self.dp_norm = 1.0
        if "tools" not in param:
            raise Exception("External tools are not set up")
        for k, v in param["tools"].items():
            setattr(self, k, v)
        for k, v in param.items():
            if hasattr(self, k):
                setattr(self, k, v)
        if (
            self.cnv_model
            and not os.path.exists(self.cnv_model)
            and not os.path.exists(os.path.dirname(self.cnv_model))
        ):
            self.cnv_model = get_data_file(self.cnv_model)
        if (
            self.lr_model
            and not os.path.exists(self.lr_model)
            and not os.path.exists(os.path.dirname(self.lr_model))
        ):
            self.lr_model = get_data_file(self.lr_model)
        if (
            self.sr_model
            and not os.path.exists(self.sr_model)
            and not os.path.exists(os.path.dirname(self.sr_model))
        ):
            self.sr_model = get_data_file(self.sr_model)
        self.threads = multiprocessing.cpu_count()
        self.logger = get_logger(self.__class__.__name__)

    def reset(self) -> None:
        self.depths = {}
        self.segments = {}
        self._depth_coords_cache.clear()  # Clear cache when resetting
        self.clipped_bam = None
        self.dp_norm = 1.0

    @staticmethod
    def rc(seq: str) -> str:
        rev = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
        return "".join([rev[s] for s in seq[-1::-1]])

    def call_variant(self, output: str, param: Dict[str, Any]) -> None:
        if self.use_existing and os.path.exists(output) and os.path.getsize(output):
            return
        if "ploidy" in param and param["ploidy"] != 2:
            param["algo"] = "Haplotyper"
        self.call_dnascope(output, param)

    def call_dnascope(self, output: str, param: Dict[str, Any]) -> None:
        if self.use_existing and os.path.exists(output) and os.path.getsize(output):
            return
        bam = self.clipped_bam if self.clipped_bam else self.bam
        driver_opt = param.get("driver_opt", "") + f" -i {bam} -r {self.ref}"
        algo_opt = param.get("algo_opt", "")
        apply_model = ""
        algo = "DNAscope" if "algo" not in param else param["algo"]
        ploidy = param.get("ploidy", 2)
        algo_opt += f" --ploidy {ploidy}"
        given = param.get("given", "")
        if given and os.path.exists(given):
            algo_opt += f" --given {given} -d {given}"
        if "model" in param and algo == "DNAscope" and ploidy == 2 and not given:
            algo_opt += f" --model {param['model']}"
            apply_model = param["model"]
        if "region" in param:
            driver_opt += f" --interval {param['region']}"
        elif "dbsnp" in param:
            algo_opt += f" -d {param['dbsnp']}"
        if apply_model:
            tmp_output = output.replace("vcf", "tmp.vcf")
        else:
            tmp_output = output
        cmd = (
            f"{self.sentieon} driver {driver_opt} --algo {algo} {algo_opt} {tmp_output}"
        )
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if result.returncode:
            raise Exception(result.stderr)
        if apply_model:
            cmd = f"{self.sentieon} driver -r {self.ref} --algo DNAModelApply -v {tmp_output} --model {apply_model} {output}"
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.returncode:
                raise Exception(result.stderr)

    def call_genotyper(self, output: str, param: Dict[str, Any]) -> None:
        if self.use_existing and os.path.exists(output) and os.path.getsize(output):
            return
        bam = self.clipped_bam if self.clipped_bam else self.bam
        driver_opt = param.get("driver_opt", "") + f" -i {bam} -r {self.ref}"
        algo_opt = param.get("algo_opt", "")
        if "region" in param:
            driver_opt += f" --interval {param['region']}"
        if "ploidy" in param:
            algo_opt += f" --ploidy {param['ploidy']}"
        if "given" in param:
            algo_opt += f" --given {param['given']} -d {param['given']}"
        elif "dbsnp" in param:
            algo_opt += f" -d {param['dbsnp']}"
        cmd = (
            f"{self.sentieon} driver {driver_opt} --algo Genotyper {algo_opt} {output}"
        )
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if result.returncode:
            raise Exception(result.stderr)

    def call_cnvscope(self, output: str, param: Dict[str, Any]) -> None:
        if self.use_existing and os.path.exists(output) and os.path.getsize(output):
            return
        driver_opt = param.get("driver_opt", "") + f" -i {self.bam} -r {self.ref}"
        algo_opt = param.get("algo_opt", "")
        if "region" in param:
            driver_opt += f" --interval {param['region']}"
        if "output_norm" in param:
            algo_opt += f" --output_normalized {param['output_norm']}"
        if "model" in param:
            algo_opt += f" --model {param['model']}"
        else:
            raise Exception("Model is not specified for CNVscope.")
        cmd = f"{self.sentieon} driver {driver_opt} --algo CNVscope {algo_opt} {output}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def phase_vcf(self, in_vcf: str, out_vcf: str, ploidy: int = 2) -> None:
        if self.use_existing and os.path.exists(out_vcf) and os.path.getsize(out_vcf):
            return
        bam = self.clipped_bam if self.clipped_bam else self.bam
        if ploidy == 2:
            cmd = f"{self.whatshap} phase -o {out_vcf} {in_vcf} {bam} --reference {self.ref}"
        else:
            cmd = f"{self.whatshap} polyphase -o {out_vcf} {in_vcf} {bam} -p {ploidy} --reference {self.ref} -B 5 -t {self.threads}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()
        cmd = f"{self.sentieon} util vcfindex {out_vcf}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def get_depth(self, regions: Union[str, IntervalList]) -> List[Any]:
        if not regions:
            return []
        if isinstance(regions, str):
            regions = IntervalList(region=regions)
        chrs = set(regions.regions.keys())
        if not self.depths:
            param = {}
            if self.cnv_model:
                param["model"] = self.cnv_model
            param["region"] = regions.to_region_str(chrs)
            output = os.path.join(f"{self.prefix}.cnv.vcf.gz")
            self.call_depth(output, param)
        depths = []
        for chr in regions.regions:
            if chr not in self.depths:
                continue

            chr_regions = regions.regions[chr]
            chr_depths = self.depths[chr]

            # Use cached coordinate lists
            if chr not in self._depth_coords_cache:
                # Fallback: compute if not cached
                starts = [d.start for d in chr_depths]
                ends = [d.end for d in chr_depths]
                self._depth_coords_cache[chr] = (starts, ends)
            else:
                starts, ends = self._depth_coords_cache[chr]

            i = 0
            while i < len(chr_regions):
                start, end = chr_regions[i], chr_regions[i + 1]

                # Binary search using pre-computed lists
                left_idx = bisect.bisect_left(ends, start)
                right_idx = bisect.bisect_right(starts, end)

                # Check only segments in this range
                for j in range(left_idx, min(right_idx, len(chr_depths))):
                    d = chr_depths[j]
                    if d.start < end and d.end > start:  # Overlap check
                        depths.append(d)

                i += 2
        return depths

    def get_segments(self, regions: Union[str, IntervalList]) -> Any:
        if not regions:
            return []
        if isinstance(regions, str):
            regions = IntervalList(region=regions)
        chrs = set(regions.regions.keys())
        if not self.segments:
            param = {}
            if self.cnv_model:
                param["model"] = self.cnv_model
            param["region"] = regions.to_region_str(chrs)
            output = f"{self.prefix}.cnv.{'.'.join(chrs)}.vcf.gz"
            self.call_depth(output, param)
        segments = set()
        for chr, reg in regions.regions.items():
            if chr not in self.segments:
                continue
            chr_segments = self.segments[chr]
            i = 0
            while i < len(reg):
                s, e = reg[i], reg[i + 1]
                segments.update(
                    seg for seg in chr_segments if seg.start < e and seg.end > s
                )
                i += 2
        return segments

    def call_depth(self, output: str, param: Dict[str, Any]) -> None:
        orig_region = None
        if "region" in param:
            orig_region = param["region"]
            chrs = set()
            for parts in param["region"].split(","):
                chrs.add(parts.split(":")[0])
            param["region"] = ",".join(chrs)
        prefix = output[: output.index(".vcf")] if ".vcf" in output else output
        output_norm = param.get("output_norm", f"{prefix}.norm.tsv")
        if orig_region:
            tmp_output = prefix + ".tmp" + output[output.index(".vcf") :]
            tmp_output_norm = prefix + ".tmp.norm.tsv"
        else:
            tmp_output = output
            tmp_output_norm = output_norm
        param["output_norm"] = tmp_output_norm
        self.call_cnvscope(tmp_output, param)
        if orig_region:
            param["region"] = orig_region
        vcf = vcflib.VCF(tmp_output)
        depth_fields = ["contig", "start", "end", "MQ", "DP", "DP0"]
        for v in vcf:
            dp = v.samples[0]["DPS"]
            dp0 = v.samples[0]["DP0S"]
            self.segments.setdefault(v.chrom, []).append(
                Segment(
                    v.chrom,
                    v.pos,
                    v.end,
                    v.info["CNT"],
                    v.samples[0]["MQS"],
                    dp[0] / self.dp_norm,
                    dp[1] / self.dp_norm,
                    dp0[0] / self.dp_norm,
                    dp0[0] / self.dp_norm,
                )
            )
        depth_df = pd.read_csv(tmp_output_norm, sep="\t")[depth_fields]
        depth_df["DP"] /= self.dp_norm
        depth_df["DP0"] /= self.dp_norm
        high_mq_depth_df = depth_df[depth_df["MQ"] > 50]
        self.overall_depth_std = high_mq_depth_df["DP"].std()
        # Optimize DepthSegment creation using vectorized operations
        self.depths = {}
        self._depth_coords_cache = {}
        for chr, group in depth_df.groupby("contig"):
            # Use vectorized numpy operations instead of itertuples
            values = group.values
            chr_depths = [DepthSegment(*row) for row in values]
            self.depths[chr] = chr_depths
            # Cache coordinates directly from numpy arrays for better performance
            self._depth_coords_cache[chr] = (
                values[:, 1].tolist(),  # start column
                values[:, 2].tolist(),  # end column
            )
        if orig_region:
            orig_region_itv = IntervalList(region=orig_region)

            # Vectorized filtering using the same logic as the earlier optimization
            mask = pd.Series(False, index=depth_df.index)

            for chr in orig_region_itv.regions:
                chr_mask = depth_df["contig"] == chr
                if not chr_mask.any():
                    continue

                chr_regions = orig_region_itv.regions[chr]
                chr_rows = depth_df[chr_mask]

                i = 0
                while i < len(chr_regions):
                    start, end = chr_regions[i], chr_regions[i + 1]

                    # Vectorized overlap check
                    overlap_mask = (chr_rows["start"] < end) & (chr_rows["end"] > start)
                    mask[chr_rows[overlap_mask].index] = True

                    i += 2

            trim_depth_df = depth_df[mask]
            trim_depth_df.to_csv(output_norm, sep="\t", index=None)
            self.trim_to_region_vcf(orig_region, tmp_output, output)

    def trim_to_region_vcf(self, region: str, in_vcf: str, out_vcf: str) -> None:
        if region.endswith(".bed"):
            region_opt = f"-T {region}"
        else:
            region_opt = f"-r {region}"
        cmd = f"{self.bcftools} view {region_opt} {in_vcf} | {self.sentieon} util vcfconvert - {out_vcf}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def init_bamh(self) -> pysam.AlignmentFile:
        if self.bamh:
            self.bamh.close()
            self.bamh = None
        if self.bam.endswith(".cram"):
            if not self.ref:
                self.logger.error("No reference file provided for CRAM input")
                sys.exit(1)
            self.bamh = pysam.AlignmentFile(self.bam, "rc", reference_filename=self.ref)
        else:
            self.bamh = pysam.AlignmentFile(self.bam, "rb")
        return self.bamh

    def _pileup_safe(
        self, bamh: pysam.AlignmentFile, region: str, truncate: bool = True
    ) -> Any:
        """Safe pileup that handles CRAM files without multiple_iterators warning"""
        if self.bam.endswith(".cram"):
            # CRAM doesn't support multiple_iterators, so don't use it
            return bamh.pileup(
                region=region, truncate=truncate, multiple_iterators=False
            )
        else:
            return bamh.pileup(region=region, truncate=truncate)

    def init_outbamh(self, out_bam: str) -> pysam.AlignmentFile:
        self.init_bamh()
        if out_bam.endswith(".cram"):
            out_bamh = pysam.AlignmentFile(
                out_bam, "wc", reference_filename=self.ref, template=self.bamh
            )
        else:
            out_bamh = pysam.AlignmentFile(out_bam, "wb", template=self.bamh)
        return out_bamh

    def get_header_field(self, field: str) -> Any:
        if not self.header:
            self.init_bamh()
            self.header = self.bamh.header
        return self.header.get(field)

    def clip_to_region(
        self, out_bam: str, regions: str, min_length: Optional[int] = 300
    ) -> None:
        if self.use_existing and os.path.exists(out_bam) and os.path.getsize(out_bam):
            return
        if self.bamh is None:
            self.init_bamh()
        out_bamh = self.init_outbamh(out_bam)
        for region in regions.split(","):
            _, _, start, end = self.bamh.parse_region(region=region)
            for read in self.bamh.fetch(region=region):
                if not read.query_sequence:
                    continue
                ref_positions = read.get_reference_positions(full_length=True)
                qual = read.query_qualities
                if read.reference_start < start:
                    for idx, pos in enumerate(ref_positions):
                        if pos and pos >= start:
                            break
                    read.query_sequence = read.query_sequence[idx:]
                    ref_positions = ref_positions[idx:]
                    qual = qual[idx:]
                    read.query_qualities = qual
                    cigar = []
                    cur_idx = 0
                    for op, length in read.cigartuples:
                        if cur_idx >= idx:
                            cigar.append((op, length))
                        elif cur_idx + length > idx:
                            cigar.append((op, length - (idx - cur_idx)))
                        if op not in (2, 5):
                            cur_idx += length
                    read.cigartuples = cigar
                    read.reference_start = ref_positions[0]
                if read.reference_end > end:
                    for idx, pos in enumerate(ref_positions[::-1]):
                        if pos and pos <= end:
                            break
                    read.query_sequence = read.query_sequence[:-idx]
                    ref_positions = ref_positions[:-idx]
                    read.query_qualities = qual[:-idx]
                    cigar = []
                    cur_idx = 0
                    for op, length in read.cigartuples[::-1]:
                        if cur_idx >= idx:
                            cigar.append((op, length))
                        elif cur_idx + length > idx:
                            cigar.append((op, length - (idx - cur_idx)))
                        if op not in (2, 5):
                            cur_idx += length
                    read.cigartuples = cigar[::-1]
                if (
                    min_length is None
                    or read.query_sequence
                    and len(read.query_sequence) > min_length
                ):
                    out_bamh.write(read)
        self.bamh.close()
        self.bamh = None
        out_bamh.close()
        pysam.index(out_bam)

    def make_fasta(self, ref_file: str, region: str) -> list[str]:
        make_ref_cmd = (
            f"{self.samtools} faidx {self.ref} {region.replace(',', ' ')} > {ref_file}"
        )
        result = subprocess.run(
            make_ref_cmd, capture_output=True, text=True, shell=True
        )
        result.check_returncode()
        pysam.faidx(ref_file)
        contig_names = []
        with open(ref_file) as f:
            for line in f:
                if line.startswith(">"):
                    contig_names.append(line[1:].strip())
        return contig_names

    @staticmethod
    def to_region(
        positions: Dict[str, Dict[int, int]], min_gap: int = 1, min_size: int = 10
    ) -> List[str]:
        regions = []
        start = 0
        for chr, pos in positions.items():
            last_pos = -1
            last_region = None
            sorted_pos = sorted(pos.keys())
            for p in sorted_pos:
                if last_pos == -1:
                    last_region = f"{chr}:{p}-"
                    start = p
                    last_pos = p
                    continue
                if (
                    last_region
                    and p > last_pos + min_gap
                    and last_pos - start + 1 >= min_size
                ):
                    last_region += str(last_pos + 1)
                    regions.append(last_region)
                    last_region = f"{chr}:{p}-"
                    start = p
                last_pos = p
            if last_pos != -1 and last_region and last_pos - start + 1 >= min_size:
                last_region += str(last_pos + 1)
                regions.append(last_region)
        return regions

    def liftover_combine(self, gene: Any, out_bam: str) -> List[Dict[str, Any]]:
        out_prefix = out_bam[:-5] if out_bam.endswith(".cram") else out_bam[:-4]
        liftover_bams = []
        results = []
        for i in range(len(gene.liftover_regions)):
            liftover_bams.append(f"{out_prefix}.{gene.liftover_region_names[i]}.bam")
            mapping_dict = gene.mapping.g2tog1 if i == 0 else gene.mapping.g1tog2
            results.append(
                self.liftover(
                    gene.liftover_regions[i],
                    gene.realign_regions[i],
                    liftover_bams[-1],
                    mapping_dict,
                    self.long_read,
                )
            )
        if self.use_existing and os.path.exists(out_bam) and os.path.getsize(out_bam):
            return results
        cmd = f"{self.sentieon} driver -r {self.ref} -i {self.bam} -i {' -i '.join(liftover_bams)} --interval {gene.chr} --algo ReadWriter {out_bam}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()
        return results

    def liftover(
        self,
        extract: str,
        target: str,
        out_bam: str,
        mapping: Dict[int, int],
        crop_first: bool = False,
        min_ratio: float = 0.3,
    ) -> Dict[str, str]:
        out_prefix = out_bam[:-5] if out_bam.endswith(".cram") else out_bam[:-4]
        ext = ".cram" if out_bam.endswith(".cram") else ".bam"
        extract = extract.replace(",", " ")
        input_bam = self.bam
        # crop reads inside extract region, mainly for long reads
        if crop_first:
            input_bam = out_prefix + ".cropped" + ext
            self.clip_to_region(input_bam, extract)
        # create ref of target region
        with tempfile.NamedTemporaryFile(
            suffix=".fa", dir=self.tmpdir, delete=False
        ) as temp_file:
            realign_ref = temp_file.name
        _ = self.make_fasta(realign_ref, target)
        chr = target.split(",")[0].split(":")[0]
        offsets = [int(t.split(":")[1].split("-")[0]) - 1 for t in target.split(",")]
        # lift over
        liftover_tmp_bam = out_prefix + ".tmp.bam"
        if (
            not self.use_existing
            or not os.path.exists(liftover_tmp_bam)
            or not os.path.getsize(liftover_tmp_bam)
        ):
            x = "map-pb" if self.long_read else "sr"
            rg_data = self.get_header_field("RG")[0]
            rg = "@RG\\tID:%s\\t%s" % (
                rg_data["ID"],
                "\\t".join([f"{k}:{v}" for k, v in rg_data.items() if k != "ID"]),
            )
            # sq_field = [
            #     "@SQ\tSN:%s\tLN:%s" % (v["SN"], v["LN"])
            #     for v in self.get_header_field("SQ")
            #     if "_" not in v["SN"] and "*" not in v["SN"]
            # ]
            # sq = "\n".join(sq_field)
            cmd = f"{self.samtools} view -F 0x100 -F 0x200 -F 0x800 -h --reference {self.ref} {input_bam} {extract} |"
            cmd += f"{self.samtools} fastq -t |"
            cmd += f'{self.sentieon} minimap2 -y -a -x {x} -R "{rg}" {realign_ref} /dev/stdin |'
            cmd += f"{self.sentieon} util sort -o {liftover_tmp_bam} --sam2bam -i -"
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            result.check_returncode()
        liftover_tmp_bamh = pysam.AlignmentFile(liftover_tmp_bam, "rb")
        out_bamh = self.init_outbamh(out_bam)
        read_names = set()
        for read in liftover_tmp_bamh.fetch(until_eof=True):
            if read.is_unmapped:
                if read.query_sequence:
                    read_seq = (
                        self.rc(read.query_sequence)
                        if read.is_reverse
                        else read.query_sequence
                    )
                    read_names.add(read.query_name.split("/")[0] + read_seq[:20])
                continue
            offset = offsets[read.reference_id]
            a = pysam.AlignedSegment(out_bamh.header)
            for attr in dir(read):
                if attr.startswith("__"):
                    continue
                value = getattr(read, attr)
                try:
                    setattr(a, attr, value)
                except Exception:
                    pass
            a.reference_start += offset
            a.reference_name = chr
            a.mapping_quality = read.mapping_quality
            a.query_qualities = read.query_qualities
            if read.has_tag("SA"):
                sa_tag = read.get_tag("SA")
                if sa_tag:
                    new_sa_tag = []
                    for sa in sa_tag.split(";"):
                        if sa:
                            flds = sa.split(",")
                            flds[0] = chr
                            flds[1] = str(int(flds[1]) + offset)
                            new_sa_tag.append(",".join(flds))
                        else:
                            new_sa_tag.append(sa)
                    a.set_tag("SA", ";".join(new_sa_tag))
            out_bamh.write(a)
        liftover_tmp_bamh.close()
        out_bamh.close()
        pysam.index(out_bam)
        # get the ref positions of all the unmapped reads
        ref_positions = {}
        self.init_bamh()
        for region in extract.split():
            for read in self.bamh.fetch(region=region):
                if not read.query_sequence:
                    continue
                cur_pos = ref_positions.setdefault(region.split(":")[0], {})
                read_seq = (
                    self.rc(read.query_sequence)
                    if read.is_reverse
                    else read.query_sequence
                )
                key = read.query_name + read_seq[:20]
                if key in read_names:
                    ref_pos = set(
                        [
                            pos
                            for pos in read.get_reference_positions(full_length=True)
                            if pos
                        ]
                    )
                    for pos in ref_pos:
                        if pos in mapping:
                            cur_pos[pos] = cur_pos.get(pos, 0) + 1
        self.bamh.close()
        self.bamh = None
        self.init_bamh()
        unmapped_regions = self.to_region(ref_positions)
        # check the percetage of unmapped reads in these positions
        for region in unmapped_regions:
            chr = region.split(":")[0]
            start, end = [int(s) for s in region.split(":")[1].split("-")]
            cur_pos = ref_positions[chr]
            for pileupcol in self._pileup_safe(self.bamh, region=region, truncate=True):
                pos = pileupcol.reference_pos
                total_cnt = pileupcol.get_num_aligned()
                if pos in cur_pos and cur_pos[pos] < min_ratio * total_cnt:
                    del ref_positions[chr][pos]
        failed_region = []
        for region in self.to_region(ref_positions, 10):
            chr, positions = region.split(":")
            start, end = [int(s) for s in positions.split("-")]
            failed_region.append(f"{chr}:{mapping[start] + 1}-{mapping[end - 1] + 2}")
        failed_region = ",".join(failed_region)
        return {"lifted_bam": out_bam, "failed_region": failed_region}

    def clip2gene(self, gene: Any) -> None:
        interval = ",".join(gene.realign_regions)
        fname = os.path.basename(self.bam)
        base, ext = os.path.splitext(fname)
        self.clipped_bam = f"{self.tmpdir}/{base}.{gene.gene_names[0]}{ext}"
        cmd = f"{self.sentieon} driver -r {self.ref} --interval {interval} --interval_padding 2000 -i {self.bam} --algo ReadWriter {self.clipped_bam}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def call_variant_all_regions(
        self,
        gene: Any,
        data: List[Dict[str, Any]],
        bam_key: str,
        params: Dict[str, Any],
    ) -> None:
        if "model" not in params:
            params["model"] = self.lr_model if self.long_read else self.sr_model
        seq_key = "long_read" if self.long_read else "short_read"
        if not self.clipped_bam:
            self.clip2gene(gene)
        # to-do: missing variants in long-del region from liftover
        for i, d in enumerate(data):
            for region_id, dd in d[bam_key].items():
                if dd["cn"] != 0 and dd["region"]:
                    params["region"] = dd["region"]
                    params["ploidy"] = dd["cn"]
                    params["dbsnp"] = gene.diff_vcf[i]
                    reg_name = "nodel" if region_id == "nodel" else "del" + region_id
                    out_vcf = f"{params['prefix']}.{gene.liftover_region_names[i]}.{bam_key}.{reg_name}.raw.vcf.gz"
                    phased_vcf = f"{params['prefix']}.{gene.liftover_region_names[i]}.{bam_key}.{reg_name}.phased.vcf.gz"
                    self.call_dnascope(out_vcf, params)
                    self.phase_vcf(out_vcf, phased_vcf, dd["cn"])
                    dd[seq_key]["raw_vcf"] = out_vcf
                    dd[seq_key]["phased_vcf"] = phased_vcf
        self.clipped_bam = None

    def calc_logprob(
        self,
        region: str,
        cn: int,
        ads: Dict[int, List[int]] = {},
        ratio: float = -1.0,
        err_rate: float = 0.01,
        w_ad: float = 0.5,  # weight for the allele depth (AD) evidence
        log_prior: float = 0.0,  # log prior for CN
    ) -> float:
        probe_depths = [2 * d.mean for d in self.get_depth(region)]
        num_probes = len(probe_depths)
        if num_probes == 0:
            return log_prior

        # Step 1: Calculate total log probabilities
        log_prob_depth_sum = sum(
            CopyNumberModel.gamma_logpdf(
                max(cn, err_rate), max(depth, err_rate), self.overall_depth_std
            )
            for depth in probe_depths
        )

        mean, std = np.mean(probe_depths), np.std(probe_depths)

        if not ads or ratio < 0:
            # No AD evidence: use depth only, scaled by size
            scaled_total_log_prob = log_prob_depth_sum + log_prior
            self.logger.debug(
                f"Region {region} (CN={cn}, mean={mean:.3f}, std={std:.3f}): logprob={scaled_total_log_prob:.6f} (depth_sum={log_prob_depth_sum:.6f}, prior={log_prior:.6f})"
            )
            return scaled_total_log_prob

        p = np.clip(ratio, err_rate, 1 - err_rate)
        # Calculate mean ADS ratio
        a_ads = np.array(list(ads.values()))
        total_first = np.sum(a_ads[:, 0])
        total_both = np.sum(a_ads)
        mean_ads_ratio = total_first / total_both if total_both > 0 else -1

        # Calculate AD log probability sum
        log_prob_ad_sum = binom.logpmf(
            k=a_ads[:, 0], n=np.sum(a_ads, axis=1), p=p
        ).sum()

        # Step 2: Hybrid approach - normalize to comparable scales
        avg_log_prob_depth = log_prob_depth_sum / num_probes
        avg_log_prob_ad = log_prob_ad_sum / len(ads)

        # Step 3: Fair weighting and scaling by segment size
        quality_score = (1 - w_ad) * avg_log_prob_depth + w_ad * avg_log_prob_ad
        scaled_total_log_prob = quality_score * num_probes + log_prior

        self.logger.debug(
            f"Region {region} (CN={cn}, mean={mean:.3f}, std={std:.3f}, ratio={ratio:.3g}, mean_ads_ratio={mean_ads_ratio:.3g}): logprob={scaled_total_log_prob:.6f} (avg_depth={avg_log_prob_depth:.6f}, avg_ad={avg_log_prob_ad:.6f}, quality={quality_score:.6f}, prior={log_prior:.6f})"
        )
        return scaled_total_log_prob


class CopyNumberModel:
    cn_priors = {
        0: 0.2,
        1: 0.3,
        2: 0.4,
        3: 0.3,
        4: 0.2,
        5: 0.01,
        6: 0.001,
        "default": 1e-4,
    }
    error_rate = 0.01
    PileUpRecord = namedtuple("PileUpRecord", "pos, MQ, AD")

    def __init__(
        self, read_data: Dict[str, Any], gene: Any, param: Dict[str, Any]
    ) -> None:
        self.read_data = read_data
        self.gene = gene
        self.params = param
        self.logger = get_logger(self.__class__.__name__)

    @staticmethod
    def get_ads(
        fname: str,
        reg: Union[str, IntervalList, None],
        filter_fn: Optional[Callable] = None,
        with_coor: bool = False,
    ) -> Union[List[List[int]], Dict[int, List[int]]]:
        if reg is None:
            return []
        if isinstance(reg, str):
            reg = IntervalList(region=reg)
        vcf = vcflib.VCF(fname)
        ads = []
        pos = []
        for v in vcf:
            if (
                (filter_fn is None or filter_fn(v))
                and reg.get(v.chrom, v.pos, v.end)
                and "AD" in v.samples[0]
            ):
                ads.append(v.samples[0]["AD"])
                pos.append(v.pos)
        vcf.close()
        if with_coor:
            return dict(zip(pos, ads))
        return ads

    def get_ads_bygiven(
        self,
        bam: "Bam",
        given: str,
        reg: Union[str, IntervalList, None],
        filter_fn: Optional[Callable] = None,
        alt: bool = False,
    ) -> Dict[int, List[int]]:
        if reg is None:
            return {}
        if isinstance(reg, str):
            reg = IntervalList(region=reg)
        vcf = vcflib.VCF(given)
        cnts = {}
        bamh = bam.init_bamh()
        for v in vcf:
            if (
                reg.get(v.chrom, v.pos, v.end)
                and len(v.ref) == 1
                and len(v.alt) == 1
                and len(v.alt[0]) == 1
            ):
                region = f"{v.chrom}:{v.pos + 1}-{v.pos + 1}"
                for pileupcol in bam._pileup_safe(bamh, region=region, truncate=True):
                    pos = pileupcol.reference_pos
                    if pos != v.pos:
                        continue
                    seq = [
                        s.upper()
                        for s in pileupcol.get_query_sequences(add_indels=True)
                    ]
                    c = Counter(seq)
                    ref_cnt = c[v.ref]
                    alt_cnt = len(seq) - ref_cnt
                    if filter_fn:
                        mqs = pileupcol.get_mapping_qualities()
                        mq = sum(mqs) / len(mqs)
                        if not filter_fn(
                            CopyNumberModel.PileUpRecord(v.pos, mq, (ref_cnt, alt_cnt))
                        ):
                            continue
                    cnts[v.pos] = [ref_cnt, alt_cnt]
                if alt and v.pos in cnts:
                    cnts[v.pos] += [0, 0]
                    pos2 = self.gene.mapping.get(v.pos)
                    region = f"{v.chrom}:{pos2 + 1}-{pos2 + 1}"
                    for pileupcol in bam._pileup_safe(
                        bamh, region=region, truncate=True
                    ):
                        pos = pileupcol.reference_pos
                        if pos != pos2:
                            continue
                        seq = [
                            s.upper()
                            for s in pileupcol.get_query_sequences(add_indels=True)
                        ]
                        c = Counter(seq)
                        ref_cnt2 = c[v.alt[0]]
                        alt_cnt2 = len(seq) - ref_cnt2
                        cnts[v.pos][2:] = [ref_cnt2, alt_cnt2]
        vcf.close()
        return cnts

    @staticmethod
    def gamma_logpdf(n: float, mean_values: float, sigma: float) -> float:
        n = max(0.01, n)
        sigma = max(0.01, sigma)
        shape = (n / sigma) ** 2  # Shape parameter
        scale = sigma**2 / n
        return gamma.logpdf(mean_values, a=shape, scale=scale)

    def detect_longdels(
        self,
        non_del_region: str,
        known_dels: List[str],
        bam: "Bam",
        longdel_names: List[str] = None,
    ) -> List[int]:
        if not known_dels:
            return []
        all_depths = [d.mean for d in bam.get_depth(non_del_region)]
        known_del_depths = []
        for del_region in known_dels:
            known_del_depths.append([d.mean for d in bam.get_depth(del_region)])
        mean, std = np.mean(all_depths) * 2, np.std(all_depths) * 2
        del_stats = [(np.mean(d) * 2, np.std(d) * 2) for d in known_del_depths]
        init_depth = int(round(mean))
        log_pdf = {
            n: self.gamma_logpdf(n, mean, std)
            for n in range(max(0, init_depth - 2), init_depth + 3)
        }
        best_depth = max(log_pdf, key=log_pdf.get)
        diffs = []
        for i, stats in enumerate(del_stats):
            m, s = stats

            # Get longdel prior if available
            longdel_log_prior = 0.0  # Default: no prior
            if longdel_names and i < len(longdel_names):
                longdel_name = longdel_names[i]
                if (
                    hasattr(self.gene, "_longdel_log_priors")
                    and longdel_name in self.gene._longdel_log_priors
                ):
                    longdel_log_prior = self.gene._longdel_log_priors[longdel_name]

            if best_depth > 1:
                # Calculate likelihood + prior for each possible deletion state
                del_posteriors = []
                for diff in range(best_depth):
                    likelihood = CopyNumberModel.gamma_logpdf(best_depth - diff, m, s)
                    # Prior: deletion present if diff > 0, absent if diff == 0
                    if diff > 0:
                        posterior = likelihood + longdel_log_prior
                    else:
                        # Use log(1 - p) for no deletion, where p is deletion probability
                        no_del_log_prior = (
                            np.log(1 - np.exp(longdel_log_prior))
                            if longdel_log_prior > -1e6
                            else 0.0
                        )
                        posterior = likelihood + no_del_log_prior
                    del_posteriors.append(posterior)

                diff = del_posteriors.index(max(del_posteriors))
                diffs.append(diff)
            else:
                # Simple case: use prior to bias decision
                if mean - m > std + s and m < 0.5:
                    # Evidence suggests deletion
                    if longdel_log_prior > np.log(0.5):  # Prior favors deletion
                        diffs.append(1)
                    else:
                        diffs.append(1 if longdel_log_prior > np.log(0.1) else 0)
                else:
                    # Weak evidence: use prior
                    diffs.append(1 if longdel_log_prior > np.log(0.5) else 0)
        return diffs

    def get_segdup_ads(self, seq_key: str) -> Dict[int, List[int]]:
        liftover_bam = self.read_data[seq_key]["liftover"]
        bam = self.read_data[seq_key]["bam"]
        excl_pos = []
        var_ads = []
        liftover_region = (
            [",".join(self.gene.liftover_target_regions)]
            if len(self.gene.diff_vcf) == 1
            else self.gene.liftover_target_regions
        )
        for i in range(len(self.gene.diff_vcf)):
            pos_ads1 = self.get_ads_bygiven(
                bam,
                self.gene.diff_vcf[i],
                self.gene.high_mq_regions[seq_key],
                filter_fn=lambda v: v.MQ > 50,
                alt=True,
            )
            excl_pos += [
                pos
                for pos, ads in pos_ads1.items()
                if (
                    ads[0] + ads[1] > 0
                    and (ads[1] / (ads[0] + ads[1]) > 0.1 or ads[1] > 5)
                )
                or (
                    ads[2] + ads[3] > 0
                    and (ads[3] / (ads[2] + ads[3]) > 0.1 or ads[3] > 5)
                )
            ]
            pos_ads = self.get_ads_bygiven(
                liftover_bam, self.gene.diff_vcf[i], liftover_region[i]
            )
            var_ads += [
                (p, ads)
                for p, ads in pos_ads.items()
                if p not in excl_pos
                and (
                    p not in pos_ads1
                    or pos_ads1[p][0] - ads[0] <= 5
                    and pos_ads1[p][2] - ads[1] <= 5
                )
            ]
        return dict(var_ads)


class Phased_vcf:
    max_mismatch = 0.25

    def __init__(self, gene: Any, read_data: Dict[str, Any]) -> None:
        self.gene = gene
        self.read_data = read_data
        self.logger = get_logger(self.__class__.__name__)

    @staticmethod
    def _annotate(v: Any, dp: Optional[int] = None) -> None:
        gt = [int(g) for g in re.split(r"\||/", v.samples[0]["GT"])]
        ac = [gt.count(i + 1) for i in range(len(v.alt))]
        af = [a / len(gt) for a in ac]
        if not dp:
            dp = int(round(v.info["DP"] / v.info["AN"] * len(gt)))
        v.info["AC"] = ac
        v.info["MLEAC"] = ac
        v.info["AN"] = len(gt)
        v.info["AF"] = af
        v.info["MLEAF"] = af
        v.info["DP"] = dp
        v.samples[0]["DP"] = dp
        del v.samples[0]["PL"]
        del v.samples[0]["AD"]

    @staticmethod
    def split_homvar(v: Any, cns: Tuple[int, int]) -> Optional[Any]:
        v1 = None
        gts = set(re.split(r"\||/", v.samples[0]["GT"]))
        if len(gts) == 1 and gts != set("0"):
            v1 = copy.deepcopy(v)
            v1.line = None
            v1.samples[0]["GT"] = "|".join(["1"] * cns[0])
            Phased_vcf._annotate(v1)
        return v1

    def proc_phased(
        self,
        phased: List[Any],
        cns: Tuple[int, int],
        matched: Dict[int, Any],
        gene_id: int,
    ) -> Tuple[List[Any], List[Any]]:
        def _proc_homvar():
            for v in phased:
                v1 = self.split_homvar(v, cns)
                if v1:
                    v1s.append(v1)

        v1s = []
        all_alleles = []
        for v in phased:
            assert sum(cns) == len(re.split(r"\||/", v.samples[0]["GT"]))
            if v.id != ".":
                all_alleles.append(
                    "_".join(
                        [
                            str(i)
                            for i, g in enumerate(v.samples[0]["GT"].split("|"))
                            if g == "1"
                        ]
                    )
                )
            elif v.pos in matched:
                m = matched[v.pos]
                if m[0]:
                    all_alleles.append(
                        "_".join(
                            [
                                str(i)
                                for i, g in enumerate(v.samples[0]["GT"].split("|"))
                                if g == "0"
                            ]
                        )
                    )
                elif m[1]:
                    all_alleles.append(
                        "_".join(
                            [
                                str(i)
                                for i, g in enumerate(v.samples[0]["GT"].split("|"))
                                if g == "1"
                            ]
                        )
                    )

        if not all_alleles:
            _proc_homvar()
            return phased, v1s
        # get most common choice
        # if there is any alleles inconsistent with CN
        allele_cnts = Counter(all_alleles)
        wrong_alleles = [a for a in set(all_alleles) if len(a.split("_")) != cns[1]]
        if not wrong_alleles:
            gene2_alleles = allele_cnts.most_common()[0]
        else:
            correct_alleles = [a for a in all_alleles if a not in wrong_alleles]
            if not correct_alleles:
                _proc_homvar()
                return phased, v1s
            correct_allele_cnts = Counter(correct_alleles)
            allele, cnt = correct_allele_cnts.most_common()[0]
            sallele = set(allele.split("_"))
            # check alleles with wrong count to see if it is consistent with the most commonly allele with correct CN
            for a in wrong_alleles:
                sa = a.split("_")
                if sallele.intersection(sa) == set(sa):  # is a subset
                    cnt += allele_cnts[a]
            gene2_alleles = [allele, cnt]
        # no consensus alleles
        if gene2_alleles[1] < (1 - self.max_mismatch) * len(all_alleles):
            _proc_homvar()
            return phased, v1s
        gene2_alleles = [int(i) for i in gene2_alleles[0].split("_")]
        gene1_alleles = [i for i in range(sum(cns)) if i not in gene2_alleles]
        for v in phased:
            label = []
            v.line = None
            alt_alleles = set(
                [
                    i
                    for i, g in enumerate(list(map(int, v.samples[0]["GT"].split("|"))))
                    if g == 1
                ]
            )
            if alt_alleles.difference(gene2_alleles):
                label.append(self.gene.liftover_region_names[gene_id])
                v1 = copy.deepcopy(v)
                v1.samples[0]["GT"] = "|".join(
                    ["1" if a in alt_alleles else "0" for a in gene1_alleles]
                )
                self._annotate(v1)
                v1s.append(v1)
            v.info["LABEL"] = ",".join(label)
        return phased, v1s

    @staticmethod
    def match(
        vars: List[Tuple[Any, Any]],
        v1s: List[Tuple[int, str, int, Any]],
        v2s: List[Tuple[int, str, int, Any]],
    ) -> Dict[int, List[Optional[Tuple[str, List[str], Any]]]]:
        matched = {}
        for pos, key, _, v in v1s:
            matched.setdefault(pos, [None, None, None])[0] = (
                key,
                re.split(r"\||/", v.samples[0]["GT"]),
                v,
            )
        for pos, key, _, v in v2s:
            matched.setdefault(pos, [None, None, None])[1] = (
                key,
                re.split(r"\||/", v.samples[0]["GT"]),
                v,
            )
        for key, v in vars:
            matched.setdefault(v.pos, [None, None, None])[2] = (
                key,
                re.split(r"\||/", v.samples[0]["GT"]),
                v,
            )
        keep = {}
        for pos, vs in matched.items():
            v1, v2, v = vs
            if v2 and v1 or not v:
                continue
            if v1 and v1[0] != v[0] or v2 and v2[0] != v[0]:
                continue
            if v1 and "0" in v1[1] or v2 and "0" in v2[1]:
                continue
            if any(x not in ("0", "1") for x in v[1]):
                continue
            if (
                v1
                and v1[1].count("1") == v[1].count("1")
                or v2
                and v2[1].count("1") == v[1].count("1")
            ):
                keep[pos] = vs
        return keep

    def concat(self, in_vcfs: List[str], out_fname: str) -> None:
        vcfs = [vcflib.VCF(vcf) for vcf in in_vcfs]
        out_vcf = vcflib.VCF(out_fname, "w")
        out_vcf.copy_header(vcfs[0])
        out_vcf.emit_header()
        vars = []
        cn = [-1, -1]
        for i, vcf in enumerate(vcfs):
            for v in vcf:
                if cn[i] < 0:
                    cn[i] = len(re.split(r"\||/", v.samples[0]["GT"]))
                heapq.heappush(vars, (v.pos, f"{v.ref}_{','.join(v.alt)}", v))
        while vars:
            _, _, v = heapq.heappop(vars)
            out_vcf.emit(v)

    @staticmethod
    def vcf_grouper(*args):
        def _key(v):
            return f"{v.ref}_{v.alt}"

        def _is_snp(v):
            return len(v.alt) == 1 and len(v.alt[0]) == 1 and len(v.ref) == 1

        def pop():
            p, ra, k, v, i = heapq.heappop(q)
            vv = next(i, None)
            if vv:
                heapq.heappush(q, (vv.pos, _key(vv), k, vv, i))
            return p, ra, k, v, i

        assert len(args) == 3
        vcfs = args[:3]
        q = []
        for k, vcf in enumerate(vcfs):
            i = iter(vcf)
            v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, _key(v), k, v, i))
        cur = [None, None, None]
        cur_pos = None
        while q:
            p, ra, k, v, i = pop()
            if cur_pos and p != cur_pos:
                if cur[2]:
                    if cur[1]:
                        cur[1].samples[0]["GT"] = cur[2].samples[0]["GT"]
                        cur[1].line = None
                        yield cur[1]
                    elif cur[0] and _key(cur[0]) == _key(cur[2]):
                        yield cur[2]
                cur = [None, None, None]
            cur_pos = p
            cur[k] = v
        if cur[2]:
            if cur[1]:
                cur[1].samples[0]["GT"] = cur[2].samples[0]["GT"]
                cur[1].line = None
                yield cur[1]
            elif cur[0] and _key(cur[0]) == _key(cur[2]):
                yield cur[2]

    def assign(
        self,
        liftover_vcf: str,
        in_vcf1: Optional[str],
        in_vcf2: Optional[str],
        params: Dict[str, Any],
    ) -> Tuple[List[Any], List[Any]]:
        vcf = vcflib.VCF(liftover_vcf)
        region1 = (
            IntervalList(region=params["region"]) if params.get("region") else None
        )
        orig_region = (
            IntervalList(region=params["orig_region"])
            if params.get("orig_region")
            else None
        )
        liftover_region = region1  # region1.subtract(orig_region) if region1 else None
        gene_id = params["gene_id"]
        orig_region2 = (
            IntervalList(region=self.gene.mapping.interval(orig_region))
            if orig_region
            else None
        )
        all_vars = [
            (f"{v.ref}_{v.alt[0]}", v)
            for v in vcf
            if liftover_region and (v.chrom, v.pos) in liftover_region
        ]
        all_v1s = []
        all_v2s = []
        if in_vcf1:
            vcf1 = vcflib.VCF(in_vcf1)
            for v1 in vcf1:
                if (
                    orig_region
                    and (v1.chrom, v1.pos) in orig_region
                    and (not v1.filter or v1.filter in (".", "PASS"))
                ):
                    heapq.heappush(all_v1s, (v1.pos, f"{v1.ref}_{v1.alt[0]}", 0, v1))
        if in_vcf2:
            vcf2 = vcflib.VCF(in_vcf2)
            to_convert = [
                v2
                for v2 in vcf2
                if orig_region2
                and (v2.chrom, v2.pos) in orig_region2
                and (not v2.filter or v2.filter in (".", "PASS"))
            ]
            direction = "backward" if gene_id == 0 else "forward"
            for v2 in self.gene.mapping.convert_vars(
                to_convert, direction, ref_diff=False
            ):
                heapq.heappush(all_v2s, (v2.pos, f"{v2.ref}_{v2.alt[0]}", 1, v2))
        cns = params["cns"]
        if len(cns) != 2:
            self.logger.error(
                "ERROR: exactly two copy numbers required for the gene and its paralog separated by comma. Example: 2,3"
            )
            sys.exit(1)
        matched = self.match(all_vars, all_v1s, all_v2s)
        phased = []
        cur_ps = None
        out_vars = []
        out_v1s = []
        for _, v in all_vars:
            ps = v.samples[0].get("PS")
            if ps:
                if cur_ps and cur_ps != ps:
                    phased, v1s = self.proc_phased(phased, cns, matched, gene_id)
                    for vv in phased:
                        heapq.heappush(out_vars, (vv.pos, f"{vv.ref}_{vv.alt[0]}", vv))
                    for vv in v1s:
                        heapq.heappush(
                            all_v1s, (vv.pos, f"{vv.ref}_{vv.alt[0]}", 1, vv)
                        )
                    phased = []
                phased.append(v)
                cur_ps = ps
            else:
                gts = set(re.split(r"\||/", v.samples[0]["GT"]))
                if len(gts) == 1 and gts != set("0"):  ## homvar
                    v1 = self.split_homvar(v, cns)
                    if v1:
                        heapq.heappush(
                            all_v1s, (v1.pos, f"{v1.ref}_{v1.alt[0]}", 1, v1)
                        )
                heapq.heappush(out_vars, (v.pos, f"{v.ref}_{v.alt[0]}", v))
        if cur_ps:
            phased, v1s = self.proc_phased(phased, cns, matched, gene_id)
            for vv in phased:
                heapq.heappush(out_vars, (vv.pos, f"{vv.ref}_{vv.alt[0]}", vv))
            for vv in v1s:
                heapq.heappush(all_v1s, (vv.pos, f"{vv.ref}_{vv.alt[0]}", 1, vv))
        last = None
        while all_v1s:
            cur = heapq.heappop(all_v1s)
            if region1 and (cur[3].chrom, cur[3].pos) not in region1:
                continue
            if last:
                if last[:2] == cur[:2]:
                    v = last[3]
                    last_gt = re.split(r"\||/", last[3].samples[0]["GT"])
                    cur_gt = re.split(r"\||/", cur[3].samples[0]["GT"])
                    if (
                        len(last_gt) == 2
                        or last_gt == cur_gt
                        or last[3].info.get("MQ", 0) > 40
                    ):
                        out_v1s.append(last[3])
                    else:
                        self.logger.debug(
                            f"{self.gene.liftover_region_names[gene_id]}: {v.chrom}:{v.pos} {v.ref}-{','.join(v.alt)}: original GT {last[3].samples[0]['GT']} replaced with {cur[3].samples[0]['GT']}."
                        )
                        out_v1s.append(cur[3])
                    last = None
                    continue
                out_v1s.append(last[3])
            last = cur
        if last:
            out_v1s.append(last[3])
        out_liftover = []
        while out_vars:
            _, _, v = heapq.heappop(out_vars)
            out_liftover.append(v)
        diff = params["cn_diff"]
        if diff < 0:
            diff = abs(diff)
            chrom = params["region"].split(":")[0]
            pos, end = [int(s) for s in params["region"].split(":")[1].split("-")]
            gt = "|".join(["0"] * cns[0] + ["1"] * diff)
            gt_liftover = "|".join(["0"] * sum(cns) + ["1"] * diff)
            ref_base = "A"  # fix me
            out_v1s.insert(
                0,
                vcflib.Variant(
                    chrom,
                    pos,
                    "",
                    ref_base,
                    ["<DEL>"],
                    None,
                    [],
                    {},
                    [{"GT": gt}],
                    end,
                    None,
                ),
            )
            out_liftover.insert(
                0,
                vcflib.Variant(
                    chrom,
                    pos,
                    "",
                    ref_base,
                    ["<DEL>"],
                    None,
                    [],
                    {},
                    [{"GT": gt_liftover}],
                    end,
                    None,
                ),
            )
            for v in out_v1s[1:]:
                gt = v.samples[0]["GT"]
                sep = "|" if "|" in gt else "/"
                v.alt.append("*")
                v.samples[0]["GT"] += sep + sep.join([str(len(v.alt))] * diff)
                if "AC" in v.info:
                    v.info["AC"].append(diff)
                if "MLEAC" in v.info:
                    v.info["MLEAC"].append(diff)
                if "PL" in v.samples[0]:
                    del v.samples[0]["PL"]
                if "AD" in v.samples[0]:
                    del v.samples[0]["AD"]
                v.line = None
            for v in out_liftover[1:]:
                gt = v.samples[0]["GT"]
                sep = "|" if "|" in gt else "/"
                v.alt.append("*")
                v.samples[0]["GT"] += sep + sep.join([str(len(v.alt))] * diff)
                if "AC" in v.info:
                    v.info["AC"].append(diff)
                if "MLEAC" in v.info:
                    v.info["MLEAC"].append(diff)
                if "PL" in v.samples[0]:
                    del v.samples[0]["PL"]
                if "AD" in v.samples[0]:
                    del v.samples[0]["AD"]
                v.line = None
        return out_liftover, out_v1s
