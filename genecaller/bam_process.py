import vcflib
import pysam
import os
import subprocess
import sys
from .util import IntervalList, get_data_file
from .logging import get_logger
import heapq
from collections import namedtuple, Counter
import statistics
from scipy.stats import binom, gamma
import tempfile
import numpy as np
import pandas as pd
import shutil
import re
import copy

Segment = namedtuple("Segment", "chrom, start, end, cnt, mq, mean, std, mean0, std0")
DepthSegment = namedtuple("DepthSegment", "chrom, start, end, mq, mean, mean0")


class Bam:
    def __init__(self, bam, ref, param):
        self.bam = bam
        self.bamh = None
        self.header = None
        self.long_read = False
        self.ref = ref
        self.depths = {}
        self.segments = {}
        self.cnv_model = ""
        self.lr_model = ""
        self.sr_model = ""
        self.prefix = None
        self.sample_name = "SAMPLE"
        self.use_existing = True
        self.clipped_bam = None
        self.outdir = None
        self.tmpdir = param["tmpdir"]
        if "tools" not in param:
            raise Exception("External tools are not set up")
        for k, v in param["tools"].items():
            setattr(self, k, v)
        for k, v in param.items():
            if hasattr(self, k):
                setattr(self, k, v)
        if self.cnv_model and not os.path.exists(self.cnv_model) and not os.path.exists(os.path.dirname(self.cnv_model)):
            self.cnv_model = get_data_file(self.cnv_model)
        if self.lr_model and not os.path.exists(self.lr_model) and not os.path.exists(os.path.dirname(self.lr_model)):
            self.lr_model = get_data_file(self.lr_model)
        if self.sr_model and not os.path.exists(self.sr_model) and not os.path.exists(os.path.dirname(self.sr_model)):
            self.sr_model = get_data_file(self.sr_model)
        self.logger = get_logger(self.__class__.__name__, "INFO")

    def reset(self):
        self.depths = {}
        self.segments = {}
        self.clipped_bam = None

    @staticmethod
    def rc(seq):
        rev = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
        return "".join([rev[s] for s in seq[-1::-1]])

    def call_dnascope(self, output, param):
        if self.use_existing and os.path.exists(output) and os.path.getsize(output):
            return
        bam = self.clipped_bam if self.clipped_bam else self.bam
        driver_opt = param.get("driver_opt", "") + f" -i {bam} -r {self.ref}"
        algo_opt = param.get("algo_opt", "")
        apply_model = ""
        if "model" in param:
            algo_opt += f" --model {param['model']}"
            apply_model = param["model"]
        if "region" in param:
            driver_opt += f" --interval {param['region']}"
        if "ploidy" in param:
            algo_opt += f" --ploidy {param['ploidy']}"
            if param["ploidy"] != 2:
                apply_model = ""
        if "given" in param:
            algo_opt += f" --given {param['given']} -d {param['given']}"
            apply_model = ""
        elif "dbsnp" in param:
            algo_opt += f" -d {param['dbsnp']}"
        if apply_model:
            tmp_output = output.replace("vcf", "tmp.vcf")
        else:
            tmp_output = output
        cmd = f"{self.sentieon} driver {driver_opt} --algo DNAscope {algo_opt} {tmp_output}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if result.returncode:
            raise Exception(result.stderr)
        if apply_model:
            cmd = f"{self.sentieon} driver -r {self.ref} --algo DNAModelApply -v {tmp_output} --model {apply_model} {output}"
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.returncode:
                raise Exception(result.stderr)

    def call_genotyper(self, output, param):
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

    def call_cnvscope(self, output, param):
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

    def phase_vcf(self, in_vcf, out_vcf, ploidy=2):
        if self.use_existing and os.path.exists(out_vcf) and os.path.getsize(out_vcf):
            return
        bam = self.clipped_bam if self.clipped_bam else self.bam
        if ploidy == 2:
            cmd = f"{self.whatshap} phase -o {out_vcf} {in_vcf} {bam} --reference {self.ref}"
        else:
            cmd = f"{self.whatshap} polyphase -o {out_vcf} {in_vcf} {bam} -p {ploidy} --reference {self.ref} -B 5"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()
        cmd = f"{self.sentieon} util vcfindex {out_vcf}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def get_depth(self, regions):
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
        depths = [
            d
            for chr in regions.regions
            for d in self.depths[chr]
            if (d.chrom, d.start) in regions or (d.chrom, d.end) in regions
        ]
        return depths

    def get_segments(self, regions):
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
            output = os.path.join(
                f"{self.tmpdir}/{self.sample_name}.{self.prefix}.cnv.{'.'.join(chrs)}.vcf.gz"
            )
            self.call_depth(output, param)
        segments = set()
        for chr, reg in regions.regions.items():
            for i in range(len(reg) // 2):
                s, e = regions.regions[chr][i * 2], regions.regions[chr][i * 2 + 1]
                segments = segments.union(
                    [seg for seg in self.segments[chr] if seg.start < e and seg.end > s]
                )
        return segments

    def call_depth(self, output, param):
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
                    dp[0],
                    dp[1],
                    dp0[0],
                    dp0[0],
                )
            )
        depth_df = pd.read_csv(tmp_output_norm, sep="\t")[depth_fields]
        self.depths = dict(
            [
                (
                    chr,
                    [
                        DepthSegment(*row)
                        for row in depth_df[depth_df["contig"] == chr].itertuples(
                            index=False
                        )
                    ],
                )
                for chr in chrs
            ]
        )
        if orig_region:
            orig_region_itv = IntervalList(region=orig_region)
            trim_depth_df = depth_df[
                depth_df.apply(
                    lambda r: (r["contig"], r["start"]) in orig_region_itv
                    or (r["contig"], r["end"]) in orig_region_itv,
                    axis=1,
                )
            ]
            trim_depth_df.to_csv(output_norm, sep="\t", index=None)
            self.trim_to_region_vcf(orig_region, tmp_output, output)

    def trim_to_region_vcf(self, region, in_vcf, out_vcf):
        if region.endswith(".bed"):
            region_opt = f"-T {region}"
        else:
            region_opt = f"-r {region}"
        cmd = f"{self.bcftools} view {region_opt} {in_vcf} | {self.sentieon} util vcfconvert - {out_vcf}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def init_bamh(self):
        if self.bamh:
            self.bamh.close()
            self.bamh = None
        if self.bam.endswith(".cram"):
            if not self.ref:
                print("No reference file provided for cram input", file=sys.stderr)
                sys.exit(1)
            self.bamh = pysam.AlignmentFile(self.bam, "rc", reference_filename=self.ref)
        else:
            self.bamh = pysam.AlignmentFile(self.bam, "rb")
        return self.bamh

    def init_outbamh(self, out_bam):
        self.init_bamh()
        if out_bam.endswith(".cram"):
            out_bamh = pysam.AlignmentFile(
                out_bam, "wc", reference_filename=self.ref, template=self.bamh
            )
        else:
            out_bamh = pysam.AlignmentFile(out_bam, "wb", template=self.bamh)
        return out_bamh

    def get_header_field(self, field):
        if not self.header:
            self.init_bamh()
            self.header = self.bamh.header
        return self.header.get(field)

    def clip_to_region(self, out_bam, regions, min_length=300):
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

    def make_fasta(self, ref_file, region):
        make_ref_cmd = f"{self.samtools} faidx {self.ref} {region} > {ref_file}"
        result = subprocess.run(
            make_ref_cmd, capture_output=True, text=True, shell=True
        )
        result.check_returncode()
        pysam.faidx(ref_file)
        with open(ref_file) as f:
            for line in f:
                contig_name = line[1:].strip()
                break
        return contig_name

    # convert 0-based to 1-based regions
    @staticmethod
    def to_region(positions, min_gap=1):
        regions = []
        for chr, pos in positions.items():
            last_pos = -1
            last_region = None
            sorted_pos = sorted(pos.keys())
            for p in sorted_pos:
                if p > last_pos + min_gap:
                    if last_region:
                        last_region += str(last_pos + 2)
                        regions.append(last_region)
                    last_region = f"{chr}:{p + 1}-"
                last_pos = p
            if last_pos and last_region:
                last_region += str(last_pos + 2)
                regions.append(last_region)
        return regions

    def liftover_combine(self, gene, out_bam):
        out_prefix = out_bam[:-5] if out_bam.endswith(".cram") else out_bam[:-4]
        liftover_bams = []
        results = []
        for i, gene_region in enumerate(gene.gene_regions):
            liftover_bams.append(f"{out_prefix}.{gene.gene_names[i]}.bam")
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
        self, extract, target, out_bam, mapping, crop_first=False, min_ratio=0.3
    ):
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
        offset = int(target.split(",")[0].split(":")[1].split("-")[0]) - 1
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
        min_pos = min(mapping.keys())
        max_pos = max(mapping.keys())
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
                        if pos >= min_pos and pos < max_pos:
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
            for pileupcol in self.bamh.pileup(region=region, truncate=True):
                pos = pileupcol.reference_pos
                total_cnt = pileupcol.get_num_aligned()
                if pos in cur_pos and cur_pos[pos] < min_ratio * total_cnt:
                    del ref_positions[chr][pos]
        failed_region = []
        for region in self.to_region(ref_positions, 10):
            chr, positions = region.split(":")
            start, end = [int(s) for s in positions.split("-")]
            failed_region.append(f"{chr}:{mapping[start]}-{mapping[end - 1] + 1}")
        failed_region = ",".join(failed_region)
        return {"lifted_bam": out_bam, "failed_region": failed_region}

    def clip2gene(self, gene):
        interval = ",".join(gene.realign_regions)
        fname = os.path.basename(self.bam)
        base, ext = os.path.splitext(fname)
        self.clipped_bam = f"{self.tmpdir}/{base}.{gene.gene_names[0]}{ext}"
        cmd = f"{self.sentieon} driver -r {self.ref} --interval {interval} --interval_padding 2000 -i {self.bam} --algo ReadWriter {self.clipped_bam}"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()

    def call_variant_all_regions(self, gene, data, bam_key, params):
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
                    out_vcf = f"{params['prefix']}.{gene.gene_names[i]}.{bam_key}.{reg_name}.raw.vcf.gz"
                    phased_vcf = f"{params['prefix']}.{gene.gene_names[i]}.{bam_key}.{reg_name}.phased.vcf.gz"
                    self.call_dnascope(out_vcf, params)
                    self.phase_vcf(out_vcf, phased_vcf, dd["cn"])
                    dd[seq_key]["raw_vcf"] = out_vcf
                    dd[seq_key]["phased_vcf"] = phased_vcf
        self.clipped_bam = None


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

    def __init__(self, read_data, gene, param):
        self.read_data = read_data
        self.gene = gene
        self.params = param
        self.logger = get_logger(self.__class__.__name__, "INFO")

    @staticmethod
    def get_ads(fname, reg, filter_fn=None, with_coor=False):
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

    def get_ads_bygiven(self, bam, given, reg, filter_fn=None, alt=False):
        if reg is None:
            return []
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
                for pileupcol in bamh.pileup(region=region, truncate=True):
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
                    pos2 = self.gene.mapping.g1tog2[v.pos]
                    region = f"{v.chrom}:{pos2 + 1}-{pos2 + 1}"
                    for pileupcol in bamh.pileup(region=region, truncate=True):
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
    def gamma_logpdf(n, mean_values, sigma):
        n = max(0.01, n)
        k = n**2 / sigma**2  # Shape parameter
        theta = sigma**2 / n  # Scale parameter
        return gamma.logpdf(mean_values, a=k, scale=theta)

    @staticmethod
    def calc_loglikelihood(ads, cn1, cn2, params={}):
        priors = params.get("priors", None)
        segs1 = params.get("segment1", None)
        segs2 = params.get("segment2", None)
        segs_sum = params.get("segment_total", None)
        err_rate = params.get("error_rate", 0.01)
        if "default" not in priors:
            priors["default"] = CopyNumberModel.cn_priors["priors"]
        if priors is None:
            log_priors = [0, 0, 0]
        else:
            total_priors = 0
            for n1, p1 in priors.items():
                for n2, p2 in priors.items():
                    try:
                        if n1 + n2 == cn1 + cn2:
                            total_priors += p1 + p2
                    except TypeError:
                        pass
            log_priors = [
                np.log(priors.get(cn1, priors["default"])),
                np.log(priors.get(cn2, priors["default"])),
                np.log(total_priors),
            ]
        ratio = min(cn1 / (cn1 + cn2) + err_rate, 1 - err_rate)
        log_cn = 0
        for segs in segs1:
            if segs:
                log_cn += CopyNumberModel.gamma_logpdf(
                    max(cn1, err_rate), segs[0], segs[1]
                )
        for segs in segs2:
            if segs:
                log_cn += CopyNumberModel.gamma_logpdf(
                    max(cn2, err_rate), segs[0], segs[1]
                )
        for segs in segs_sum:
            if segs:
                log_cn += CopyNumberModel.gamma_logpdf(
                    max(cn1 + cn2, err_rate), segs[0], segs[1]
                ) + sum(log_priors)
        if ads:
            a_ads = np.array(ads)
            return (
                log_cn * len(ads)
                + binom.logpmf(a_ads[:, 0], np.sum(a_ads, axis=1), ratio).sum()
            )
        return log_cn

    def detect_longdels(self, non_del_region, known_dels):
        if not known_dels:
            return []
        bam = self.read_data["short_read"]["liftover"]
        all_depths = [d.mean for d in bam.get_depth(non_del_region)]
        known_del_depths = []
        for del_region in known_dels:
            known_del_depths.append([d.mean for d in bam.get_depth(del_region)])
        mean, std = statistics.mean(all_depths) * 2, statistics.stdev(all_depths) * 2
        del_stats = [
            (statistics.mean(d) * 2, statistics.stdev(d) * 2) for d in known_del_depths
        ]
        init_depth = int(round(mean))
        log_pdf = dict(
            [
                (n, self.gamma_logpdf(n, mean, std))
                for n in range(max(0, init_depth - 2), init_depth + 3)
            ]
        )
        best_depth = max(log_pdf, key=log_pdf.get)
        diffs = []
        for i, stats in enumerate(del_stats):
            m, s = stats
            if best_depth > 1:
                del_pdf = [
                    CopyNumberModel.gamma_logpdf(best_depth - i, m, s)
                    for i in range(best_depth)
                ]
                diff = del_pdf.index(max(del_pdf))
                diffs.append(diff)
            else:
                if mean - m > std + s and m < 0.5:
                    diffs.append(1)
                else:
                    diffs.append(0)
        return diffs

    def call_cn(self, data, seq_keys, params):
        if "cn_priors" not in params:
            params["priors"] = self.cn_priors
        if "error_rate" not in params:
            params["error_rate"] = self.error_rate
        total_cn_ub = 8
        cn1_ub = 8
        cn2_ub = 8
        for region_id, d in data[0]["orig"].items():
            self.logger.debug(f"Call copy number for {region_id}")
            segs1 = []
            segs2 = []
            segs_total = []
            var_ads = []
            for seq_key in seq_keys:
                liftover = data[0]["liftover"][region_id][seq_key]
                orig = d[seq_key]
                orig2 = data[1]["orig"][region_id][seq_key]
                liftover_region = liftover["region"]
                high_mq_region = orig["region"]
                high_mq_region2 = orig2["region"]
                dps = [
                    [2 * d.mean for d in bam.get_depth(region)]
                    for bam, region in (
                        (self.read_data[seq_key]["bam"], high_mq_region),
                        (self.read_data[seq_key]["bam"], high_mq_region2),
                        (self.read_data[seq_key]["liftover"], liftover_region),
                    )
                ]
                read_stats = [
                    (statistics.mean(dp), statistics.stdev(dp)) if dp else None
                    for dp in dps
                ]
                segs1.append(read_stats[0])
                segs2.append(read_stats[1])
                segs_total.append(read_stats[2])
                liftover_bam = self.read_data[seq_key]["liftover"]
                bam = self.read_data[seq_key]["bam"]
                excl_pos = []
                pos_ads1 = self.get_ads_bygiven(
                    bam,
                    self.gene.diff_vcf[0],
                    data[0]["orig"][region_id]["region"],
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
                    liftover_bam, self.gene.diff_vcf[0], liftover_region
                )
                var_ads += [
                    ads
                    for p, ads in pos_ads.items()
                    if p not in excl_pos
                    and (
                        p not in pos_ads1
                        or pos_ads1[p][0] - ads[0] <= 5
                        and pos_ads1[p][2] - ads[1] <= 5
                    )
                ]
            params["segment1"] = segs1
            params["segment2"] = segs2
            params["segment_total"] = segs_total
            self.logger.debug(f"Depth for Gene1 {segs1}")
            self.logger.debug(f"Depth for Gene2 {segs2}")
            self.logger.debug(f"Depth for liftover in Gene1 {segs_total}")
            if var_ads:
                self.logger.debug(
                    f"AD stats {statistics.mean([a[0] / sum(a) if sum(a) > 0 else 0 for a in var_ads])}, {statistics.stdev([a[0] / sum(a) if sum(a) > 0 else 0 for a in var_ads])}"
                )
            else:
                self.logger.debug("AD stats: None")
            probs = {}
            best_probs = -1e9
            best_total = -1
            for total_cn in range(1, total_cn_ub + 1):
                probs1 = {}
                for cn1 in range(max(0, total_cn - cn2_ub), min(cn1_ub, total_cn) + 1):
                    probs1[cn1] = self.calc_loglikelihood(
                        var_ads, cn1, total_cn - cn1, params
                    )
                max_prob = max(probs1.values())
                if best_total == -1 or max_prob > best_probs:
                    best_total = total_cn
                    best_probs = max_prob
                probs[total_cn] = probs1
            probs1 = probs[best_total]
            if len(probs1) == 1 or max(probs1.values()) != min(probs1.values()):
                best_cn1 = max(probs1, key=probs1.get)
                best_cn2 = best_total - best_cn1
            else:
                best_cn1 = best_cn2 = -1
                print("Insufficient information to determine the copy numbers.")
            d["total_cn"] = best_total
            d["cn"] = best_cn1
            data[1]["orig"][region_id]["total_cn"] = best_total
            data[1]["orig"][region_id]["cn"] = best_cn2
            if region_id == "nodel":
                total_cn_ub = best_total
                cn1_ub = best_cn1
                cn2_ub = best_cn2
            if region_id in data[0]["liftover"]:
                data[0]["liftover"][region_id]["cn"] = best_total
            if region_id in data[1]["liftover"]:
                data[1]["liftover"][region_id]["cn"] = best_total
        for i, d in enumerate(data):
            for region_id, dd in d["liftover"].items():
                if "cn" not in dd:
                    dd["cn"] = data[1 - i]["orig"]["nodel"]["cn"]
        return data


class Phased_vcf:
    max_mismatch = 0.25

    def __init__(self, data, gene, read_data):
        self.data = data
        self.gene = gene
        self.read_data = read_data
        self.debug = False
        self.logger = get_logger(self.__class__.__name__, "INFO")

    @staticmethod
    def _annotate(v, dp=None):
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
    def split_homvar(v, cns):
        v1 = None
        gts = set(re.split(r"\||/", v.samples[0]["GT"]))
        if len(gts) == 1 and gts != set("0"):
            v1 = copy.deepcopy(v)
            v1.line = None
            v1.samples[0]["GT"] = "|".join(["1"] * cns[0])
            Phased_vcf._annotate(v1)
        return v1

    def proc_phased(self, phased, cns, matched, gene_id):
        def _proc_homvar():
            for v in phased:
                v1 = self.split_homvar(v, cns)
                if v1:
                    v1s.append(v1)

        v1s = []
        all_alleles = []
        for v in phased:
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
                label.append(self.gene.gene_names[gene_id])
                v1 = copy.deepcopy(v)
                v1.samples[0]["GT"] = "|".join(
                    ["1" if a in alt_alleles else "0" for a in gene1_alleles]
                )
                self._annotate(v1)
                v1s.append(v1)
            v.info["LABEL"] = ",".join(label)
        return phased, v1s

    @staticmethod
    def match(vars, v1s, v2s):
        matched = {}
        for pos, key, _, v in v1s:
            matched.setdefault(pos, [None, None, None])[0] = (
                key,
                re.split("\||/", v.samples[0]["GT"]),
                v,
            )
        for pos, key, _, v in v2s:
            matched.setdefault(pos, [None, None, None])[1] = (
                key,
                re.split("\||/", v.samples[0]["GT"]),
                v,
            )
        for key, v in vars:
            matched.setdefault(v.pos, [None, None, None])[2] = (
                key,
                re.split("\||/", v.samples[0]["GT"]),
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
            if (
                v1
                and v1[1].count("1") == v[1].count("1")
                or v2
                and v2[1].count("1") == v[1].count("1")
            ):
                keep[pos] = vs
        return keep

    @staticmethod
    def save_vars(vcf_template, out_vcf, vars, extra_header=None, sort_vars=False):
        in_vcf = vcflib.VCF(vcf_template)
        out_vcf = vcflib.VCF(out_vcf, "w")
        out_vcf.copy_header(in_vcf, update=extra_header)
        out_vcf.emit_header()
        if sort_vars:
            vars_list = [(v.pos, i, v) for i, v in enumerate(vars)]
            vars = [v[2] for v in sorted(vars_list)]
        for v in vars:
            out_vcf.emit(v)
        out_vcf.close()
        in_vcf.close()

    def concat_regional_vcf(self):
        for seq_key in self.read_data:
            params = self.read_data[seq_key]["params"]
            for i, d in enumerate(self.data):
                for bam_key, dd in d.items():
                    vcfs = []
                    for region_id, reg_data in dd.items():
                        if region_id != "all" and "result_vcf" in reg_data[seq_key]:
                            vcfs.append(reg_data[seq_key]["result_vcf"])
                    if vcfs:
                        out_fname = f"{params['prefix']}.{self.gene.gene_names[i]}.{bam_key}.vcf.gz"
                        self.concat(vcfs, out_fname)
                        dd.setdefault("all", {})[seq_key] = {
                            "result_vcf": out_fname,
                            "region": self.gene.gene_regions[i],
                        }

    def concat(self, in_vcfs, out_fname):
        vcfs = [vcflib.VCF(vcf) for vcf in in_vcfs]
        out_vcf = vcflib.VCF(out_fname, "w")
        out_vcf.copy_header(vcfs[0])
        out_vcf.emit_header()
        vars = []
        cn = [-1, -1]
        for i, vcf in enumerate(vcfs):
            for v in vcf:
                if cn[i] < 0:
                    cn[i] = len(re.split("\||/", v.samples[0]["GT"]))
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

    def process(self):
        label_header = (
            '##INFO=<ID=LABEL,Number=1,Type=String,Description="LABEL for GENE1 or GENE2">',
        )
        for i, d in enumerate(self.data):
            result_vcfs = {}
            result_liftover_vcfs = {}
            for region_id, dd in d["orig"].items():
                for seq_key, data in self.read_data.items():
                    seq_data = dd[seq_key]
                    params = data["params"]
                    if seq_key not in result_vcfs:
                        result_vcfs[seq_key] = []
                        result_liftover_vcfs[seq_key] = []
                    liftover_vcf = d["liftover"][region_id][seq_key].get("phased_vcf")
                    if not liftover_vcf:
                        continue
                    vcf1 = seq_data.get("phased_vcf")
                    vcf2 = (
                        self.data[1 - i]["orig"]
                        .get(region_id, self.data[1 - i]["orig"]["nodel"])[seq_key]
                        .get("phased_vcf")
                    )
                    cn1 = dd["cn"]
                    cn2 = dd["total_cn"] - cn1
                    if cn1 == 0:
                        continue
                    run_params = {
                        "region": dd["region"],
                        "orig_region": dd[seq_key]["region"],
                        "gene_id": i,
                        "seq": seq_key,
                        "cns": (cn1, cn2),
                        "cn_diff": d["orig"]["nodel"]["cn"] - cn1,
                        "mismatch": params.get("max_mismatch", self.max_mismatch),
                    }
                    out_liftover, out_vcf1 = self.assign(
                        liftover_vcf, vcf1, vcf2, run_params
                    )
                    reg_name = "nodel" if region_id == "nodel" else "del" + region_id
                    out_liftover_fname = f"{params['prefix']}.{self.gene.gene_names[i]}.{reg_name}.liftover.labeled.vcf.gz"
                    out_vcf1_fname = f"{params['prefix']}.{self.gene.gene_names[i]}.{reg_name}.vcf.gz"
                    # result_vcfs[seq_key].append(out_vcf1_fname)
                    # result_liftover_vcfs[seq_key].append(out_liftover_fname)
                    vcf_tmpl = vcf1 if vcf1 else (vcf2 if vcf2 else liftover_vcf)
                    self.save_vars(vcf_tmpl, out_vcf1_fname, out_vcf1)
                    self.save_vars(
                        liftover_vcf,
                        out_liftover_fname,
                        out_liftover,
                        extra_header=label_header,
                    )
                    d["liftover"][region_id][seq_key]["result_vcf"] = out_liftover_fname
                    seq_data["result_vcf"] = out_vcf1_fname
            # concat the vcf for either long and short reads from different regions
        self.concat_regional_vcf()
        for i, d in enumerate(self.data):
            if d["orig"]["nodel"]["cn"] == 0:
                continue
            # combine different seq_keys
            out_vcf1_fname = f"{params['outdir']}/{params['sample_name']}.{self.gene.gene_names[i]}.result.vcf.gz"
            if len(self.read_data) > 1:
                params = self.read_data["short_read"]["params"]
                long_vcf1 = vcflib.VCF(d["orig"]["all"]["long_read"]["result_vcf"])
                long_vcf2 = vcflib.VCF(
                    self.data[1 - i]["orig"]["all"]["long_read"]["result_vcf"]
                )
                short_vcf = vcflib.VCF(d["liftover"]["all"]["short_read"]["result_vcf"])
                short_vcf1 = vcflib.VCF(d["orig"]["all"]["short_read"]["result_vcf"])
                short_vcf2 = vcflib.VCF(
                    self.data[1 - i]["orig"]["all"]["short_read"]["result_vcf"]
                )
                convert_fn = (
                    self.gene.mapping.convert21
                    if i == 0
                    else self.gene.mapping.convert12
                )
                converted_long_vcf2 = convert_fn([v for v in long_vcf2], ref_diff=False)
                converted_short_vcf2 = convert_fn(
                    [v for v in short_vcf2], ref_diff=False
                )
                out_vcf1 = []
                for v1 in self.vcf_grouper(short_vcf, short_vcf1, long_vcf1):
                    out_vcf1.append(v1)
                self.save_vars(
                    d["orig"]["all"]["short_read"]["result_vcf"],
                    out_vcf1_fname,
                    out_vcf1,
                )
                del d["orig"]["all"]["short_read"]["result_vcf"]
                del d["orig"]["all"]["long_read"]["result_vcf"]
            else:
                shutil.copy(
                    d["orig"]["all"]["short_read"]["result_vcf"], out_vcf1_fname
                )
                ext = ".tbi" if out_vcf1_fname.endswith(".gz") else ".idx"
                shutil.copy(
                    d["orig"]["all"]["short_read"]["result_vcf"] + ext,
                    out_vcf1_fname + ext,
                )
                del d["orig"]["all"]["short_read"]["result_vcf"]
            d["orig"]["all"]["result_vcf"] = out_vcf1_fname
        return self.data

    def assign(self, liftover_vcf, in_vcf1, in_vcf2, params):
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
        map_fn = (
            self.gene.mapping.interval12
            if gene_id == 0
            else self.gene.mapping.interval21
        )
        orig_region2 = IntervalList(region=map_fn(orig_region)) if orig_region else None
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
            convert_fn = (
                self.gene.mapping.convert21
                if gene_id == 0
                else self.gene.mapping.convert12
            )
            for v2 in convert_fn(to_convert, ref_diff=False):
                heapq.heappush(all_v2s, (v2.pos, f"{v2.ref}_{v2.alt[0]}", 1, v2))
        cns = params["cns"]
        if len(cns) != 2:
            print(
                "ERROR: exactly two copy numbers required for the gene and its paralog seperated by comma. Example: 2,3",
                file=sys.stderr,
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
                        if self.debug:
                            print(
                                f"{self.gene.gene_names[gene_id]}: {v.chrom}:{v.pos} {v.ref}-{','.join(v.alt)}: original GT {last[3].samples[0]['GT']} replaced with {cur[3].samples[0]['GT']}."
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
        if diff:
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
