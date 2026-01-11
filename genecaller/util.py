import vcflib
import io
import bisect
import collections
import re
import heapq
import copy
import sys
import pysam
import gzip
from typing import List, Optional, Dict, Tuple, Union, Any, Callable
import subprocess
import shutil
from packaging.version import Version
from importlib import resources
import os

from .logging import get_logger

logger = get_logger(__name__)


def get_cmd(
    cmd: str,
    version: Optional[Version] = None,
    return_version: bool = False,
) -> Optional[Union[str, Tuple[str, str]]]:
    cmd_list: List[str] = cmd.split()
    exec_file = shutil.which(cmd_list[0])
    if not exec_file:
        logger.error("Error: no '%s' found in the PATH", cmd)
        return None

    if version is None and not return_version:
        return exec_file

    cmd_list.append("--version")
    cmd_version_str = (
        subprocess.check_output(cmd_list).decode("utf-8", "ignore").strip()
    )
    if cmd_list[0] == "sentieon":
        cmd_version_str = cmd_version_str.split("-")[-1]
    elif cmd_list[0] == "pbsv":
        cmd_version_str = cmd_version_str.split(" ")[1]
    elif cmd_list[0] == "hificnv":
        cmd_version_str = cmd_version_str.split(" ")[1].split("-")[0]
    else:
        cmd_version_str = cmd_version_str.split("\n")[0].split()[-1].split("-")[0]

    if version is not None:
        cmd_version = Version(cmd_version_str)
        if cmd_version < version:
            logger.error(
                "Error: the pipeline requires %s version '%s' or later "
                "but %s '%s' was found in the PATH",
                cmd,
                version,
                cmd,
                cmd_version,
            )
            return None

    if return_version:
        return exec_file, cmd_version_str
    return exec_file


def get_software_versions(calling_min_versions: Dict) -> Dict[str, str]:
    versions = {}
    for cmd_full in calling_min_versions.keys():
        tool_name = cmd_full.split()[0]
        try:
            result = get_cmd(cmd_full, return_version=True)
            if result:
                _, version_str = result
                versions[tool_name] = version_str
            else:
                versions[tool_name] = "unknown"
        except Exception as e:
            logger.warning(f"Could not determine version for {tool_name}: {e}")
            versions[tool_name] = "unknown"
    return versions


def get_data_file(input_filename: str) -> Optional[str]:
    """
    Finds a specific file. First processes bundle logic to get the correct path,
    then checks if the path exists as-is, then looks within the package's 'data' directory.
    Returns the path to the file.
    """

    # Handle bundle logic on the original path: only when second-to-last component ends with .bundle
    parts = input_filename.split("/")
    if not parts[-1]:
        parts = parts[:-1]
    bundle = False
    bundle_index = -1
    for i, p in enumerate(parts):
        if p.endswith(".bundle") and i < len(parts) - 1:
            bundle = True
            bundle_index = i
            break

    if bundle:
        # For bundle files, the actual path is up to the bundle directory
        actual_path = "/".join(parts[: bundle_index + 1])
        model_file = "/".join(parts[bundle_index + 1 :])
    else:
        actual_path = input_filename
        model_file = None

    # Check if the actual path exists as a regular file path
    if os.path.exists(actual_path):
        if bundle and model_file:
            return os.path.join(actual_path, model_file)
        else:
            return actual_path

    # If not found, look in the package data directory
    # Strip leading "data/" for package lookup
    package_parts = actual_path.split("/")
    if package_parts[0] == "data":
        package_parts = package_parts[1:]
        actual_path = "/".join(package_parts)

    try:
        data_dir_traversable = resources.files("genecaller").joinpath("data")
        file_path_traversable = data_dir_traversable.joinpath(actual_path)

        with resources.as_file(file_path_traversable) as concrete_file_path:
            # concrete_file_path is a pathlib.Path object
            # This path is guaranteed to exist on the filesystem.
            # If the package is zipped, the file will be temporarily extracted.
            if bundle and model_file:
                return str(concrete_file_path.joinpath(model_file))
            return str(concrete_file_path)

    except FileNotFoundError:
        logger.error(f"Error: Resource {actual_path} not found in genecaller.data")
        return None
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        return None


class IntervalList(object):
    def __init__(self, region: str = "", path: str = "") -> None:
        if path:
            self.regions = self.load_bed(path)
        else:
            self.regions = self.load_region(region)

    def size(self) -> int:
        total = 0
        for _, cur_regions in self.regions.items():
            i = 0
            while i < len(cur_regions):
                start, end = cur_regions[i], cur_regions[i + 1]
                total += end - start
                i += 2
        return total

    def load_region(self, input_region: str) -> Dict[str, List[int]]:
        if not input_region:
            return {}
        regions = {}
        for rgn in input_region.split(","):
            chr, parts = rgn.split(":")
            start, end = [int(p) for p in parts.split("-")]
            regions.setdefault(chr, []).append([start - 1, end])
        for chrom in regions:
            v = []
            s0 = e0 = None
            for s, e in sorted(regions[chrom]):
                if e0 is None:
                    s0 = s
                    e0 = e
                elif s > e0:
                    v.extend((s0, e0))
                    s0 = s
                    e0 = e
                else:
                    e0 = e
            if e0 is not None:
                v.extend((s0, e0))
            regions[chrom] = v
        return regions

    def load_bed(self, path: str) -> Dict[str, List[int]]:
        regions = {}
        if path.endswith(".gz"):
            fp = vcflib.bgzf.open(path, "rb")
        else:
            fp = io.open(path, "rb")
        for line in fp:
            if line.startswith(b"#") or line.startswith(b"track"):
                continue
            cols = line.rstrip().decode().split("\t")
            if len(cols) < 3:
                continue
            chrom, start, end = cols[:3]
            regions.setdefault(chrom, []).append([int(start), int(end)])
        fp.close()
        for chrom in regions:
            v = []
            s0 = e0 = None
            for s, e in sorted(regions[chrom]):
                if e0 is None:
                    s0 = s
                    e0 = e
                elif s > e0:
                    v.extend((s0, e0))
                    s0 = s
                    e0 = e
                else:
                    e0 = e
            if e0 is not None:
                v.extend((s0, e0))
            regions[chrom] = v
        return regions

    def get(self, c: str, s: int, e: int) -> List[Tuple[int, int]]:
        r = self.regions.get(c)
        v = []
        if r is not None:
            i = bisect.bisect_right(r, s)
            if i % 2:
                v.append((r[i - 1], r[i]))
                i += 1
            while i < len(r) and r[i] < e:
                v.append((r[i], r[i + 1]))
                i += 2
        return v

    def intersection(self, c: str, s: int, e: int) -> List[Tuple[int, int]]:
        regs = self.get(c, s, e)
        itsc = []
        for reg_start, reg_end in regs:
            reg_start = max(reg_start, s)
            reg_end = min(reg_end, e)
            itsc.append((reg_start, reg_end))
        return itsc

    def __contains__(self, cs: Union[str, Tuple[str, int]]) -> bool:
        if isinstance(cs, str):
            return cs in self.regions
        c, s = cs
        r = self.regions.get(c)
        return bool(r and bisect.bisect_right(r, s) % 2 != 0)

    # region can be {c: [(s, e)]}, or [(c, s, e)] or c:s-e,c:s-e
    @staticmethod
    def parse_region(
        region: Union[
            "IntervalList",
            Dict[str, List[Tuple[int, int]]],
            List[Tuple[str, int, int]],
            str,
        ],
    ) -> Dict[str, List[Tuple[int, int]]]:
        if isinstance(region, IntervalList):
            res = {}
            for chr, locs in region.regions.items():
                cur_regions = []
                for i in range(len(locs) // 2):
                    cur_regions.append((locs[i * 2], locs[i * 2 + 1]))
                res[chr] = cur_regions
            return res
        if isinstance(region, dict):
            return region
        if isinstance(region, list):
            res = {}
            for c, s, e in region:
                res.setdefault(c, []).append([s, e])
            return res
        res = {}
        for rgn in region.split(","):
            chr, parts = rgn.split(":")
            start, end = [int(p) for p in parts.split("-")]
            res.setdefault(chr, []).append((start - 1, end))
        return res

    @staticmethod
    def merge_regions(regino_list: List[Tuple[int, int]]) -> List[int]:
        v = []
        s0 = e0 = None
        for s, e in sorted(regino_list):
            if e0 is None:
                s0 = s
                e0 = e
            elif s > e0:
                v.extend((s0, e0))
                s0 = s
                e0 = e
            else:
                e0 = e
        if e0 is not None:
            v.extend((s0, e0))
        return v

    def intersect(
        self, region: Union[str, "IntervalList", Dict, List]
    ) -> "IntervalList":
        if not region:
            return self
        result_interval = copy.deepcopy(self)
        regions = {}
        input_region = self.parse_region(region)
        for chr, rgns in input_region.items():
            if chr not in self.regions:
                continue
            res = []
            for r in rgns:
                rs = self.get(chr, r[0], r[1])
                if rs:
                    if rs[0][0] < r[0]:
                        rs[0] = (r[0], rs[0][1])
                    if rs[-1][1] > r[1]:
                        rs[-1] = (rs[-1][0], r[1])
                res += rs
            if res:
                regions[chr] = self.merge_regions(res)
        result_interval.regions = regions
        return result_interval

    def union(self, region: Union[str, "IntervalList", Dict, List]) -> "IntervalList":
        if not region:
            return self
        input_regions = self.parse_region(region)
        regions = {}

        # Process all chromosomes from both self and input regions
        all_chrs = set(self.regions.keys()) | set(input_regions.keys())

        for chr in all_chrs:
            rgn = input_regions.get(chr, [])
            if chr in self.regions:
                i = 0
                cur_regions = self.regions[chr]
                while i < len(cur_regions):
                    rgn.append((cur_regions[i], cur_regions[i + 1]))
                    i += 2
            regions[chr] = self.merge_regions(rgn)

        result_interval = copy.deepcopy(self)
        result_interval.regions = regions
        return result_interval

    def subtract(
        self, region: Union[str, "IntervalList", Dict, List]
    ) -> "IntervalList":
        if not region:
            return self
        result_interval = copy.deepcopy(self)
        result_interval.regions = {}
        input_region = self.parse_region(region)
        for chr, rgns in input_region.items():
            if chr not in self.regions:
                continue
            i = 0
            rs = self.regions[chr]
            cur_regions = []
            while i < len(rs):
                cur_regions.append((rs[i], rs[i + 1]))
                i += 2
            for rr in rgns:
                subtract_intervals = self.get(chr, rr[0], rr[1])
                if subtract_intervals:
                    if subtract_intervals[0][0] < rr[0]:
                        subtract_intervals[0] = (rr[0], subtract_intervals[0][1])
                    if subtract_intervals[-1][1] > rr[1]:
                        subtract_intervals[-1] = (subtract_intervals[-1][0], rr[1])

                    subtracted = []
                    for current_region in cur_regions:
                        remaining_parts = [current_region]
                        for subtract_interval in subtract_intervals:
                            new_remaining = []
                            for part in remaining_parts:
                                if (
                                    part[1] <= subtract_interval[0]
                                    or part[0] >= subtract_interval[1]
                                ):
                                    new_remaining.append(part)
                                else:
                                    if part[0] < subtract_interval[0]:
                                        new_remaining.append(
                                            (part[0], subtract_interval[0])
                                        )
                                    if part[1] > subtract_interval[1]:
                                        new_remaining.append(
                                            (subtract_interval[1], part[1])
                                        )
                            remaining_parts = new_remaining
                        subtracted.extend(remaining_parts)
                    cur_regions = subtracted
            result_interval.regions[chr] = self.merge_regions(cur_regions)
        return result_interval

    def pad(self, padding: int) -> "IntervalList":
        for chr, regions in self.regions.items():
            cur_regions = []
            i = 0
            while i < len(regions):
                cur_regions.append(
                    [max(0, regions[i] - padding), regions[i + 1] + padding]
                )
                i += 2
            self.regions[chr] = self.merge_regions(cur_regions)
        return self

    def to_region_str(self, chrs: Optional[set] = None) -> str:
        res = []
        for chr, regions in self.regions.items():
            i = 0
            if chrs and chr not in chrs:
                continue
            cur_regions = self.regions[chr]
            while i < len(cur_regions):
                start, end = cur_regions[i], cur_regions[i + 1]
                res.append(self.interval_str(chr, start, end))
                i += 2
        return ",".join(res)

    @staticmethod
    def interval_str(chr: str, s: int, e: int) -> str:
        return f"{chr}:{s + 1}-{e}"


Contig = collections.namedtuple("Contig", "length offset width skip")


class Reference(object):
    def __init__(self, path: str) -> None:
        self.path = path
        self.index = collections.OrderedDict()
        with io.open(self.path + ".fai", "rb") as fp:
            for line in fp:
                flds = line.rstrip().decode().split("\t")
                if len(flds) != 5:
                    raise RuntimeError("Corrupted index")
                self.index[flds[0]] = Contig._make(map(int, flds[1:]))
        self.fp = io.open(path, "rb")

    def __getstate__(self) -> Dict[str, Any]:
        odict = self.__dict__.copy()
        odict["fp"] = self.fp.tell()
        return odict

    def __setstate__(self, ndict: Dict[str, Any]) -> None:
        path = ndict["path"]
        fp = io.open(path, "rb")
        fp.seek(ndict.pop("fp"))
        ndict["fp"] = fp
        self.__dict__.update(ndict)

    def get(self, c: str, s: int, e: int) -> str:
        ci = self.index[c]
        e = min(e, ci.length)
        if s > e:
            raise ValueError
        seq = b""
        while s < e:
            o = ci.offset + s // ci.width * ci.skip + s % ci.width
            n = ci.width - s % ci.width
            n = min(n, e - s)
            self.fp.seek(o)
            seq += self.fp.read(n)
            s += n
        return seq.decode()

    def __iter__(self) -> Any:
        return iter(self.index.items())


TRANS_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(seq: str) -> str:
    return seq.translate(TRANS_TABLE)[::-1]


def complement(seq: str) -> str:
    return seq.translate(TRANS_TABLE)


MappingBlock = collections.namedtuple("MappingBlock", "s1, e1, s2, e2, r1, r2")


class HeapItem:
    __slots__ = ("start", "end", "is_variant", "data", "uid")

    def __init__(self, s, e, is_v, d, u):
        self.start, self.end, self.is_variant, self.data, self.uid = s, e, is_v, d, u

    def __lt__(self, o):
        if self.start != o.start:
            return self.start < o.start
        if self.is_variant != o.is_variant:
            return not self.is_variant
        if self.end != o.end:
            return self.end < o.end
        return self.uid < o.uid


class GeneMapping:
    def __init__(
        self, mapping_file: str, combine: bool = False, ref_genome: Optional[Any] = None
    ) -> None:
        self.intervals1, self.intervals2 = [], []
        self.reverse, self.ref_genome = False, ref_genome
        self.g1tog2, self.g2tog1 = {}, {}

        with gzip.open(mapping_file, "rt") as f:
            init_s2 = 0
            for line in f:
                flds = line.strip().split("\t")
                if len(flds) < 8:
                    continue
                c1, c2 = map(int, flds[:2])
                r1, r2 = flds[2:4]
                s1, e1, s2, e2 = map(int, flds[4:])

                if not init_s2:
                    init_s2 = s2
                elif init_s2 > 0:
                    self.reverse = s2 < init_s2
                    init_s2 = -1
                if c1 == s1 and c2 == s2:
                    blk = MappingBlock(s1, e1, s2, e2, r1, r2)
                    self.intervals1.append(blk)
                    self.intervals2.append(blk)
                self.g1tog2[c1] = c2
                self.g2tog1[c2] = c1
                if combine:
                    self.g1tog2[c1] = c2
                    self.g2tog1[c2] = c1

        self.intervals1.sort(key=lambda x: x.s1)
        self.keys1 = [x.s1 for x in self.intervals1]
        self.intervals2.sort(key=lambda x: x.s2)
        self.keys2 = [x.s2 for x in self.intervals2]

    def convert21(
        self,
        variants: List[Any],
        ref_diff: bool = True,
        annotate_fn: Optional[Callable] = None,
    ):
        return self.convert_vars(variants, True, ref_diff, annotate_fn)

    def convert12(
        self,
        variants: List[Any],
        ref_diff: bool = True,
        annotate_fn: Optional[Callable] = None,
    ):
        return self.convert_vars(variants, False, ref_diff, annotate_fn)

    def convert_vars(
        self,
        variants: List[Any],
        to_gene1: bool,
        ref_diff: bool,
        annotate_fn: Optional[Callable] = None,
    ):
        blocks, use_s2 = (
            (self.intervals2, True) if to_gene1 else (self.intervals1, False)
        )
        pq, uid = [], 0
        cn = 2
        if variants:
            cn = len(re.split(r"\||/", variants[0].samples[0]["GT"]))
            chr = variants[0].chrom

        for v in variants:
            heapq.heappush(pq, HeapItem(v.pos, v.pos + len(v.ref), True, v, uid))
            uid += 1

        if ref_diff:
            for b in blocks:
                if b.r1 != b.r2:
                    s, e = (b.s2, b.e2) if use_s2 else (b.s1, b.e1)
                    heapq.heappush(pq, HeapItem(s, e, False, b, uid))
                    uid += 1

        proc, max_end = [], -1
        while pq:
            item = heapq.heappop(pq)
            if proc and item.start >= max_end:
                yield from self._process_block(
                    proc, to_gene1, ref_diff, annotate_fn, (cn, chr)
                )
                proc, max_end = [], -1
            proc.append(item)
            max_end = max(max_end, item.end)

        if proc:
            yield from self._process_block(
                proc, to_gene1, ref_diff, annotate_fn, (cn, chr)
            )

    def _process_block(
        self,
        items: List[HeapItem],
        to_g1: bool,
        rd: bool,
        ann: Optional[Callable],
        params: Any,
    ):
        cn, chr = params
        vars_in = [x.data for x in items if x.is_variant]
        if not vars_in:
            if not rd:
                return
            elif to_g1:  # all ref-diff variants can be output as is
                for v in items:
                    var = self._mk_var(
                        chr, v.data.s1, v.data.r1, [v.data.r2], None, [1] * cn
                    )
                    yield ann(var) if ann else var
                return

        min_s, max_e = min(x.start for x in items), max(x.end for x in items)
        intervals, keys = (
            (self.intervals2, self.keys2) if to_g1 else (self.intervals1, self.keys1)
        )
        if not to_g1 and self.reverse:
            max_e += 1
        bg = self._get_bg(min_s, max_e, intervals, keys, to_g1)
        if not bg:
            return

        sb_s = (
            min(bg, key=lambda b: b.s2).s2 if to_g1 else min(bg, key=lambda b: b.s1).s1
        )
        r2_seq = "".join(b.r2 for b in bg)
        r1_seq = "".join(b.r1 for b in bg)
        in_ctx, tgt_ctx = (r2_seq, r1_seq) if to_g1 else (r1_seq, r2_seq)
        if self.reverse and to_g1:
            in_ctx = rc(in_ctx)  # ref with diff variants
        if self.reverse and not to_g1:
            tgt_ctx = rc(tgt_ctx)
        src_v = None
        var_types = []
        alleles_pos = []
        for _ in range(cn):
            alleles_pos.append({p: b for p, b in enumerate(in_ctx, sb_s)})
            var_types.append([])
        stars = []
        if vars_in:
            vars_in.sort(key=lambda v: v.pos)
            src_v = vars_in[-1]
            for v in vars_in:
                gt = list(map(int, re.split(r"\||/", v.samples[0]["GT"])))
                for i, g in enumerate(gt):
                    if not g:  # reference allele
                        continue
                    a = v.alt[g - 1]
                    if a == "*":
                        stars.append(i)
                        continue
                    if len(a) == len(v.ref):  # snp
                        for ii, aa in enumerate(a):
                            if v.pos + ii in alleles_pos[i]:
                                alleles_pos[i][v.pos + ii] = aa
                    elif len(a) > len(v.ref):  # ins
                        if v.pos in alleles_pos[i]:
                            alleles_pos[i][v.pos] += a[1:]
                    else:  # del. remove everything after the first base
                        for ii in range(1, len(v.ref)):
                            if v.pos + ii in alleles_pos[i]:
                                del alleles_pos[i][v.pos + ii]
        final_in = ["".join([ap[k] for k in sorted(ap.keys())]) for ap in alleles_pos]
        lifted = [rc(s) if self.reverse else s for s in final_in]
        new_ref = tgt_ctx
        new_alts = list(set([a for a in lifted if a != new_ref]))
        if not new_alts:
            return
        gt = []
        for seq in lifted:
            if seq == new_ref:
                gt.append(0)
            else:
                gt.append(new_alts.index(seq) + 1)

        fin_ref, fin_alts = new_ref, new_alts
        if to_g1:
            t_pos = min(bg, key=lambda b: b.s1).s1
        else:
            t_pos = min(bg, key=lambda b: b.s2).s2

        while True:
            op = False
            # decompose multiple snps in a row
            if (
                max([len(a) for a in fin_alts])
                == min([len(a) for a in fin_alts])
                == len(fin_ref)
            ):
                if len(fin_ref) > 1:
                    for k, r in enumerate(fin_ref):
                        all_alts = [a[k] for a in fin_alts]
                        if any([r != a for a in all_alts]):
                            v = self._mk_var(chr, t_pos + k, r, all_alts, src_v, gt)
                            yield ann(v) if ann else v
                    return
            elif any(
                [a[0] != fin_ref[0] for a in fin_alts]
            ):  # add anchor base for INDEL
                t_pos -= 1
                if t_pos in self.keys1:
                    anc = self.intervals1[self.keys1.index(t_pos)].r1
                elif self.ref_genome:
                    anc = self.ref_genome.get(chr, t_pos, t_pos + 1)[0]
                else:
                    raise Exception(
                        f"Failed to get anchor base for indel at {chr}:{t_pos}"
                    )
                fin_alts = [anc + a for a in fin_alts]
                fin_ref = anc + fin_ref
                op = True

            # remove trailing bases
            while (
                len(fin_ref) > 1
                and all(len(a) > 1 for a in fin_alts)
                and all(a[-1] == fin_ref[-1] for a in fin_alts)
            ):
                fin_ref = fin_ref[:-1]
                fin_alts = [a[:-1] for a in fin_alts]
                op = True

            # remove leading bases
            while (
                len(fin_ref) > 1
                and all(len(a) > 1 for a in fin_alts)
                and all(a[:2] == fin_ref[:2] for a in fin_alts)
            ):
                fin_ref = fin_ref[1:]
                fin_alts = [a[1:] for a in fin_alts]
                t_pos += 1
                op = True

            # Remove trailing SNPs
            while (
                all([fin_ref[0] == a[0] for a in fin_alts])
                and len(fin_ref) > 1
                and min([len(a) for a in fin_alts]) > 1
            ):
                v = self._mk_var(
                    chr,
                    t_pos + len(fin_ref) - 1,
                    fin_ref[-1],
                    [a[-1] for a in fin_alts],
                    src_v,
                    gt,
                )
                yield ann(v) if ann else v
                fin_ref = fin_ref[:-1]
                fin_alts = [a[:-1] for a in fin_alts]

            if not op:
                break

        if len(fin_ref) == 0 and all(len(a) == 0 for a in fin_alts):
            return

        if self.ref_genome and (
            len(fin_ref) == 0 or any(len(a) == 0 for a in fin_alts)
        ):
            anc = self.ref_genome.get(chr, t_pos - 1, t_pos)
            fin_ref = anc + fin_ref
            fin_alts = [anc + a for a in fin_alts]
            t_pos -= 1

        if stars:
            stars = set(stars)
            for star in stars:
                if gt[star] == 0:
                    new_alts.append("*")
                    gt[star] = len(new_alts)
                else:
                    new_alts[gt[star] - 1] = "*"
                if len(new_alts) == 1:
                    return

        v = self._mk_var(chr, t_pos, fin_ref, fin_alts, src_v, gt)
        yield ann(v) if ann else v

    def _mk_var(self, c, p, r, a, src, gt):
        if src:
            nv = copy.deepcopy(src)
            nv.chrom, nv.pos, nv.id, nv.ref, nv.alt = c, p, ".", r, a
            nv.end = p + len(r)
            if len(gt) > 1:
                nv.samples[0]["GT"] = nv.samples[0]["GT"][1].join(map(str, gt))
            else:
                nv.samples[0]["GT"] = str(gt[0])
            nv.line = None
            return nv
        else:
            s = [{"GT": "|".join(map(str, gt))}]
            return vcflib.Variant(c, p, ".", r, a, None, [], {}, s, p + len(r), None)

    def _get_bg(self, s, e, ints, keys, use_s2):
        idx = bisect.bisect_right(keys, s) - 1
        idx = 0 if idx < 0 else idx
        res = []
        for i in range(idx, len(ints)):
            b = ints[i]
            bs, be = (b.s2, b.e2) if use_s2 else (b.s1, b.e1)
            if bs >= e:
                break
            if be > s and bs < e:
                res.append(b)
        res.sort(key=lambda x: x.s1)  # sorted in the original
        return res

    def get(self, s, direction="auto"):
        if direction == "forward":
            mapping = self.g1tog2
        elif direction == "backward":
            mapping = self.g2tog1
        else:  # "auto detection"
            if s in self.g1tog2 or s + 5 in self.g1tog2 or s - 5 in self.g1tog2:
                mapping = self.g1tog2
            elif s in self.g2tog1 or s + 5 in self.g2tog1 or s - 5 in self.g2tog1:
                mapping = self.g2tog1
            else:
                raise ValueError(f"position {s} not in mapping")
        if s in mapping:
            return mapping[s]
        if s - 5 in mapping:
            return mapping[s - 5] + 5
        if s + 5 in mapping:
            return mapping[s + 5] - 5
        raise ValueError(f"position {s} not in mapping")

    def interval(self, region: str, direction="auto") -> str:
        if not region:
            return ""

        def _func(s):
            return self.get(s, direction)

        out_regions = []
        for chr, regions in IntervalList.parse_region(region).items():
            for s, e in regions:
                s1, e1, m1 = _func(s), _func(e - 1), _func((s + e) // 2)
                if self.reverse and s1 < m1 or not self.reverse and s1 > m1:
                    n, ss = 0, s
                    while True:
                        n, ss = n + 1, ss + 1
                        ss1 = _func(ss)
                        if self.reverse and ss1 >= m1 or not self.reverse and ss1 <= m1:
                            s1 = ss1 + n if self.reverse else ss1 - n
                            break
                if self.reverse and e1 > m1 or not self.reverse and e1 < m1:
                    n, ee = 0, e
                    while True:
                        n, ee = n + 1, ee - 1
                        ee1 = _func(ee)
                        if self.reverse and ee1 <= m1 or not self.reverse and ee1 >= m1:
                            e1 = ee1 - n if self.reverse else ee1 + n
                            break
                if self.reverse:
                    s1, e1 = e1, s1
                out_regions.append(IntervalList.interval_str(chr, s1, e1 + 1))
        return ",".join(out_regions)

    def interval12(self, region: str) -> str:
        return self.interval(region, "forward")

    def interval21(self, region: str) -> str:
        return self.interval(region, "backward")


def load_bam(in_bam: str, in_ref: Optional[str] = None) -> pysam.AlignmentFile:
    if in_bam.endswith(".cram"):
        if not in_ref:
            logger.error("No reference file provided for CRAM input")
            sys.exit(1)
        bamh = pysam.AlignmentFile(in_bam, "rc", reference_filename=in_ref)
    else:
        bamh = pysam.AlignmentFile(in_bam, "rb")
    return bamh
