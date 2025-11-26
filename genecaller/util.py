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


Mapping = collections.namedtuple("Mapping", "pos, m_pos, ref, m_ref")


class GeneMapping:
    def __init__(self, mapping_file: str, combine: bool = False) -> None:
        self.g1tog2 = {}
        self.g2tog1 = {}
        self.map12 = []
        self.map21 = []
        self.reverse = False
        with gzip.open(mapping_file, "rt") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue
                c1, c2, r1, r2, s1, e1, s2, e2 = (
                    [int(s) - 1 for s in fields[:2]]
                    + fields[2:4]
                    + [int(s) - 1 for s in fields[4:]]
                )
                self.reverse = e2 < s2  # type: ignore
                if c1 == s1 and c2 == s2:
                    if self.reverse:
                        s2 = e2 + 1  # type: ignore
                    self.map12.append(Mapping(s1, s2, r1, r2))
                    self.map21.append(Mapping(s2, s1, r2, r1))
                    if combine:
                        self.map21.append(Mapping(s1, s2, r1, r2))
                        self.map12.append(Mapping(s2, s1, r2, r1))
                self.g1tog2[c1] = c2
                self.g2tog1[c2] = c1
                if combine:
                    self.g1tog2[c2] = c1
                    self.g2tog1[c1] = c2
        self.map21 = sorted(self.map21)

    @staticmethod
    def rc(seq: str) -> str:
        rev = {"A": "T", "C": "G", "T": "A", "G": "C"}
        return "".join([rev[s] for s in seq[-1::-1]])

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
        func = lambda s: self.get(s, direction)

        out_regions = []
        for chr, regions in IntervalList.parse_region(region).items():
            for s, e in regions:
                s1, e1, m1 = func(s), func(e - 1), func((s + e) // 2)
                if self.reverse and s1 < m1 or not self.reverse and s1 > m1:
                    n, ss = 0, s
                    while True:
                        n, ss = n + 1, ss + 1
                        ss1 = func(ss)
                        if self.reverse and ss1 >= m1 or not self.reverse and ss1 <= m1:
                            s1 = ss1 + n if self.reverse else ss1 - n
                            break
                if self.reverse and e1 > m1 or not self.reverse and e1 < m1:
                    n, ee = 0, e
                    while True:
                        n, ee = n + 1, ee - 1
                        ee1 = func(ee)
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

    def convert_vars(
        self,
        v1s: List[Any],
        direction: str,
        ref_diff: bool = True,
        annotate_fn: Optional[Callable] = None,
    ) -> Any:
        if direction == "forward":
            yield from self.convert12(v1s, ref_diff, annotate_fn)
        elif direction == "backward":
            yield from self.convert21(v1s, ref_diff, annotate_fn)
        else:
            raise ValueError("No auto conversion direction is available")

    def convert12(
        self,
        v1s: List[Any],
        ref_diff: bool = True,
        annotate_fn: Optional[Callable] = None,
    ) -> Any:
        def _pop():
            pos, i, v = heapq.heappop(vars)
            vv = next(iters[i], None)
            if vv:
                heapq.heappush(vars, (vv.pos, i, vv))
            return pos, i, v

        if not v1s:
            return
        cn = len(re.split(r"\||/", v1s[0].samples[0]["GT"]))
        chr = v1s[0].chrom
        vars = []
        iters = [iter(v1s), iter(self.map12)]
        for i, it in enumerate(iters):
            v = next(it, None)
            if v:
                heapq.heappush(vars, (v.pos, i, v))
        while vars:
            to_proc = []
            last_loc = 0, 0
            while vars:
                pos, i, v = _pop()
                if v.pos in self.g2tog1:
                    yield v
                    continue
                if last_loc[0] and (v.pos >= last_loc[1] or v.pos < last_loc[0]):
                    heapq.heappush(vars, (pos, i, v))
                    break
                to_proc.append((i, v))
                last_loc = min(v.pos, last_loc[0]), max(last_loc[1], v.pos + len(v.ref))
            if len(to_proc) == 1 and to_proc[0][0] == 1:
                # if this is a diff site
                v = to_proc[0][1]
                if ref_diff and v.ref != v.m_ref:
                    v = (
                        chr,
                        v.m_pos,
                        "",
                        v.m_ref,
                        [v.ref],
                        None,
                        [],
                        {},
                        [{"GT": "/".join(["1"] * cn)}],
                        v.m_pos + len(v.m_ref),
                        None,
                    )
                    yield (
                        annotate_fn(vcflib.Variant(*v))
                        if annotate_fn
                        else vcflib.Variant(*v)
                    )
            else:
                mm = [v for i, v in to_proc if i == 1]
                if not mm:
                    if to_proc:
                        variant_pos = to_proc[0][1].pos
                        logger.warning(f"variant {variant_pos} cannot be matched")
                    continue
                r1 = "".join([m.ref for m in mm])
                r2 = "".join([m.m_ref for m in mm])
                pos = mm[0].pos
                init_pos = pos
                alleles = [""] * cn
                vv = [v for i, v in to_proc if i == 0]
                if not vv:
                    continue
                star = set()
                for v in vv:
                    gt = [int(g) for g in re.split(r"\||/", v.samples[0]["GT"])]
                    if "*" in v.alt:
                        star_idx = v.alt.index("*") + 1
                        star = star.union(
                            set([i for i, g in enumerate(gt) if g == star_idx])
                        )
                    if v.pos > pos:
                        alleles = [
                            a + r1[pos - init_pos : v.pos - init_pos] for a in alleles
                        ]
                    for i, g in enumerate(gt):
                        if g == 0:
                            alleles[i] += v.ref
                        else:
                            alleles[i] += v.alt[g - 1]
                    pos = v.pos + len(v.ref)
                if pos != mm[-1].pos + len(mm[-1].ref):
                    alleles = [a + r1[pos - init_pos :] for a in alleles]
                for i in star:
                    alleles[i] = "*"
                alts = list(set([a for a in alleles if a != r2]))
                if not alts:
                    continue
                gt = ["0" if a == r2 else str(alts.index(a) + 1) for a in alleles]
                v = copy.deepcopy(vv[0])
                v.pos = mm[0].m_pos
                min_len = min([len(a) for a in [r2] + alts if a != "*"])
                for i in range(min_len - 1):
                    last = set([a[-1] for a in [r2] + alts if a != "*"])
                    if len(last) == 1:  # trim the last bit
                        r2 = r2[:-1]
                        alts = [a[:-1] if a != "*" else a for a in alts]
                    else:
                        break
                min_len = min([len(a) for a in [r2] + alts if a != "*"])
                for i in range(min_len - 1):
                    first = set([a[0] for a in [r2] + alts if a != "*"])
                    if len(first) == 1:  # trim the last bit
                        r2 = r2[1:]
                        alts = [a[1:] if a != "*" else a for a in alts]
                        v.pos += 1
                    else:
                        break
                v.ref = r2
                v.alt = alts
                sep = "/" if "/" in v.samples[0]["GT"] else "|"
                v.samples[0]["GT"] = sep.join(gt)
                v.line = None
                yield annotate_fn(v) if annotate_fn else v

    def convert21(
        self,
        v1s: List[Any],
        ref_diff: bool = True,
        annotate_fn: Optional[Callable] = None,
    ) -> Any:
        def _pop():
            pos, i, v = heapq.heappop(vars)
            vv = next(iters[i], None)
            if vv:
                heapq.heappush(vars, (vv.pos, i, vv))
            return pos, i, v

        if not v1s:
            return
        cn = len(re.split(r"\||/", v1s[0].samples[0]["GT"]))
        chr = v1s[0].chrom
        vars = []
        iters = [iter(v1s), iter(self.map21)]
        for i, it in enumerate(iters):
            v = next(it, None)
            if v:
                heapq.heappush(vars, (v.pos, i, v))
        while vars:
            to_proc = []
            last_end = 0
            while vars:
                pos, i, v = _pop()
                if v.pos in self.g1tog2:
                    yield v
                    continue
                if last_end and v.pos >= last_end:
                    heapq.heappush(vars, (pos, i, v))
                    break
                to_proc.append((i, v))
                last_end = max(last_end, v.pos + len(v.ref))
            if len(to_proc) == 1 and to_proc[0][0] == 1:
                # if this is a diff site
                v = to_proc[0][1]
                if ref_diff and v.ref != v.m_ref:
                    v = (
                        chr,
                        v.m_pos,
                        "",
                        v.m_ref,
                        [v.ref],
                        None,
                        [],
                        {},
                        [{"GT": "/".join(["1"] * cn)}],
                        v.m_pos + len(v.m_ref),
                        None,
                    )
                    yield (
                        annotate_fn(vcflib.Variant(*v))
                        if annotate_fn
                        else vcflib.Variant(*v)
                    )
            else:
                mm = [v for i, v in to_proc if i == 1]
                r1 = "".join([m.ref for m in mm])
                r2 = "".join([m.m_ref for m in mm])
                pos = mm[0].pos
                init_pos = pos
                alleles = [""] * cn
                vv = [v for i, v in to_proc if i == 0]
                star = set()
                for v in vv:
                    gt = [int(g) for g in re.split(r"\||/", v.samples[0]["GT"])]
                    if "*" in v.alt:
                        star_idx = v.alt.index("*") + 1
                        star = star.union(
                            set([i for i, g in enumerate(gt) if g == star_idx])
                        )
                    if v.pos > pos:
                        alleles = [
                            a + r1[pos - init_pos : v.pos - init_pos] for a in alleles
                        ]
                    for i, g in enumerate(gt):
                        if g == 0:
                            alleles[i] += v.ref
                        else:
                            alleles[i] += v.alt[g - 1]
                    pos = v.pos + len(v.ref)
                if pos != mm[-1].pos + len(mm[-1].ref):
                    alleles = [a + r1[pos - init_pos :] for a in alleles]
                for i in star:
                    alleles[i] = "*"
                alts = list(set([a for a in alleles if a != r2 and a != "*"]))
                if star:
                    alts.append("*")
                if not alts:
                    continue
                gt = ["0" if a == r2 else str(alts.index(a) + 1) for a in alleles]
                v = copy.deepcopy(vv[0])
                v.pos = mm[0].m_pos
                min_len = min([len(a) for a in [r2] + alts if a != "*"])
                for i in range(min_len - 1):
                    last = set([a[-1] for a in [r2] + alts if a != "*"])
                    if len(last) == 1:  # trim the last bit
                        r2 = r2[:-1]
                        alts = [a[:-1] if a != "*" else a for a in alts]
                    else:
                        break
                min_len = min([len(a) for a in [r2] + alts if a != "*"])
                for i in range(min_len - 1):
                    first = set([a[0] for a in [r2] + alts if a != "*"])
                    if len(first) == 1:  # trim the last bit
                        r2 = r2[1:]
                        alts = [a[1:] if a != "*" else a for a in alts]
                        v.pos += 1
                    else:
                        break
                v.ref = r2
                v.alt = alts
                sep = "/" if "/" in v.samples[0]["GT"] else "|"
                v.samples[0]["GT"] = sep.join(gt)
                v.line = None
                yield annotate_fn(v) if annotate_fn else v


def load_bam(in_bam: str, in_ref: Optional[str] = None) -> pysam.AlignmentFile:
    if in_bam.endswith(".cram"):
        if not in_ref:
            logger.error("No reference file provided for CRAM input")
            sys.exit(1)
        bamh = pysam.AlignmentFile(in_bam, "rc", reference_filename=in_ref)
    else:
        bamh = pysam.AlignmentFile(in_bam, "rb")
    return bamh
