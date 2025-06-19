from .util import IntervalList, GeneMapping, get_data_file
from .logging import get_logger
import vcflib
import copy


class Gene:
    def __init__(self, cfg):
        self.gene_regions = cfg["gene_regions"]
        if "extract_regions" in cfg:
            self.liftover_regions = cfg["extract_regions"]
        else:
            self.liftover_regions = [
                ",".join([g for j, g in enumerate(self.gene_regions) if j != i])
                for i in range(len(self.gene_regions))
            ]
        if "realign_regions" in cfg:
            self.realign_regions = cfg["realign_regions"]
        else:
            self.realign_regions = [
                IntervalList(region=g).pad(cfg.get("padding", 200)).to_region_str()
                for g in self.gene_regions
            ]
        self.diff_vcf = [get_data_file(p) for p in cfg["diff_vcf"]]
        self.mapping = GeneMapping(get_data_file(cfg["map"]))
        self.gene_names = cfg["gene_names"]
        self.known_longdel = (
            cfg["longdel"].split(",") if "longdel" in cfg else []
        )  # list of intervals in gene1
        self.set_longdels([1] * len(self.known_longdel))
        self.chr = cfg["gene_regions"][0].split(",")[0].split(":")[0]
        self.logger = get_logger(self.__class__.__name__, "INFO")

    def set_longdels(self, dels):
        self.longdels = [self.known_longdel[i] for i, d in enumerate(dels) if d]
        self.longdel_region = (
            ",".join(self.longdels),
            ",".join([self.mapping.interval12(d) for d in self.longdels]),
        )
        self.nodel_region = [
            IntervalList(region=g).subtract(self.longdel_region[i]).to_region_str()
            for i, g in enumerate(self.gene_regions)
        ]
        self.liftover_longdel = (
            self.get_liftover_longdel(0),
            self.get_liftover_longdel(1),
        )
        self.liftover_nodel = [
            IntervalList(region=g).subtract(self.liftover_longdel[i]).to_region_str()
            for i, g in enumerate(self.gene_regions)
        ]

    @staticmethod
    def region2bed(region, fname):
        with open(fname, "w") as f:
            for r in region.split(","):
                chr, parts = r.split(":")
                start, end = [int(p) for p in parts]
                f.write(f"{chr}:{start - 1}:{end}\n")

    def get_liftover_longdel(self, which):
        cur_rgn = self.longdel_region[which]
        vcf = vcflib.VCF(self.diff_vcf[which])
        long_dels = [
            v for v in vcf if len(v.ref) - min([len(a) for a in v.alt]) > 100
        ]
        if not long_dels:
           return cur_rgn
        itv = IntervalList(region=cur_rgn)
        if cur_rgn:
            chrs = list(itv.regions.keys())
        else:
            chrs = [long_dels[0].chrom]
        for chr in chrs:
            cur_regions = []
            region = itv.regions.get(chr, [])
            i = 0
            while i < len(region):
                cur_regions.append([region[i], region[i + 1]])
                i += 2
            for v in long_dels:
                if v.chrom == chr:
                    cur_regions.append([v.pos, v.pos + len(v.ref)])
            itv.regions[chr] = itv.merge_regions(cur_regions)
        return itv.to_region_str()

    def init_result_data(self, read_data):
        data = []
        for i in range(len(self.gene_names)):
            data.append(
                {
                    "orig": {"nodel": {"region": self.nodel_region[i]}},
                    "liftover": {"nodel": {"region": self.liftover_nodel[i]}},
                }
            )
            for name, d in read_data.items():
                if "liftover_valid" in d:
                    data[-1]["liftover"]["nodel"][name] = {
                        "region": IntervalList(region=d["liftover_valid"][i])
                        .intersect(self.liftover_nodel[i])
                        .to_region_str()
                    }
                else:
                    data[-1]["liftover"]["nodel"][name] = {
                        "region": self.liftover_nodel[i]
                    }
            del_regions = dict(
                [
                    (del_region, str(j))
                    for j, del_region in enumerate(self.longdel_region[i].split(","))
                    if del_region
                ]
            )
            for del_region, j in del_regions.items():
                data[-1]["orig"][j] = {"region": del_region}
                for name, d in read_data.items():
                    data[-1]["orig"][j][name] = {}
            for name, d in read_data.items():
                data[-1]["orig"]["nodel"][name] = {}
            for del_region in self.liftover_longdel[i].split(","):
                if not del_region:
                    continue
                j = del_regions.get(del_region, del_region)
                data[-1]["liftover"][j] = {"region": del_region}
                for name, d in read_data.items():
                    if "liftover_valid" in d:
                        data[-1]["liftover"][j][name] = {
                            "region": IntervalList(region=d["liftover_valid"][i])
                            .intersect(del_region)
                            .to_region_str()
                        }
                    else:
                        data[-1]["liftover"][j][name] = {"region": del_region}
        return data

    def prepare_output(self, data):
        result_data = {}
        for i, d in enumerate(data):
            dd = copy.deepcopy(d["orig"]["all"]) if "all" in d["orig"] else {}
            dd["CN"] = {"overall": d["orig"]["nodel"]["cn"]}
            for r, _d in d["orig"].items():
                if r not in ("all", "nodel"):
                    dd["CN"][_d["region"]] = _d["cn"]
            result_data[self.gene_names[i]] = dd
        return result_data


class SMN1(Gene):
    def prepare_output(self, data):
        return super().prepare_output(data)
