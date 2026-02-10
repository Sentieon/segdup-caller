"""Statistical models for gene copy number analysis."""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Any
from dataclasses import dataclass
import numpy as np
import scipy.stats
from scipy.special import gammaln
import logging
from genecaller.logging import get_logger


class GenePriors(ABC):
    """Base class for gene copy number priors."""

    def __init__(self, config: Dict[str, Any], cn_regions: Dict[int, Dict[str, Any]]):
        self.config = config
        self.cn_regions = cn_regions
        self.logger = logging.getLogger(self.__class__.__name__)

    @abstractmethod
    def get_log_prior(self, cn: int, region_name: Optional[str] = None) -> float:
        """Get log prior probability for a copy number state."""
        pass


class GaussianPriors(GenePriors):
    """Gaussian-based copy number priors centered at diploid (CN=2)."""

    def __init__(self, config: Dict[str, Any], cn_regions: Dict[int, Dict[str, Any]]):
        super().__init__(config, cn_regions)
        self.cn_prior_std = config.get("cn_prior_std", 2.0)

        # Pre-calculate log priors for CN 0-8
        self._cached_log_priors = {}
        for cn in range(9):
            self._cached_log_priors[cn] = self._calculate_log_prior(cn)

        self.logger.debug(f"Using Gaussian CN priors with std={self.cn_prior_std}")

    def _calculate_log_prior(self, cn: int) -> float:
        if cn <= 0:
            upper = scipy.stats.norm.cdf(0.5, loc=2, scale=self.cn_prior_std)
            prob = upper
        else:
            upper = scipy.stats.norm.cdf(cn + 0.5, loc=2, scale=self.cn_prior_std)
            lower = scipy.stats.norm.cdf(cn - 0.5, loc=2, scale=self.cn_prior_std)
            prob = upper - lower
        return np.log(max(float(prob), 1e-10))

    def get_log_prior(self, cn: int, region_name: Optional[str] = None) -> float:
        if cn in self._cached_log_priors:
            return self._cached_log_priors[cn]
        return self._calculate_log_prior(cn)


class CategoricalPriors(GenePriors):
    """Categorical (discrete) copy number priors.

    Supports both gene-level and region-level configurations:
    - Gene-level: cn_priors: {0: 0.01, 1: 0.1, 2: 0.89}
    - Region-level: cn_priors: {CFH: {0: 0.001, 1: 0.01, 2: 0.98}, CFH31: {0: 0.05, 1: 0.15, 2: 0.79}}
    """

    def __init__(self, config: Dict[str, Any], cn_regions: Dict[int, Dict[str, Any]]):
        super().__init__(config, cn_regions)

        priors_config = config["cn_priors"]

        # Detect format based on value type
        first_value = next(iter(priors_config.values()))
        self.is_per_region = isinstance(first_value, dict)

        if self.is_per_region:
            self._init_region_level_priors(priors_config)
        else:
            self._init_gene_level_priors(priors_config)

    def _init_gene_level_priors(self, priors_config: Dict[int, float]) -> None:
        priors = dict(priors_config)

        total_prob = sum(priors.values())
        if not np.isclose(total_prob, 1.0, rtol=1e-6):
            self.logger.warning(
                f"CN priors sum to {total_prob:.6f}, not 1.0. Normalizing..."
            )
            priors = {cn: prob / total_prob for cn, prob in priors.items()}

        self._log_priors = {}
        for cn, prob in priors.items():
            if prob <= 0:
                self.logger.warning(
                    f"CN prior for CN={cn} is {prob} <= 0. Setting to very small value."
                )
                self._log_priors[cn] = -1e9
            else:
                self._log_priors[cn] = np.log(prob)

        self.logger.debug(f"Using gene-level categorical CN priors: {priors}")

    def _init_region_level_priors(
        self, priors_config: Dict[str, Dict[int, float]]
    ) -> None:
        # Build set of valid region names
        valid_region_names = set()
        for region_info in self.cn_regions.values():
            valid_region_names.add(region_info["name"])

        self._log_priors = {}

        for region_name, priors in priors_config.items():
            if region_name not in valid_region_names:
                self.logger.warning(
                    f"CN prior specified for '{region_name}' but no matching cn_region found. Skipping."
                )
                continue

            total_prob = sum(priors.values())
            if not np.isclose(total_prob, 1.0, rtol=1e-6):
                self.logger.warning(
                    f"CN priors for {region_name} sum to {total_prob:.6f}, not 1.0. Normalizing..."
                )
                priors = {cn: prob / total_prob for cn, prob in priors.items()}

            self._log_priors[region_name] = {}
            for cn, prob in priors.items():
                if prob <= 0:
                    self.logger.warning(
                        f"CN prior for {region_name} CN={cn} is {prob} <= 0. Setting to very small value."
                    )
                    self._log_priors[region_name][cn] = -1e9
                else:
                    self._log_priors[region_name][cn] = np.log(prob)

            self.logger.debug(
                f"Using categorical CN priors for {region_name}: {priors}"
            )

    def get_log_prior(self, cn: int, region_name: Optional[str] = None) -> float:
        if self.is_per_region:
            if region_name is not None and region_name in self._log_priors:
                if cn in self._log_priors[region_name]:
                    return self._log_priors[region_name][cn]
            self.logger.debug(
                f"No prior found for CN={cn}, region_name={region_name}. Using default small value."
            )
            return -1e9
        else:
            if cn in self._log_priors:
                return self._log_priors[cn]
            self.logger.debug(f"No prior found for CN={cn}. Using default small value.")
            return -1e9


class GammaPriors(GenePriors):
    """Gamma-like distribution priors for copy number."""

    def __init__(self, config: Dict[str, Any], cn_regions: Dict[int, Dict[str, Any]]):
        super().__init__(config, cn_regions)

        self.gamma_shape = config.get("gamma_shape", 4.0)
        self.gamma_scale = config.get("gamma_scale", 0.5)

        # Pre-calculate log priors for CN 0-8
        self._cached_log_priors = {}
        for cn in range(9):
            self._cached_log_priors[cn] = self._calculate_log_prior(cn)

        mean = self.gamma_shape * self.gamma_scale
        self.logger.debug(
            f"Using Gamma CN priors with shape={self.gamma_shape}, "
            f"scale={self.gamma_scale}, mean={mean:.2f}"
        )

    def _calculate_log_prior(self, cn: int) -> float:
        if cn <= 0:
            prob = scipy.stats.gamma.cdf(
                0.5, a=self.gamma_shape, scale=self.gamma_scale
            )
        else:
            upper = scipy.stats.gamma.cdf(
                cn + 0.5, a=self.gamma_shape, scale=self.gamma_scale
            )
            lower = scipy.stats.gamma.cdf(
                cn - 0.5, a=self.gamma_shape, scale=self.gamma_scale
            )
            prob = upper - lower
        return np.log(max(float(prob), 1e-10))

    def get_log_prior(self, cn: int, region_name: Optional[str] = None) -> float:
        if cn in self._cached_log_priors:
            return self._cached_log_priors[cn]
        return self._calculate_log_prior(cn)


def create_gene_priors(
    config: Dict[str, Any], cn_regions: Dict[int, Dict[str, Any]]
) -> GenePriors:
    """Factory method to create appropriate GenePriors instance."""
    prior_type = config.get("prior_type", "gaussian")

    if "cn_priors" in config:
        prior_type = "categorical"

    if prior_type == "categorical":
        return CategoricalPriors(config, cn_regions)
    elif prior_type == "gamma":
        return GammaPriors(config, cn_regions)
    else:
        return GaussianPriors(config, cn_regions)


@dataclass
class Signal:
    """Signal value at a genomic position with quality weight"""

    position: int
    id: str
    ads: list[int]
    ratio: float
    logodds: float
    zscore: float
    weight: float
    segment_alt_allele_count: Optional[int] = None

    def is_conversion_site(self):
        return self.id and self.id != "." and len(self.id.split("_")) != 3


def calc_conversion_signals(records) -> list[Signal]:
    """Calculate conversion signals from variants using multiple methods."""
    logger = get_logger(__name__)
    logger.debug(f"Calculating conversion signals from {len(records)} records")

    alt_counts = []
    for rec in records:  # rec is a tuple (id, pos, ads)
        ads = rec[2]
        alt_counts.append(ads[1])

    if not alt_counts:
        logger.warning("No valid alt counts found in variants")
        return []

    mean = np.mean(alt_counts)
    std = np.std(alt_counts) + 1e-6
    logger.debug(f"  Alt count statistics: mean={mean:.2f}, std={std:.2f}")

    signals = []
    for rec in records:
        ads = rec[2]
        if len(ads) != 2:
            continue
        ratio = ads[1] / (ads[0] + ads[1] + 1e-6)
        log_odds = np.log((ads[1] + 1) / (ads[0] + 1))
        z_score = (ads[1] - mean) / std
        signals.append(
            Signal(
                position=rec[1],
                id=rec[0],
                ads=ads,
                ratio=ratio,
                logodds=log_odds,
                zscore=z_score,
                weight=1.0,
            )
        )

    logger.debug(f"  Generated {len(signals)} signals")
    if signals:
        logger.debug(
            f"  Ratio range: [{min(s.ratio for s in signals):.3f}, {max(s.ratio for s in signals):.3f}]"
        )
        logger.debug(
            f"  Z-score range: [{min(s.zscore for s in signals):.3f}, {max(s.zscore for s in signals):.3f}]"
        )

    return signals


@dataclass
class GeneConversionSegment:
    start: int
    end: int
    allele_counts: list[int]
    mean_signal: float
    signals: list[Signal]
    conversion_allele_count: Optional[int] = None
    is_fusion_candidate: bool = False

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def num_sites(self) -> int:
        return len(self.signals)


class GeneSegmenter(ABC):
    segments: list[GeneConversionSegment]

    def __init__(self, signals):
        self.signals = signals
        self.breakpoints = []
        self.segments = []
        self.total_cn = -1
        self.logger = get_logger(self.__class__.__name__)

    def segment(self):
        # Implement segmentation logic here
        pass

    def get_breakpoints(self):
        return self.breakpoints

    def get_segments(self):
        if len(self.breakpoints) < 2:
            return []

        segments = []
        for i in range(len(self.breakpoints) - 1):
            start = self.breakpoints[i]
            end = self.breakpoints[i + 1]

            seg_signals = [s for s in self.signals if start <= s.position < end]
            if not seg_signals:
                continue

            alt_allele_count = (
                seg_signals[0].segment_alt_allele_count
                if seg_signals[0].segment_alt_allele_count is not None
                else -1
            )
            mean_ratio = np.mean([s.ratio for s in seg_signals])

            segment = GeneConversionSegment(
                start=start,
                end=end,
                allele_counts=[self.total_cn - alt_allele_count, alt_allele_count],
                mean_signal=mean_ratio,  # type: ignore
                signals=seg_signals,
                conversion_allele_count=None,
            )
            segments.append(segment)

        return segments

    def _merge_adjacent_identical_segments(self):
        """
        Merge adjacent segments with identical allele counts.

        This prevents post-processing steps (like _merge_small_segments in HMM)
        from creating adjacent segments with identical states after recalculation.

        This method should be called by all segmentation algorithms after initial
        segmentation to ensure consistent behavior across CBS, PELT, and HMM.
        """
        if not self.segments or len(self.segments) <= 1:
            return

        merged = True
        while merged:
            merged = False
            i = 0
            while i < len(self.segments) - 1:
                current = self.segments[i]
                next_seg = self.segments[i + 1]

                # Check if adjacent segments have identical allele counts
                # Use allele_counts[1] (alt allele count) for comparison
                if current.allele_counts[1] == next_seg.allele_counts[1]:
                    self.logger.debug(
                        f"  Merging adjacent identical segments: "
                        f"[{current.start}-{current.end}, AC={current.allele_counts[1]}] + "
                        f"[{next_seg.start}-{next_seg.end}, AC={next_seg.allele_counts[1]}] "
                        f"-> [{current.start}-{next_seg.end}]"
                    )

                    # Merge the two segments
                    merged_signals = current.signals + next_seg.signals
                    merged_seg = GeneConversionSegment(
                        start=current.start,
                        end=next_seg.end,
                        allele_counts=current.allele_counts,  # Same as next_seg.allele_counts
                        mean_signal=float(np.mean([s.ratio for s in merged_signals])),
                        signals=merged_signals,
                        conversion_allele_count=None,
                    )

                    # Replace current segment with merged segment and remove next
                    self.segments[i] = merged_seg
                    self.segments.pop(i + 1)
                    merged = True
                    # Don't increment i - check the new merged segment against its new neighbor
                else:
                    i += 1

    def _assign_cn_state(self, signals):
        best_state = 0
        best_ll = -np.inf

        alt_counts = np.array([s.ads[1] for s in signals])
        total_counts = np.array([s.ads[0] + s.ads[1] for s in signals])

        valid_mask = total_counts > 0
        if not np.any(valid_mask):
            return 0

        alt_counts = alt_counts[valid_mask]
        total_counts = total_counts[valid_mask]

        for state in range(self.total_cn + 1):
            # Use 0.01 epsilon to model ~1% sequencing/mapping error rate
            # This prevents single outlier signals from dominating the likelihood
            p = np.clip(state / self.total_cn, 0.01, 0.99)
            all_lls = scipy.stats.binom.logpmf(alt_counts, total_counts, p)
            ll = np.sum(all_lls)

            if ll > best_ll:
                best_ll = ll
                best_state = state

        return best_state

    def _extract_breakpoints(self):
        if len(self.breakpoints) < 2:
            self.segments = []
            return

        segments = []
        for i in range(len(self.breakpoints) - 1):
            start_pos = self.breakpoints[i]
            end_pos = self.breakpoints[i + 1]

            seg_signals = [s for s in self.signals if start_pos <= s.position < end_pos]
            if not seg_signals:
                continue

            best_state = self._assign_cn_state(seg_signals)
            mean_ratio = np.mean([s.ratio for s in seg_signals])

            segments.append(
                {
                    "start": start_pos,
                    "end": end_pos,
                    "signals": seg_signals,
                    "cn_state": best_state,
                    "mean_ratio": mean_ratio,
                }
            )

        while True:
            merge_candidates = []
            for i in range(len(segments) - 1):
                if segments[i]["cn_state"] == segments[i + 1]["cn_state"]:
                    ratio_diff = abs(
                        segments[i]["mean_ratio"] - segments[i + 1]["mean_ratio"]
                    )
                    merge_candidates.append((i, ratio_diff))

            if not merge_candidates:
                break

            merge_idx, _ = min(merge_candidates, key=lambda x: x[1])

            merged_signals = (
                segments[merge_idx]["signals"] + segments[merge_idx + 1]["signals"]
            )
            merged_state = self._assign_cn_state(merged_signals)
            merged_ratio = np.mean([s.ratio for s in merged_signals])

            new_segment = {
                "start": segments[merge_idx]["start"],
                "end": segments[merge_idx + 1]["end"],
                "signals": merged_signals,
                "cn_state": merged_state,
                "mean_ratio": merged_ratio,
            }

            segments = segments[:merge_idx] + [new_segment] + segments[merge_idx + 2 :]

        self.breakpoints = [segments[0]["start"]]
        for seg in segments:
            self.breakpoints.append(seg["end"])
        self.breakpoints = sorted(list(set(self.breakpoints)))

        for seg in segments:
            for signal in seg["signals"]:
                signal.segment_alt_allele_count = seg["cn_state"]

        # Convert to GeneConversionSegment objects and store
        self.segments = []
        for seg in segments:
            alt_allele_count = seg["cn_state"]
            segment = GeneConversionSegment(
                start=seg["start"],
                end=seg["end"],
                allele_counts=[self.total_cn - alt_allele_count, alt_allele_count],
                mean_signal=seg["mean_ratio"],
                signals=seg["signals"],
                conversion_allele_count=None,
            )
            self.segments.append(segment)


class CBSGeneSegmenter(GeneSegmenter):
    def __init__(self, signals: list[Signal], params={}):
        super().__init__(signals)
        self.alpha = params.get("alpha", 0.01)
        self.min_width = params.get("min_width", 2)
        self.value_type = params.get("value_type", "ratio")
        self.total_cn = params.get("total_cn", -1)

    def segment(self):
        if not self.signals:
            self.logger.warning("CBS: No signals to segment")
            self.breakpoints = []
            self.segments = []
            return

        self.logger.info(
            f"CBS segmentation: {len(self.signals)} signals, alpha={self.alpha}, min_width={self.min_width}"
        )

        self.positions = [s.position for s in self.signals]
        values = np.array(
            [s.ratio if self.value_type == "ratio" else s.logodds for s in self.signals]
        )
        weights = np.array([s.weight for s in self.signals])

        n = len(values)
        self.cum_W = np.zeros(n + 1)
        self.cum_S = np.zeros(n + 1)
        self.cum_S2 = np.zeros(n + 1)

        for i in range(n):
            self.cum_W[i + 1] = self.cum_W[i] + weights[i]
            self.cum_S[i + 1] = self.cum_S[i] + weights[i] * values[i]
            self.cum_S2[i + 1] = self.cum_S2[i] + weights[i] * values[i] ** 2

        self.breakpoints = [self.positions[0]]
        self._recursive_segment(0, n)
        self.breakpoints.append(self.positions[-1] + 1)
        self.breakpoints = sorted(set(self.breakpoints))

        self.logger.debug(
            f"  Found {len(self.breakpoints)} breakpoints: {self.breakpoints[:10]}..."
            if len(self.breakpoints) > 10
            else f"  Found {len(self.breakpoints)} breakpoints: {self.breakpoints}"
        )

        if self.total_cn > 0:
            self.logger.debug(f"  Extracting breakpoints with total_cn={self.total_cn}")
            self._extract_breakpoints()
        else:
            self.segments = self.get_segments()

        # Merge adjacent segments with identical allele counts
        self._merge_adjacent_identical_segments()

        self.logger.info(f"CBS result: {len(self.segments)} segments")
        for i, seg in enumerate(self.segments):
            self.logger.debug(
                f"    Segment {i + 1}: {seg.start}-{seg.end}, allele_counts={seg.allele_counts}, mean_signal={seg.mean_signal:.3f}, {seg.num_sites} sites"
            )

    def _recursive_segment(self, start, end):
        n = end - start
        if n < 2 * self.min_width:
            return

        best_split = None
        best_tstat = 0

        for i in range(start + self.min_width, end - self.min_width):
            W_left = self.cum_W[i] - self.cum_W[start]
            S_left = self.cum_S[i] - self.cum_S[start]
            S2_left = self.cum_S2[i] - self.cum_S2[start]

            W_right = self.cum_W[end] - self.cum_W[i]
            S_right = self.cum_S[end] - self.cum_S[i]
            S2_right = self.cum_S2[end] - self.cum_S2[i]

            if W_left == 0 or W_right == 0:
                continue

            mean_left = S_left / W_left
            mean_right = S_right / W_right

            var_left = max(0, (S2_left / W_left) - mean_left**2)
            var_right = max(0, (S2_right / W_right) - mean_right**2)

            pooled_var = (W_left * var_left + W_right * var_right) / (W_left + W_right)

            if pooled_var == 0:
                continue

            t_stat = abs(mean_left - mean_right) / np.sqrt(
                pooled_var * (1 / W_left + 1 / W_right)
            )

            if t_stat > best_tstat:
                best_tstat = t_stat
                best_split = i

        if best_split is not None:
            df = n - 2
            p_value = 2 * (1 - scipy.stats.t.cdf(best_tstat, df))

            if p_value < self.alpha:
                self.breakpoints.append(self.positions[best_split])
                self._recursive_segment(start, best_split)
                self._recursive_segment(best_split, end)


class HMMSEGGeneSegmenter(GeneSegmenter):
    def __init__(self, signals, params={}):
        super().__init__(signals)
        self.total_cn = params.get("total_cn", -1)
        if self.total_cn <= 0:
            raise ValueError("HMMSEG requires 'total_cn' parameter > 0")
        self.n_states = self.total_cn + 1
        self.min_num_signals = params.get("min_num_signals", 1)

        # Get transition rates - support both legacy single rate and new dual rates
        legacy_transition_rate = float(params.get("transition_rate", 0.01))
        self.transition_rate_background = float(params.get("transition_rate_background", legacy_transition_rate))
        self.transition_rate_hotspot = float(params.get("transition_rate_hotspot", legacy_transition_rate))

        # Store for logging
        self.transition_rate = legacy_transition_rate

        # Parse recombination regions if provided
        self.recombination_regions = params.get("recombination_regions", None)
        self.hotspot_intervals = None

        if self.recombination_regions:
            # Parse recombination regions into IntervalList for efficient lookup
            from .util import IntervalList

            # Convert list of region strings to comma-separated format
            if isinstance(self.recombination_regions, list):
                region_str = ",".join(self.recombination_regions)
            else:
                region_str = self.recombination_regions
            self.hotspot_intervals = IntervalList(region=region_str)
            self.logger.debug(
                f"  Recombination hotspots configured: {self.recombination_regions}"
            )
            self.logger.debug(
                f"  Transition rates: background={self.transition_rate_background:.2e}, hotspot={self.transition_rate_hotspot:.2e}"
            )

        # Precompute log transition probability matrices
        self.log_trans_matrix_cold = self._build_trans_matrix(
            self.transition_rate_background
        )
        self.log_trans_matrix_hot = self._build_trans_matrix(
            self.transition_rate_hotspot
        )

        # For backward compatibility, keep single matrix reference
        self.log_trans_matrix = self.log_trans_matrix_cold

    def _build_trans_matrix(self, transition_rate):
        """Build a transition matrix for given transition rate."""
        log_stay = np.log(1 - transition_rate + 1e-10)
        log_jump = np.log(transition_rate / (self.n_states - 1) + 1e-10)

        # This matrix holds log_prob(from_state, to_state)
        # Rows (axis 0) = from_state, Columns (axis 1) = to_state
        matrix = np.full((self.n_states, self.n_states), log_jump)
        np.fill_diagonal(matrix, log_stay)
        return matrix

    def _is_hotspot(self, position):
        """Check if a position falls within a recombination hotspot."""
        if self.hotspot_intervals is None:
            return False
        # IntervalList.__contains__ with tuple checks if position is in any interval
        # Need to extract chromosome from signals - assuming all signals are from same chromosome
        if not self.signals:
            return False
        # IntervalList stores regions by chromosome, we need to know which chromosome
        # For now, check all chromosomes in the hotspot_intervals
        for chrom in self.hotspot_intervals.regions:
            if (chrom, position) in self.hotspot_intervals:
                return True
        return False

    def segment(self):
        if not self.signals:
            self.logger.warning("HMMSEG: No signals to segment")
            self.breakpoints = []
            self.segments = []
            return

        # Validate recombination regions overlap with signal range
        if self.hotspot_intervals is not None:
            signal_start = min(s.position for s in self.signals)
            signal_end = max(s.position for s in self.signals)

            # Check if any hotspot overlaps with signal range
            has_overlap = False
            for chrom in self.hotspot_intervals.regions:
                intervals = self.hotspot_intervals.get(chrom, signal_start, signal_end)
                if intervals:
                    has_overlap = True
                    break

            if not has_overlap:
                self.logger.warning(
                    f"Recombination regions configured but none overlap with signal range "
                    f"[{signal_start}-{signal_end}]. Check coordinate system (hg19 vs hg38)."
                )

        self.logger.info(
            f"HMMSEG segmentation: {len(self.signals)} signals, {self.n_states} states, transition_rate={self.transition_rate}, min_num_signals={self.min_num_signals}"
        )

        state_sequence = self._viterbi(self.signals)

        # Log state distribution
        state_counts = np.bincount(state_sequence, minlength=self.n_states)
        self.logger.debug(f"  State distribution: {dict(enumerate(state_counts))}")

        for i, state in enumerate(state_sequence):
            self.signals[i].segment_alt_allele_count = state

        positions = [s.position for s in self.signals]
        self.breakpoints = self._extract_breakpoints(positions, state_sequence)
        self.segments = self.get_segments()

        self.logger.info(
            f"HMMSEG initial result: {len(self.segments)} segments, {len(self.breakpoints)} breakpoints"
        )

        # Merge adjacent segments with identical allele counts
        self._merge_adjacent_identical_segments()
        self.logger.info(
            f"HMMSEG after merging adjacent: {len(self.segments)} segments"
        )

        # Apply post-processing to merge small segments
        if self.min_num_signals > 1:
            self._merge_small_segments()
            self.logger.info(
                f"HMMSEG after merging small: {len(self.segments)} segments"
            )

            # Merge adjacent identical segments again after _merge_small_segments
            # because the recalculation during small segment merging can create
            # new adjacent segments with identical allele counts
            self._merge_adjacent_identical_segments()
            self.logger.info(f"HMMSEG after final merge: {len(self.segments)} segments")

        for i, seg in enumerate(self.segments):
            self.logger.debug(
                f"    Segment {i + 1}: {seg.start}-{seg.end}, allele_counts={seg.allele_counts}, mean_signal={seg.mean_signal:.3f}, {seg.num_sites} sites"
            )

    def _viterbi(self, signals):
        n_obs = len(signals)
        n_states = self.n_states

        log_V = np.zeros((n_obs, n_states))
        path = np.zeros((n_obs, n_states), dtype=int)

        # initialization
        first_emission_log_probs = np.array(
            [self._emission_log_prob(signals[0], s) for s in range(n_states)]
        )
        log_V[0, :] = -np.log(n_states) + first_emission_log_probs

        # Track hotspot transitions for diagnostic logging
        hotspot_transitions = 0
        background_transitions = 0

        # Recursion
        for t in range(1, n_obs):
            # Calculate all emission probs for this step
            emission_ll = np.array(
                [self._emission_log_prob(signals[t], s) for s in range(n_states)]
            )

            # Select transition matrix based on current position
            if self._is_hotspot(signals[t].position):
                current_trans_matrix = self.log_trans_matrix_hot
                hotspot_transitions += 1
            else:
                current_trans_matrix = self.log_trans_matrix_cold
                background_transitions += 1

            scores_matrix = log_V[t - 1, :, np.newaxis] + current_trans_matrix

            # Now, find the best path *into* each current state 's'
            # (i.e., find max along the 'prev_s' axis 0)
            path[t, :] = np.argmax(scores_matrix, axis=0)
            log_V[t, :] = np.max(scores_matrix, axis=0) + emission_ll

        # Backtracking
        states = np.zeros(n_obs, dtype=int)
        states[-1] = np.argmax(log_V[-1, :])

        for t in range(n_obs - 2, -1, -1):
            states[t] = path[t + 1, states[t + 1]]

        # Log diagnostic information
        if self.hotspot_intervals is not None:
            self.logger.debug(
                f"  Viterbi: {hotspot_transitions} hotspot transitions, {background_transitions} background transitions"
            )

        return states

    def _emission_log_prob(self, signal: Signal, state: int):
        """
        Calculates the Binomial log-likelihood of observing the
        read counts given the hidden state (integer allele number).
        """
        alt_count = signal.ads[1]
        total_count = signal.ads[0] + signal.ads[1]

        if total_count == 0:
            return 0.0  # No evidence, neutral log-prob

        # The expected ratio p given the state k
        expected_ratio = state / self.total_cn

        # Use 0.01 epsilon to model ~1% sequencing/mapping error rate
        # This prevents single outlier signals from dominating the likelihood
        p = np.clip(expected_ratio, 0.01, 0.99)

        return scipy.stats.binom.logpmf(alt_count, total_count, p)

    def _extract_breakpoints(self, positions, states):
        breakpoints = [positions[0]]
        for i in range(1, len(states)):
            if states[i] != states[i - 1]:
                breakpoints.append(positions[i])
        breakpoints.append(positions[-1] + 1)
        return breakpoints

    def _merge_small_segments(self):
        """Merge segments with fewer than min_num_signals with their neighbors."""
        if not self.segments or len(self.segments) <= 1:
            return

        merged = True
        while merged:
            merged = False
            i = 0
            while i < len(self.segments):
                seg = self.segments[i]
                if seg.num_sites < self.min_num_signals:
                    self.logger.debug(
                        f"  Merging small segment at {seg.start}-{seg.end} ({seg.num_sites} sites)"
                    )

                    # Decide which neighbor to merge with
                    # Prefer merging with neighbor that has most similar allele count
                    merge_with_prev = False
                    if i == 0:
                        # First segment, merge with next
                        merge_with_prev = False
                    elif i == len(self.segments) - 1:
                        # Last segment, merge with previous
                        merge_with_prev = True
                    else:
                        # Middle segment, choose neighbor with most similar allele count
                        prev_seg = self.segments[i - 1]
                        next_seg = self.segments[i + 1]
                        diff_prev = abs(
                            seg.allele_counts[1] - prev_seg.allele_counts[1]
                        )
                        diff_next = abs(
                            seg.allele_counts[1] - next_seg.allele_counts[1]
                        )
                        merge_with_prev = diff_prev <= diff_next

                    if merge_with_prev and i > 0:
                        # Merge with previous
                        prev_seg = self.segments[i - 1]
                        merged_signals = prev_seg.signals + seg.signals
                        merged_allele_count = self._assign_cn_state(merged_signals)
                        merged_seg = GeneConversionSegment(
                            start=prev_seg.start,
                            end=seg.end,
                            allele_counts=[
                                self.total_cn - merged_allele_count,
                                merged_allele_count,
                            ],
                            mean_signal=float(
                                np.mean([s.ratio for s in merged_signals])
                            ),
                            signals=merged_signals,
                            conversion_allele_count=None,
                        )
                        self.segments[i - 1] = merged_seg
                        self.segments.pop(i)
                        merged = True
                        self.logger.debug(
                            f"    Merged with previous: {merged_seg.start}-{merged_seg.end} ({merged_seg.num_sites} sites)"
                        )
                        # Don't increment i - check the segment that shifted down to position i
                    elif not merge_with_prev and i < len(self.segments) - 1:
                        # Merge with next
                        next_seg = self.segments[i + 1]
                        merged_signals = seg.signals + next_seg.signals
                        merged_allele_count = self._assign_cn_state(merged_signals)
                        merged_seg = GeneConversionSegment(
                            start=seg.start,
                            end=next_seg.end,
                            allele_counts=[
                                self.total_cn - merged_allele_count,
                                merged_allele_count,
                            ],
                            mean_signal=float(
                                np.mean([s.ratio for s in merged_signals])
                            ),
                            signals=merged_signals,
                            conversion_allele_count=None,
                        )
                        self.segments[i] = merged_seg
                        self.segments.pop(i + 1)
                        merged = True
                        self.logger.debug(
                            f"    Merged with next: {merged_seg.start}-{merged_seg.end} ({merged_seg.num_sites} sites)"
                        )
                        # Don't increment i - check the merged segment at position i again
                    else:
                        # Can't merge (shouldn't happen), move to next
                        i += 1
                else:
                    # Segment is large enough, move to next
                    i += 1

        # Update breakpoints after merging
        if self.segments:
            self.breakpoints = [self.segments[0].start]
            for seg in self.segments:
                self.breakpoints.append(seg.end)
            self.breakpoints = sorted(list(set(self.breakpoints)))


class PELTGeneSegmenter(GeneSegmenter):
    def __init__(self, signals, params={}):
        super().__init__(signals)
        self.penalty = params.get("penalty", 2 * np.log(max(1, len(signals))))
        self.min_width = params.get("min_width", 2)
        self.total_cn = params.get("total_cn", -1)

    def segment(self):
        if not self.signals:
            self.logger.warning("PELT: No signals to segment")
            self.breakpoints = []
            self.segments = []
            return

        n = len(self.signals)
        self.logger.info(
            f"PELT segmentation: {n} signals, penalty={self.penalty:.2f}, min_width={self.min_width}"
        )

        # 1. Pre-calculate cumulative sums
        cum_Alt = np.zeros(n + 1)
        cum_Total = np.zeros(n + 1)

        # This is the new, third cumulative sum for the log-coefficient
        cum_LogCoeff = np.zeros(n + 1)

        for i in range(n):
            alt_c = self.signals[i].ads[1]
            total_c = self.signals[i].ads[0] + self.signals[i].ads[1]

            cum_Alt[i + 1] = cum_Alt[i] + alt_c
            cum_Total[i + 1] = cum_Total[i] + total_c

            if total_c > 0:
                # log(n!) = gammaln(n+1)
                # log(n C k) = log(n!) - log(k!) - log((n-k)!)
                log_coeff = (
                    gammaln(total_c + 1)
                    - gammaln(alt_c + 1)
                    - gammaln(total_c - alt_c + 1)
                )
                cum_LogCoeff[i + 1] = cum_LogCoeff[i] + log_coeff

        # 2. Initialize PELT data structures
        F = np.full(n + 1, np.inf)
        F[0] = -self.penalty
        cp = np.zeros(n + 1, dtype=int)
        candidates = [0]

        # 3. Run the PELT algorithm
        for t in range(1, n + 1):
            min_cost = np.inf
            best_tau = 0

            for tau in candidates:
                if t - tau < self.min_width:
                    continue

                # This cost calculation is now O(1)
                cost = (
                    F[tau]
                    + self._segment_cost(tau, t, cum_Alt, cum_Total, cum_LogCoeff)
                    + self.penalty
                )
                if cost < min_cost:
                    min_cost = cost
                    best_tau = tau

            F[t] = min_cost
            cp[t] = best_tau

            # --- Pruning Step ---
            pruned_candidates = []
            for s in candidates:
                # This is also O(1)
                cost_s = F[s] + self._segment_cost(
                    s, t, cum_Alt, cum_Total, cum_LogCoeff
                )
                if cost_s <= F[t]:
                    pruned_candidates.append(s)

            pruned_candidates.append(t)
            candidates = pruned_candidates

        self.logger.debug(f"  PELT optimization complete, final cost={F[n]:.2f}")

        # 4. Backtrack to find changepoints
        changepoints = []
        current = n
        while current > 0:
            changepoints.append(current)
            current = cp[current]
        changepoints.reverse()

        self.logger.debug(f"  Found {len(changepoints)} changepoints")

        # 5. Convert indices back to genomic positions
        positions = [s.position for s in self.signals]
        self.breakpoints = [positions[0]]
        for i in changepoints[:-1]:
            if 0 < i < n:
                self.breakpoints.append(positions[i])
        self.breakpoints.append(positions[-1] + 1)
        self.breakpoints = sorted(list(set(self.breakpoints)))

        self.logger.debug(f"  Converted to {len(self.breakpoints)} genomic breakpoints")

        if self.total_cn > 0:
            self.logger.debug(f"  Extracting breakpoints with total_cn={self.total_cn}")
            self._extract_breakpoints()
        else:
            self.segments = self.get_segments()

        # Merge adjacent segments with identical allele counts
        self._merge_adjacent_identical_segments()

        self.logger.info(f"PELT result: {len(self.segments)} segments")
        for i, seg in enumerate(self.segments):
            self.logger.debug(
                f"    Segment {i + 1}: {seg.start}-{seg.end}, allele_counts={seg.allele_counts}, mean_signal={seg.mean_signal:.3f}, {seg.num_sites} sites"
            )

    def _segment_cost(self, start: int, end: int, cum_Alt, cum_Total, cum_LogCoeff):
        """
        Calculates the O(1) Negative Log-Likelihood cost for a segment.
        Cost = -max(log(L))
        """

        # Get all components in O(1)
        sum_alt = cum_Alt[end] - cum_Alt[start]
        sum_total = cum_Total[end] - cum_Total[start]
        sum_log_coeff = cum_LogCoeff[end] - cum_LogCoeff[start]

        if sum_total == 0:
            return 0.0

        sum_ref = sum_total - sum_alt

        # p_hat = sum_alt / sum_total
        # We need log(p_hat) and log(1-p_hat), with checks for p=0 or p=1

        # Use logaddexp pattern to avoid -inf
        log_p_hat = np.log(sum_alt + 1e-10) - np.log(sum_total)
        log_1_minus_p_hat = np.log(sum_ref + 1e-10) - np.log(sum_total)

        # This is the maximized log-likelihood for the segment
        max_log_likelihood = (
            sum_log_coeff + sum_alt * log_p_hat + sum_ref * log_1_minus_p_hat
        )

        # Cost is the negative log-likelihood
        return -max_log_likelihood
