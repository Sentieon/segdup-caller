"""Statistical models for gene copy number analysis."""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Any
import numpy as np
import scipy.stats
import logging


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
