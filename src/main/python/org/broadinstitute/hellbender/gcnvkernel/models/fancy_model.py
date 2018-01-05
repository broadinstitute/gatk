import pymc3 as pm

import logging
from pymc3 import Model
from typing import Optional, Set, Dict

_logger = logging.getLogger(__name__)


class GeneralizedContinuousModel(Model):
    """An extension of PyMC3 `Model` class with the added functionality of labeling RVs
    as either global or sample-specific (for the purpose of I/O, and structured optimization).
    """
    def __init__(self):
        self.approx: Optional[pm.MeanField] = None
        self.global_var_registry: Set[str] = set()
        self.sample_specific_var_registry: Dict[str, int] = dict()
        super().__init__()

    @staticmethod
    def _get_var_name(var) -> str:
        assert hasattr(var, 'name')
        name = var.name
        if hasattr(var, 'transformed'):
            name = var.transformed.name
        return name

    def _assert_var_is_unannotated(self, var):
        name = self._get_var_name(var)
        assert name not in self.sample_specific_var_registry, \
            "Variable ({0}) is already registered as sample-specific".format(name)
        assert name not in self.global_var_registry, \
            "Variable ({0}) is already registered as global".format(name)

    def register_as_global(self, var):
        """Register a variable as global.

        Args:
            var: a PyMC3 free variable

        Returns:
            None
        """
        self._assert_var_is_unannotated(var)
        name = self._get_var_name(var)
        self.global_var_registry.add(name)

    def register_as_sample_specific(self, var, sample_axis: int):
        """Register a variable as sample-specific.

        Args:
            var: a PyMC3 free variable
            sample_axis: axis corresponding to sample index (it is used for slicing `var` to obtain single-sample
                parameters)

        Returns:
            None
        """
        self._assert_var_is_unannotated(var)
        name = self._get_var_name(var)
        self.sample_specific_var_registry[name] = sample_axis

    def verify_var_registry(self):
        """Verifies that all variables are registered as either as global or sample-specific."""
        _logger.info("Global model variables: " + str(self.global_var_registry))
        _logger.info("Sample-specific model variables: " + str(set(self.sample_specific_var_registry.keys())))
        model_var_set = {self._get_var_name(var) for var in self.vars}
        unannotated_vars = model_var_set \
            .difference(self.global_var_registry) \
            .difference(set(self.sample_specific_var_registry.keys()))
        assert len(unannotated_vars) == 0, \
            "The following model variables have not been registered " \
            "either as global or sample-specific: {0}".format(unannotated_vars)
