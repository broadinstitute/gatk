import numpy as np
import pymc3 as pm
from typing import Callable, Dict
from ..io import io_commons


class VariationalParameterTrackerConfig:
    """Configuration for `VariationalParameterTracker`."""
    def __init__(self):
        self.var_names = []
        self.inv_trans_list = []
        self.inv_trans_var_names = []

    def add(self,
            var_name: str,
            inv_trans: Callable[[np.ndarray], np.ndarray],
            inv_trans_var_name: str):
        """Adds a new parameter to be tracked.

        Args:
            var_name: name of the variable (must be a free variable name in the model)
            inv_trans: inverse transformation
            inv_trans_var_name: name of the variable after inverse transformation

        Returns:
            None
        """
        self.var_names.append(var_name)
        self.inv_trans_list.append(inv_trans)
        self.inv_trans_var_names.append(inv_trans_var_name)


class VariationalParameterTracker:
    """Keeps track of specified variational parameter."""
    def __init__(self, tracker_config: VariationalParameterTrackerConfig):
        self._var_trans_dict = {}
        self.tracked_var_values_dict = {}
        for var_name, inv_trans, inv_trans_var_name in zip(
                tracker_config.var_names,
                tracker_config.inv_trans_list,
                tracker_config.inv_trans_var_names):
            self._var_trans_dict[var_name] = (inv_trans, inv_trans_var_name)
            self.tracked_var_values_dict[inv_trans_var_name] = []

    def _extract_param_mean(self, approx: pm.approximations.MeanField) -> Dict[str, np.ndarray]:
        mu_flat_view = approx.mean.get_value(borrow=True)
        vmap_list = io_commons.get_var_map_list_from_mean_field_approx(approx)
        out = dict()
        for vmap in vmap_list:
            param_name = vmap.var
            if param_name in self._var_trans_dict.keys():
                bare_param_mean = mu_flat_view[vmap.slc].reshape(vmap.shp).astype(vmap.dtyp)
                inv_trans = self._var_trans_dict[param_name][0]
                inv_trans_param_name = self._var_trans_dict[param_name][1]
                if inv_trans is None:
                    out[inv_trans_param_name] = bare_param_mean
                else:
                    out[inv_trans_param_name] = inv_trans(bare_param_mean)
        return out

    def record(self, approx, _loss, _i):
        out = self._extract_param_mean(approx)
        for key in self.tracked_var_values_dict.keys():
            self.tracked_var_values_dict[key].append(out[key])

    __call__ = record

    def clear(self):
        for key in self.tracked_var_values_dict.keys():
            self.tracked_var_values_dict[key] = []

    def __getitem__(self, key):
        return self.tracked_var_values_dict[key]
