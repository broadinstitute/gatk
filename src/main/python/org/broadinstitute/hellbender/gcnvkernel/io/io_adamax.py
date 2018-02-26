import logging
import os

import numpy as np

from . import io_commons
from . import io_consts
from .. import config, types
from ..inference.fancy_optimizers import FancyAdamax

_logger = logging.getLogger(__name__)


class AdamaxStateWriter:
    """Writes the state of adamax optimizer to disk."""
    def __init__(self,
                 fancy_adamax: FancyAdamax,
                 output_path: str):
        self.fancy_adamax = fancy_adamax
        self.output_path = output_path

    def __call__(self):
        _logger.info("Writing adamax moment estimates...")
        io_commons.assert_output_path_writable(self.output_path)
        io_commons.write_gcnvkernel_version(self.output_path)

        mu_m = self.fancy_adamax.get_mu_m().get_value(borrow=True)
        np.save(os.path.join(self.output_path, "mu_" + io_consts.default_adamax_m_filename), mu_m)

        rho_m = self.fancy_adamax.get_rho_m().get_value(borrow=True)
        np.save(os.path.join(self.output_path, "rho_" + io_consts.default_adamax_m_filename), rho_m)

        mu_u = self.fancy_adamax.get_mu_u().get_value(borrow=True)
        np.save(os.path.join(self.output_path, "mu_" + io_consts.default_adamax_u_filename), mu_u)

        rho_u = self.fancy_adamax.get_rho_u().get_value(borrow=True)
        np.save(os.path.join(self.output_path, "rho_" + io_consts.default_adamax_u_filename), rho_u)

        if not self.fancy_adamax.disable_bias_correction:
            res = self.fancy_adamax.get_res_tensor().get_value(borrow=True)
            np.save(os.path.join(self.output_path, io_consts.default_adamax_res_filename), res)


class AdamaxStateReader:
    """Import the state of adamax optimizer from disk."""
    def __init__(self,
                 fancy_adamax: FancyAdamax,
                 input_path: str):
        self.fancy_adamax = fancy_adamax
        self.input_path = input_path

    @staticmethod
    def _assert_shape(imported_ndarray, expected_shared_tensor):
        assert imported_ndarray.shape == expected_shared_tensor.get_value(borrow=True).shape, \
            "The exported adamax state has a different shape (shape={0}) than the currently instantiated adamax " \
            "(shape={1}). This can occur if the exported adamax state corresponds to a run with different number of " \
            "samples, intervals, or model configuration (i.e. w/ or w/o explicit GC bias modeling or bias " \
            "factor discovery)".format(imported_ndarray.shape, expected_shared_tensor.get_value(borrow=True).shape)

    def __call__(self):
        _logger.info("Importing adamax moment estimates...")
        io_commons.check_gcnvkernel_version_from_path(self.input_path)

        imported_mu_m = np.load(
            os.path.join(self.input_path, "mu_" + io_consts.default_adamax_m_filename))
        self._assert_shape(imported_mu_m, self.fancy_adamax.get_mu_m())
        self.fancy_adamax.get_mu_m().set_value(imported_mu_m, borrow=config.borrow_numpy)

        imported_mu_u = np.load(
            os.path.join(self.input_path, "mu_" + io_consts.default_adamax_u_filename))
        self._assert_shape(imported_mu_u, self.fancy_adamax.get_mu_u())
        self.fancy_adamax.get_mu_u().set_value(imported_mu_u, borrow=config.borrow_numpy)

        imported_rho_m = np.load(
            os.path.join(self.input_path, "rho_" + io_consts.default_adamax_m_filename))
        self._assert_shape(imported_rho_m, self.fancy_adamax.get_rho_m())
        self.fancy_adamax.get_rho_m().set_value(imported_rho_m, borrow=config.borrow_numpy)

        imported_rho_u = np.load(
            os.path.join(self.input_path, "rho_" + io_consts.default_adamax_u_filename))
        self._assert_shape(imported_rho_u, self.fancy_adamax.get_rho_u())
        self.fancy_adamax.get_rho_u().set_value(imported_rho_u, borrow=config.borrow_numpy)

        if not self.fancy_adamax.disable_bias_correction:
            res_filename = os.path.join(self.input_path, io_consts.default_adamax_res_filename)
            if not os.path.exists(res_filename):
                _logger.warning("Could not find adamax bias correction residue tensor in \"{0}\"; setting the bias "
                                "correction residue to zero and proceeding".format(self.input_path))
                self.fancy_adamax.get_res_tensor().set_value(
                    np.asarray([0], dtype=types.floatX), borrow=config.borrow_numpy)
            else:
                imported_res = np.load(res_filename)
                self.fancy_adamax.get_res_tensor().set_value(imported_res, borrow=config.borrow_numpy)
