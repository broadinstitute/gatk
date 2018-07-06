from abc import abstractmethod
from collections import OrderedDict
from functools import partial
from typing import List

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
from pymc3.variational.updates import get_or_compute_grads

from .. import types
from ..io import io_commons
from ..models.fancy_model import GeneralizedContinuousModel


class FancyStochasticOptimizer:
    """The base class of stochastic optimizers equipped with the functionality of saving and loading the
    optimizer state to and from disk (e.g. for stateful optimizers such as ADAM and ADAMAX), and the
    possibility of utilizing the extra attributes of `GeneralizedContinuousModel` to perform structured
    parameter updates, e.g. updating only sample-specific variables while keeping global variables intact
    (see `FancyAdamax` for a concrete implementation).
    """
    @abstractmethod
    def get_optimizer(self,
                      model: GeneralizedContinuousModel=None,
                      approx: pm.MeanField=None):
        """

        Args:
            model: a generalized continuous PyMC3 model
            approx: an instance of PyMC3 mean-field approximation

        Returns:
            A callable function that upon providing `loss_or_grads` and `params`, returns an
            `OrderedDict` of shared theano tensor updates (for example, see `FancyAdamax.get_optimizer`).
        """
        raise NotImplementedError

    @staticmethod
    def get_call_kwargs(_locals_):
        _locals_ = _locals_.copy()
        _locals_.pop('loss_or_grads')
        _locals_.pop('params')
        return _locals_

    @abstractmethod
    def save(self, output_path: str):
        raise NotImplementedError

    @abstractmethod
    def load(self, input_path: str):
        raise NotImplementedError


class FancyAdamax(FancyStochasticOptimizer):
    """Adamax optimizer with saving/loading functionality and sample-specific-only update mode."""
    def __init__(self,
                 learning_rate: float = 0.002,
                 beta1: float = 0.9,
                 beta2: float = 0.999,
                 epsilon: float = 1e-8,
                 sample_specific_only: bool = False,
                 disable_bias_correction: bool = False):
        """Initializer.

        Args:
            learning_rate: learning rate
            beta1: first moment forgetting factor
            beta2: second moment forgetting factor
            epsilon: a small float for avoiding division-by-zero
            sample_specific_only: only update sample-specific variables (as specified in the generalized model)
            disable_bias_correction: disable moment estimation bias correction
        """
        self.learning_rate = learning_rate
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon
        self.sample_specific_only = sample_specific_only
        self.disable_bias_correction = disable_bias_correction

        # placeholder for first (m) and second (u) moments
        # in mean-field type approximation, ``mu`` and ``rho`` each have their own tensors
        # the list elements correspond to ``mu`` and ``rho`` moments, respectively
        self.m_tensors: List[types.TensorSharedVariable] = []
        self.u_tensors: List[types.TensorSharedVariable] = []

        # placeholder for the state of moment estimation bias corrector
        self.res_tensor: types.TensorSharedVariable = None

    def _assert_shared_tensors_available(self):
        m_u_available = len(self.m_tensors) == 2 and len(self.u_tensors) == 2
        res_available = self.disable_bias_correction or self.res_tensor is not None
        assert m_u_available and res_available, "Adamax tensors are not available yet"

    def get_mu_m(self):
        self._assert_shared_tensors_available()
        return self.m_tensors[0]

    def get_rho_m(self):
        self._assert_shared_tensors_available()
        return self.m_tensors[1]

    def get_mu_u(self):
        self._assert_shared_tensors_available()
        return self.u_tensors[0]

    def get_rho_u(self):
        self._assert_shared_tensors_available()
        return self.u_tensors[1]

    def get_res_tensor(self):
        return self.res_tensor

    @staticmethod
    def structured_adamax(loss_or_grads=None,
                          params=None,
                          model: GeneralizedContinuousModel=None,
                          approx: pm.MeanField=None,
                          learning_rate=0.002, beta1=0.9,
                          beta2=0.999, epsilon=1e-8,
                          sample_specific_only=False,
                          disable_bias_correction=False,
                          base_class: 'FancyAdamax' = None):
        """Adamax stochastic optimizer with partial sample-specific-only update functionality.

        Args:
            loss_or_grads: symbolic loss function or gradients
            params: variational parameter bundle
            model: an instance of generalized model
            approx: an instance of variational approximation for the model
            learning_rate: global learning rate
            beta1: first moment estimation forgetting factor
            beta2: second moment estimation forgetting factor
            epsilon: a small float to avoid division-by-zero
            sample_specific_only: only update parameters registered in the generalized model as sample-specific
            disable_bias_correction: disable moment estimation bias correction
            base_class: a reference to the base class to store a reference to the shared tensors (for I/O)

        Returns:
            returns the function itself if `loss_or_grads` and `params` are not given;
            otherwise, returns an ordered dict of shared tensor updates (to be used in pymc3 for compiling
            the step function)
        """
        if loss_or_grads is None and params is None:
            return partial(FancyAdamax.structured_adamax,
                           **FancyStochasticOptimizer.get_call_kwargs(locals()))
        elif loss_or_grads is None or params is None:
            raise ValueError('Please provide both `loss_or_grads` and `params` to get updates')
        assert model is not None, 'Please provide `model` to get updates'
        assert approx is not None, 'Please provide `approx` to get updates'

        all_grads = get_or_compute_grads(loss_or_grads, params)
        updates = OrderedDict()

        # indices of sample-specific vars
        if sample_specific_only:
            vmap_list = io_commons.get_var_map_list_from_mean_field_approx(approx)
            sample_specific_indices = []
            for vmap in vmap_list:
                if vmap.var in model.sample_specific_var_registry:
                    sample_specific_indices += [idx for idx in range(vmap.slc.start, vmap.slc.stop)]
            update_indices = th.shared(np.asarray(sample_specific_indices, dtype=np.int))
            num_dof = len(sample_specific_indices)

        # Using theano constant to prevent upcasting of float32
        one = tt.constant(1)

        if disable_bias_correction:
            a_t = learning_rate
        else:
            res_prev = th.shared(pm.theanof.floatX(beta1))
            res = beta1 * res_prev
            a_t = learning_rate / (one - res)
            updates[res_prev] = res
            if base_class is not None:
                base_class.res_tensor = res_prev

        for param, g_t in zip(params, all_grads):
            if sample_specific_only:
                g_t_view = g_t[update_indices]
                m_prev = th.shared(np.zeros((num_dof,), dtype=types.floatX),
                                   broadcastable=(False,))
                u_prev = th.shared(np.zeros((num_dof,), dtype=types.floatX),
                                   broadcastable=(False,))
            else:
                g_t_view = g_t
                value = param.get_value(borrow=True)
                m_prev = th.shared(np.zeros(value.shape, dtype=types.floatX),
                                   broadcastable=(False,))
                u_prev = th.shared(np.zeros(value.shape, dtype=types.floatX),
                                   broadcastable=(False,))

            # save a reference to m and u in the base class
            if base_class is not None:
                base_class.m_tensors.append(m_prev)
                base_class.u_tensors.append(u_prev)

            m_t = beta1 * m_prev + (one - beta1) * g_t_view
            u_t = tt.maximum(beta2 * u_prev, abs(g_t_view))
            step = a_t * m_t / (u_t + epsilon)

            if sample_specific_only:
                new_param = tt.inc_subtensor(param[update_indices], -step)
            else:
                new_param = param - step

            updates[m_prev] = m_t
            updates[u_prev] = u_t
            updates[param] = new_param

        return updates

    def get_optimizer(self,
                      model: GeneralizedContinuousModel=None,
                      approx: pm.MeanField=None):
        return FancyAdamax.structured_adamax(
            model=model,
            approx=approx,
            beta1=self.beta1,
            beta2=self.beta2,
            learning_rate=self.learning_rate,
            epsilon=self.epsilon,
            sample_specific_only=self.sample_specific_only,
            disable_bias_correction=self.disable_bias_correction,
            base_class=self)

    def save(self, output_path: str) -> None:
        """Saves the state of the optimizer to disk.

        Args:
            output_path: output path (must be writable directory)
        """
        from ..io import io_adamax  # lazy import to break import cycle
        io_adamax.AdamaxStateWriter(self, output_path)()

    def load(self, input_path: str):
        """Loads the state of the optimizer from disk.

        Args:
            input_path: input path (must be a readable directory)
        """
        from ..io import io_adamax  # lazy import to break import cycle
        io_adamax.AdamaxStateReader(self, input_path)()

    def initialize_state_from_instance(self, instance: 'FancyAdamax'):
        self.get_mu_m().set_value(instance.get_mu_m().get_value())
        self.get_mu_u().set_value(instance.get_mu_u().get_value())
        self.get_rho_m().set_value(instance.get_rho_m().get_value())
        self.get_rho_u().set_value(instance.get_rho_u().get_value())
        if not self.disable_bias_correction:
            self.get_res_tensor().set_value(instance.get_res_tensor().get_value())