import pymc3 as pm
from .. import types

Operator = pm.operators.Operator
Inference = pm.Inference
MeanField = pm.MeanField


class KLThermal(Operator):
    """Kullback-Leibler divergence operator with finite temperature."""
    def __init__(self,
                 approx: pm.approximations.Approximation,
                 temperature: types.TensorSharedVariable):
        """Initializer.

        Args:
            approx: an instance of PyMC3 approximation
            temperature: a scalar shared theano tensor variable
        """
        super().__init__(approx)
        assert temperature is not None
        self.temperature = temperature

    def apply(self, f):
        z = self.input
        return self.temperature * self.logq_norm(z) - self.logp_norm(z)


class ADVIDeterministicAnnealing(Inference):
    """ADVI with deterministic annealing functionality.

    Note:
        Temperature is not updated automatically by this class. This task is delegated to the ADVI step
        function. This can be done by including a temperature update in `more_updates`; refer to
        `pymc3.opvi.ObjectiveFunction.step_function` for more information.

    """
    def __init__(self,
                 local_rv=None,
                 model=None,
                 cost_part_grad_scale=1,
                 scale_cost_to_minibatch=False,
                 random_seed=None, start=None,
                 temperature=None):

        assert temperature is not None, "Temperature (a scalar theano shared tensor) is not provided"
        super().__init__(
            KLThermal, MeanField, None,
            local_rv=local_rv,
            model=model,
            cost_part_grad_scale=cost_part_grad_scale,
            scale_cost_to_minibatch=scale_cost_to_minibatch,
            random_seed=random_seed,
            start=start,
            op_kwargs={'temperature': temperature})
