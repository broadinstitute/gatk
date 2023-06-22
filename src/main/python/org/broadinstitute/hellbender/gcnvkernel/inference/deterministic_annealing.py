import pymc3 as pm
from .. import types
from pymc3.variational.approximations import MeanFieldGroup
from pymc3.variational import Group
from pymc3.variational import Approximation
from pymc3.variational.opvi import GroupError
import theano
from pymc3.model import modelcontext
import numpy as np

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
        return (self.temperature * self.logq_norm - self.logp_norm)[0]


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
        approx = DeterministicMeanField(local_rv=local_rv,
                                        model=model,
                                        cost_part_grad_scale=cost_part_grad_scale,
                                        scale_cost_to_minibatch=scale_cost_to_minibatch,
                                        random_seed=random_seed,
                                        start=start)
        super().__init__(
            KLThermal, approx, None, temperature=temperature)

class DeterministicApproximation(Approximation):

    def __init__(self, groups, model=None):
        self._scale_cost_to_minibatch = theano.shared(np.int8(1))
        model = modelcontext(model)
        if not model.free_RVs:
            raise TypeError('Model does not have FreeRVs')
        self.groups = list()
        seen = set()
        rest = None
        for g in groups:
            if g.group is None:
                if rest is not None:
                    raise GroupError('More than one group is specified for '
                                     'the rest variables')
                else:
                    rest = g
            else:
                if set(g.group) & seen:
                    raise GroupError('Found duplicates in groups')
                seen.update(g.group)
                self.groups.append(g)
        if set(model.free_RVs) - seen:
            if rest is None:
                raise GroupError('No approximation is specified for the rest variables')
            else:
                unseen_free_RVs = [var for var in model.free_RVs if var not in seen]
                rest.__init_group__(unseen_free_RVs)
                self.groups.append(rest)
        self.model = model



class SingleGroupApproximation(DeterministicApproximation):
    """Base class for Single Group Approximation"""
    _group_class = None

    def __init__(self, *args, **kwargs):
        local_rv = kwargs.get('local_rv')
        groups = [self._group_class(None, *args, **kwargs)]
        if local_rv is not None:
            groups.extend([Group([v], params=p, local=True, model=kwargs.get('model'))
                           for v, p in local_rv.items()])
        super(SingleGroupApproximation, self).__init__(groups, model=kwargs.get('model'))

    def __getattr__(self, item):
        return getattr(self.groups[0], item)

    def __dir__(self):
        d = set(super(SingleGroupApproximation, self).__dir__())
        d.update(self.groups[0].__dir__())
        return list(sorted(d))
    
class DeterministicMeanField(SingleGroupApproximation):
    __doc__ = """**Single Group Mean Field Approximation**

    """ + str(MeanFieldGroup.__doc__)
    _group_class = MeanFieldGroup
