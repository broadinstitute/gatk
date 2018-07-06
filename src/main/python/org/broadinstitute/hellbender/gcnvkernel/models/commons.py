import numpy as np
import logging
import theano.tensor as tt
import pymc3 as pm
import theano as th
import pymc3.distributions.dist_math as pm_dist_math
from typing import Tuple, Generator

_logger = logging.getLogger(__name__)

_log_2_pi = 1.837877066409345  # np.log(2 * np.pi)
_10_inv_log_10 = 4.342944819032518  # 10 / np.log(10)


def get_normalized_prob_vector(prob_vector: np.ndarray, prob_sum_tol: float) -> np.ndarray:
    """Normalizes the probability vector of a categorical RV to unity.

    Args:
        prob_vector: input probability vector a categorical RV of choice
        prob_sum_tol: tolerated amount of deviation from unity before performing normalization

    Returns:
        A new and normalized probability vector if it deviates from unity more than `prob_sum_tol`
        Otherwise, `prob_vector` is returned unchanged.
    """
    assert all(prob_vector >= 0), "Probabilities must be non-negative"
    prob_sum = np.sum(prob_vector)
    if np.abs(prob_sum - 1.0) < prob_sum_tol:
        return prob_vector
    else:
        _logger.warning("The given probability vector ({0}) is not normalized to unity within the provided "
                        "tolerance ({1}); sum = {2}. Normalization was enforced. However, please check the inputs "
                        "for unintentional errors.".format(
            prob_vector, prob_sum_tol, prob_sum))
        return prob_vector / prob_sum


def negative_binomial_logp(mu, alpha, value):
    """Generates symbolic negative binomial logp.

    Args:
        mu: negative binomial mean tensor
        alpha: negative binomial over-dispersion tensor
        value: negative binomial counts

    Note:
        `mu`, `alpha` and `value` must be either shape-compatible or be commensurately broadcastable.

    Returns:
        symbolic negative binomial logp
    """
    return pm_dist_math.bound(pm_dist_math.binomln(value + alpha - 1, value)
                              + pm_dist_math.logpow(mu / (mu + alpha), value)
                              + pm_dist_math.logpow(alpha / (mu + alpha), alpha),
                              mu > 0, value >= 0, alpha > 0)


def negative_binomial_gaussian_approx_logp(mu, alpha, value):
    """Generates symbolic Gaussian approximation to negative binomial logp.

    Args:
        mu: negative binomial mean tensor
        alpha: negative binomial over-dispersion tensor
        value: negative binomial counts

    Note:
        `mu`, `alpha` and `value` must be either shape-compatible or be commensurately broadcastable.

    Returns:
        symbolic approximate negative binomial logp
    """
    tau = alpha / (mu * (alpha + mu))  # precision
    return pm_dist_math.bound(0.5 * (tt.log(tau) - _log_2_pi - tau * tt.square(value - mu)),
                              mu > 0, value >= 0, alpha > 0)


# todo
def negative_binomial_smart_approx_logp(mu, alpha, value):
    """Generates symbolic negative binomial logp with conditional switching to Gaussian approximation
    if the approximation is valid.

    Args:
        mu: negative binomial mean tensor
        alpha: negative binomial over-dispersion tensor
        value: negative binomial counts

    Returns:
        symbolic approximate negative binomial logp
    """
    raise NotImplementedError


def centered_heavy_tail_logp(mu, value):
    """This distribution is obtained by taking X ~ Exp and performing a Bose transformation
    Y = (exp(X) - 1)^{-1}. The result is:

        p(y) = (1 + 2 \mu) y^{2\mu} (1 + y)^{-2(1 + \mu)}

    It is a heavy-tail distribution with non-existent first moment.

    Args:
        mu: exponential parameter of X
        value: values of Y

    Returns:
        symbolic logp
    """
    return pm_dist_math.bound(tt.log(1.0 + 2.0 * mu) + 2.0 * mu * tt.log(value)
                              - 2.0 * (1.0 + mu) * tt.log(1.0 + value),
                              mu >= 0, value > 0)


def safe_logaddexp(a, b):
    """Symbolic log(exp(a) + exp(b)). The edge case where `a` - `b` is undefined is handled by
    setting the difference to 0. This occurs if both `a` and `b` are +inf or -inf.

    Returns:
        symbolic log(exp(a) + exp(b))
    """
    diff = b - a
    safe_diff = tt.switch(tt.isnan(diff), 0, diff)
    return tt.switch(safe_diff >= 0,
                     b + tt.log1p(tt.exp(-safe_diff)),
                     a + tt.log1p(tt.exp(safe_diff)))


def get_jensen_shannon_divergence(log_p_1, log_p_2):
    """Symbolic Jensen-Shannon distance (symmetric KL divergence) between two discrete distributions.

    Args:
        log_p_1: first discrete probability distribution in log space
        log_p_2: second discrete probability distribution in log space

    Returns:
        Symbolic Jensen-Shannon distance
    """
    p_1 = tt.exp(log_p_1)
    p_2 = tt.exp(log_p_2)
    diff_12 = p_1 - p_2
    log_diff_12 = log_p_1 - log_p_2
    safe_log_diff_12 = tt.switch(tt.isnan(log_diff_12), 0, log_diff_12)
    return 0.5 * tt.sum(diff_12 * safe_log_diff_12, axis=-1)


def get_hellinger_distance(log_p_1, log_p_2):
    """Symbolic Hellinger distance between two discrete distributions.

    Args:
        log_p_1: first discrete probability distribution in log space
        log_p_2: second discrete probability distribution in log space

    Returns:
        Symbolic Hellinger distance
    """
    p_1 = tt.exp(log_p_1)
    p_2 = tt.exp(log_p_2)
    return tt.sqrt(tt.sum(tt.square(tt.sqrt(p_1) - tt.sqrt(p_2)), axis=-1)) / tt.sqrt(2)


def perform_genotyping(log_p: np.ndarray) -> Tuple[int, float]:
    """Takes a vector of probabilities in log space and perform genotyping.

    Args:
        log_p: a vector probabilities in log scape

    Note:
        log_p must be properly normalized, i.e. np.sum(np.exp(log_p)) == 1
        (this is not explicitly asserted)

    Returns:
        A tuple (most likely genotype index, phred-scaled genotyping quality)
    """
    assert log_p.ndim == 1, "The log_p array is not a vector"
    assert log_p.size >= 2, "At least two states are required for genotyping"
    sorted_log_p = sorted(enumerate(log_p), key=lambda x: -x[1])
    max_likely_genotype_idx = sorted_log_p[0][0]
    phred_genotype_quality = _10_inv_log_10 * (sorted_log_p[0][1] - sorted_log_p[1][1])
    return max_likely_genotype_idx, phred_genotype_quality


def stochastic_node_mean_symbolic(approx: pm.MeanField, node, size=100,
                                  more_replacements=None, shape=None, dtype=None):
    """Symbolic mean of a given PyMC3 stochastic node with respect to a given variational
    posterior approximation.

    Args:
        approx: an instance of PyMC3 approximation
        node: stochastic node
        size: the number of samples to use for calculating the mean
        more_replacements: (optional) an ordered dictionary of node replacements to be
            applied to the computational graph before sampling
        shape: (optional) shape of the node
        dtype: (optional) dtype of the node

    Raises:
        ValueError: If `node.tag.test_value` is not present and `shape` and `dtype` are
            not provided

    Returns:
        Symbolic approximate mean of the stochastic node
    """

    assert size > 0

    if shape is not None and dtype is not None:
        cum_sum = tt.zeros(shape, dtype)
    elif hasattr(node.tag, 'test_value') and node.tag.test_value is not None:
        cum_sum = tt.zeros(node.tag.test_value.shape, node.tag.test_value.dtype)
    else:
        raise ValueError("Can not determine the shape of the node to be sampled")

    if more_replacements is not None:
        node = th.clone(node, more_replacements, strict=False)
    posterior_samples = approx.random(size)
    node = approx.to_flat_input(node)

    def add_sample_to_cum_sum(posterior_sample, _cum_sum):
        new_sample = th.clone(node, {approx.input: posterior_sample}, strict=False)
        return _cum_sum + tt.patternbroadcast(new_sample, _cum_sum.broadcastable)

    outputs, _ = th.scan(add_sample_to_cum_sum,
                         sequences=[posterior_samples],
                         outputs_info=[cum_sum],
                         n_steps=size)

    return outputs[-1] / size


def get_sampling_generator_for_model_approximation(model_approx: pm.MeanField, model_var_name: str,
                                                   num_samples: int = 250) -> Generator:
    """Get a generator that returns samples of a precomputed model approximation for a specific variable in that model

    Args:
        model_approx: an instance of PyMC3 mean-field approximation
        model_var_name: a stochastic node in the model
        num_samples: number of samples to draw

    Returns:
        A generator that will yield `num_samples` samples from an approximation to a posterior
    """
    return (model_approx.sample()[model_var_name] for _ in range(num_samples))
