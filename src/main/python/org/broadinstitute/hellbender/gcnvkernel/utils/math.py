import sys

import numpy as np

# 10 / ln(10)
INV_LN_10_TIMES_10 = 4.342944819032518

# ln(1/2)
LN_HALF = -0.6931471805599453

# maximum possible phred-scaled value
MAX_PHRED = -INV_LN_10_TIMES_10 * np.log(sys.float_info.min)


def logp_to_phred(logp: float, complement: bool = False) -> float:
    """Converts probabilities from natural log scale to phred scale.

    Args:
        logp: a probability in the natural log scale
        complement: if True, returns the result for the complement of 'logp'

    Returns:
        phred-scaled probability
    """
    logp_zero_capped = min(0., logp)
    value = -INV_LN_10_TIMES_10 * (logp_zero_capped if not complement else logp_complement(logp_zero_capped))
    return min(value, MAX_PHRED)


def logp_complement(logp: float) -> float:
    """Calculates the complement of a probability in the natural log scale:

        log(1 - exp(logp)),

    in a numerically stable fashion.

    Args:
        logp: a probability in the natural log scale

    Returns:
        complement of the the probability in the natural log scale
    """
    logp_zero_capped = min(0., logp)
    if logp_zero_capped >= LN_HALF:
        return np.log(-np.expm1(logp_zero_capped))
    else:
        return np.log1p(-np.exp(logp_zero_capped))


def logsumexp_double_complement(a: np.ndarray, rel_tol: float = 1e-3) -> float:
    """Calculates the following expression in a numerically stable fashion:

        log(1 - (1 - exp(a_0)) x (1 - exp(a_1)) x ...)

    where a_i are the entries of `a` and assumed to be non-positive. The algorithm is as follows:

    We define:

        exp(x_n) = 1 - \prod_{i=0}^n (1 - exp(a_n)),

    Thus, we have x_0 = a_0 and the recursion relation:

        exp(x_{n+1}) = exp(x_n) + exp(b_{n+1}),

    where

        b_{n+1} = a_{n+1} + log(1 - exp(x_n)).

    We sort `a` in the descending order and update `x` term by term. It is easy to show that x_{n} is monotonically
    increasing and that |x_{N} - x_{n}| < (N - n) |x_{n} - x_{n-1}|. We use the last inequality to bound the error
    for early stopping.

    Args:
        a: a float array
        rel_tol: relative error tolerance for early stopping of calculation

    Returns:
        a float scalar
    """
    try:
        assert isinstance(a, np.ndarray)
        a = np.asarray(a.copy(), dtype=np.float)
    except AssertionError:
        try:
            a = np.asarray(a, dtype=np.float)
        except ValueError:
            raise ValueError("The input argument must be castable to a float ndarray.")
    assert len(a) > 0
    assert 0. <= rel_tol < 1.0

    # enforce all entries of a to be negative or zero
    a[a > 0.] = 0.

    if len(a) == 1:
        return np.asscalar(a)
    else:
        a = np.sort(a.flatten())[::-1]
        x = a[0]
        sz = len(a)
        for i, entry in enumerate(a[1:]):
            x_new = np.logaddexp(x, entry + logp_complement(x))
            if np.abs(x_new - x) * (sz - i - 1) < rel_tol * np.abs(x):
                return x_new
            else:
                x = x_new
        return x
