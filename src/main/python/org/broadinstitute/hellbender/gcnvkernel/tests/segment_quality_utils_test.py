import unittest

import numpy as np
import scipy.stats as stats


class HMMSegmentationTestDataGenerator:
    """Generates test data for a Markov chain with one true segment. The emission probabilities are drawn from a
    negative binomial distribution:

        p_emission(s_trial_t) ~ NB(x_t | mu_t[s_trial_t], r_t)

    where `t` indicates the position, `s_trial_t` is the emission trial state, `mu_t` is given as:

        mu_t[s] = depth * s * (1 - error_rate) + error_rate * depth,

    where `depth` is the emission depth (higher depth implies stronger signal), `r_t` is the number of failures
    until the experiment is stopped and is related to "over-dispersion" `psi_t` as:

        1 / r_t = exp(psi_t) - 1,

    finally, `x_t` is the observed data which we set to the mode of the distribution (calculated for s_truth_t)
    with an additive white noise:

        x_t = max(floor(mu_t[s_truth_t] + noise_std * normal_random), 0).

    For noise_std = 0, it is guaranteed that `p_emission(s_trial_t)` is maximized for s_trial_t = s_truth_t.

    In order to test situations where event start/end determination qualities are different, we allow the
    over-dispersion to linear change from `start_over_dispersion` at t = 0 to `end_over_dispersion` at
    t = num_positions - 1. Note that higher over-dispersion implies reduced emission strength.

    The HMM which we use for detecting the event is specified with a prior that is peaked at the `baseline_state`
    with the respect of the probability mass distributed equality among other states:

        prior_prob(s) = | (1 - prob_event)               if s = baseline_state,
                        | prob_event / (num_states - 1)  if s != baseline_state.

    The transition matrix is the same at all positions:

        prob_trans_ab = prior_prob(b) * [1 - exp(-1 / coherence_length)] +
                        delta(a, b) * exp(-1 / coherence_length).
    """

    def __init__(self,
                 emission_depth: int = 100.0,
                 error_rate: float = 0.01,
                 num_positions: int = 100,
                 event_start: int = 10,
                 event_end: int = 20,
                 event_state: int = 1,
                 start_over_dispersion: float = 0.001,
                 end_over_dispersion: float = 0.001,
                 baseline_state: int = 2,
                 num_states: int = 5,
                 prior_prob_event: float = 1e-3,
                 coherence_length: float = 0.1):
        self.emission_depth = emission_depth
        self.error_rate = error_rate
        self.num_positions = num_positions
        self.event_start = event_start
        self.event_end = event_end
        self.event_state = event_state
        self.start_over_dispersion = start_over_dispersion
        self.end_over_dispersion = end_over_dispersion
        self.baseline_state = baseline_state
        self.num_states = num_states
        self.prior_prob_event = prior_prob_event
        self.coherence_length = coherence_length

    def get_log_trans_tss(self) -> np.ndarray:
        pass

    def get_log_prior_s(self) -> np.ndarray:
        pass

    def get_log_emission_ts(self) -> np.ndarray:
        pass

    def get_x_t(self) -> np.ndarray:
        pass

    def get_s_truth_t(self) -> np.ndarray:
        pass

    @staticmethod
    def _psi_to_alpha(psi: float):
        return 1. / (np.exp(psi) - 1.)

    @staticmethod
    def _get_neg_binom_lpmf(n: int, mu: float, psi: float):
        alpha = HMMSegmentationTestDataGenerator._psi_to_alpha(psi)
        return stats.nbinom.logpmf(n, n=alpha, p=alpha/(mu+alpha))

    @staticmethod
    def _get_neg_binom_lpmf_mode(mu: float, psi: float):
        alpha = HMMSegmentationTestDataGenerator._psi_to_alpha(psi)
        mode = np.floor(mu * (alpha - 1.) / alpha)
        return HMMSegmentationTestDataGenerator._get_neg_binom_lpmf(mode, mu, psi)


class TestHMMSegmentationQualityCalculator(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestHMMSegmentationQualityCalculator, self).__init__(*args, **kwargs)

    def test_quality_increased_with_compatible_prior(self):
        pass

    def test_quality_increased_with_larger_events(self):
        pass

    def test_quality_start_end(self):
        pass

    def test_quality_all_called(self):
        pass

    def test_quality_some_called(self):
        pass
