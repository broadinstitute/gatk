import collections
import unittest

import numpy as np
import scipy.stats as stats
from typing import Optional
from gcnvkernel.models.theano_hmm import TheanoForwardBackward
from gcnvkernel.postprocess.segment_quality_utils import HMMSegmentationQualityCalculator


# default values for the test data generator
DEFAULT_EMISSION_DEPTH = 100.0
DEFAULT_ERROR_RATE = 0.01
DEFAULT_NUM_POSITIONS = 100
DEFAULT_START_OVER_DISPERSION = 0.001
DEFAULT_END_OVER_DISPERSION = 0.001
DEFAULT_NOISE_STD = 0.
DEFAULT_BASELINE_STATE = 2
DEFAULT_NUM_STATES = 5
DEFAULT_PRIOR_EVENT_PROB = 1e-3
DEFAULT_COHERENCE_LENGTH = 0.1
DEFAULT_RANDOM_SEED = 1984

# a few typical intervals in [0, 100) for generating test events
TEST_INTERVALS = [(0, 10), (89, 99), (10, 30), (25, 75), (70, 90)]
DEFAULT_NON_BASELINE_STATES = [0, 1, 3, 4]


class HMMSegmentationTestData:
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
                 event_start: int,
                 event_end: int,
                 event_state: int,
                 emission_depth: float = DEFAULT_EMISSION_DEPTH,
                 error_rate: float = DEFAULT_ERROR_RATE,
                 num_positions: int = DEFAULT_NUM_POSITIONS,
                 start_over_dispersion: float = DEFAULT_START_OVER_DISPERSION,
                 end_over_dispersion: float = DEFAULT_END_OVER_DISPERSION,
                 noise_std: float = DEFAULT_NOISE_STD,
                 baseline_state: int = DEFAULT_BASELINE_STATE,
                 num_states: int = DEFAULT_NUM_STATES,
                 prior_prob_event: float = DEFAULT_PRIOR_EVENT_PROB,
                 prior_s: Optional[np.ndarray] = None,
                 coherence_length: float = DEFAULT_COHERENCE_LENGTH,
                 random_seed: int = DEFAULT_RANDOM_SEED):
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

        # generate s_truth_t
        self.s_truth_t = baseline_state * np.ones((num_positions,), dtype=np.int)
        self.s_truth_t[event_start:event_end + 1] = event_state

        # generate psi_t
        self.psi_t = np.asarray([start_over_dispersion +
                                 t * (end_over_dispersion - start_over_dispersion) / (num_positions - 1)
                                 for t in range(num_positions)])

        # generate x_t
        np.random.seed(random_seed)
        noise_t = noise_std * np.random.randn(num_positions)
        mu_for_s_truth_t = np.asarray([self._get_mu(s) for s in self.s_truth_t])
        mode_for_s_truth_t = mu_for_s_truth_t * (2. - np.exp(self.psi_t))
        self.x_t = np.maximum(np.floor(mode_for_s_truth_t + noise_t), 0.).astype(np.int)

        # generate log emission probabilities
        self.log_emission_ts = np.asarray(
            [[self._get_neg_binom_lpmf(x, self._get_mu(s_trial), psi)
              for s_trial in range(num_states)] for x, psi in zip(self.x_t, self.psi_t)])

        # generate HMM specs
        if prior_s is None:
            prior_s = np.ones((num_states,)) * prior_prob_event / (num_states - 1)
            prior_s[baseline_state] = 1. - prior_prob_event
        prob_stay = np.exp(-1. / coherence_length)
        trans_ss = prior_s[np.newaxis, :] * (1. - prob_stay) + prob_stay * np.eye(num_states)
        self.log_prior_s = np.log(prior_s)
        self.log_trans_tss = np.tile(np.log(trans_ss), (num_positions - 1, 1, 1))

    def get_log_trans_tss(self) -> np.ndarray:
        return self.log_trans_tss

    def get_log_prior_s(self) -> np.ndarray:
        return self.log_prior_s

    def get_log_emission_ts(self) -> np.ndarray:
        return self.log_emission_ts

    def get_x_t(self) -> np.ndarray:
        return self.x_t

    def get_s_truth_t(self) -> np.ndarray:
        return self.s_truth_t

    def _get_mu(self, s):
        return self.emission_depth * s * (1. - self.error_rate) + self.error_rate * self.emission_depth

    @staticmethod
    def _psi_to_alpha(psi: float):
        return 1. / (np.exp(psi) - 1.)

    @staticmethod
    def _get_neg_binom_lpmf(n: int, mu: float, psi: float):
        alpha = HMMSegmentationTestData._psi_to_alpha(psi)
        return stats.nbinom.logpmf(n, n=alpha, p=alpha/(mu+alpha))

    @staticmethod
    def _get_neg_binom_lpmf_mode(mu: float, psi: float):
        mode = np.floor(mu * (2. - np.exp(psi)))
        return HMMSegmentationTestData._get_neg_binom_lpmf(mode, mu, psi)


class TestHMMSegmentationQualityCalculator(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestHMMSegmentationQualityCalculator, self).__init__(*args, **kwargs)
        self.fb_theano = TheanoForwardBackward(include_alpha_beta_output=True)

    SegmentQualityBundle = collections.namedtuple(
        'SegmentQualityBundle',
        ['quality_some_called', 'quality_all_called', 'quality_start', 'quality_end'])

    def _get_test_data_segment_quality(self,
                                       test_data: HMMSegmentationTestData,
                                       segment_start: int,
                                       segment_end: int,
                                       segment_call: int):
        """Calculates the segmentation quality metrics for a given segmentation and test data."""
        fb_result = self.fb_theano.perform_forward_backward(
            test_data.log_prior_s, test_data.log_trans_tss, test_data.log_emission_ts)
        segment_quality_calculator = HMMSegmentationQualityCalculator(
            test_data.log_emission_ts, test_data.log_trans_tss, fb_result)

        return TestHMMSegmentationQualityCalculator.SegmentQualityBundle(
            segment_quality_calculator.get_segment_quality_some_called(segment_start, segment_end, segment_call),
            segment_quality_calculator.get_segment_quality_all_called(segment_start, segment_end, segment_call),
            segment_quality_calculator.get_segment_quality_start(segment_start, segment_call),
            segment_quality_calculator.get_segment_quality_end(segment_end, segment_call))

    def test_quality_increased_with_better_prior(self):
        """Tests that all segmentation quality factors increase if a favorable prior is used."""
        for event_start, event_end in TEST_INTERVALS:
            for event_state in DEFAULT_NON_BASELINE_STATES:

                test_data_default_prior = HMMSegmentationTestData(
                    event_start, event_end, event_state)

                # make a event_state-favoring prior
                favorable_prior_s = np.ones((DEFAULT_NUM_STATES,))
                favorable_prior_s[event_state] = 10.0
                favorable_prior_s /= np.sum(favorable_prior_s)

                test_data_favorable_prior = HMMSegmentationTestData(
                    event_start, event_end, event_state, prior_s=favorable_prior_s)

                qualities_default_prior = self._get_test_data_segment_quality(
                    test_data_default_prior, event_start, event_end, event_state)
                qualities_favorable_prior = self._get_test_data_segment_quality(
                    test_data_favorable_prior, event_start, event_end, event_state)

                self.assertTrue(qualities_favorable_prior.quality_some_called >
                                qualities_default_prior.quality_some_called)
                self.assertTrue(qualities_favorable_prior.quality_all_called >
                                qualities_default_prior.quality_all_called)
                self.assertTrue(qualities_favorable_prior.quality_start >
                                qualities_default_prior.quality_start)
                self.assertTrue(qualities_favorable_prior.quality_end >
                                qualities_default_prior.quality_end)

    def test_quality_start_end_decreased_with_bad_breakpoints(self):
        """Tests that poorly determined breakpoints reduces all quality metrics."""
        for event_start, event_end in TEST_INTERVALS:
            for event_state in DEFAULT_NON_BASELINE_STATES:
                test_data = HMMSegmentationTestData(event_start, event_end, event_state)
                qualities_exact_event = self._get_test_data_segment_quality(
                    test_data, event_start, event_end, event_state)

                for left_jitter in range(-2, 2):
                    for right_jitter in range(-2, 2):

                        test_segment_start = min(max(0, left_jitter + event_start), DEFAULT_NUM_POSITIONS - 1)
                        test_segment_end = min(max(0, right_jitter + event_end), DEFAULT_NUM_POSITIONS - 1)

                        qualities_jittery_event = self._get_test_data_segment_quality(
                            test_data, test_segment_start, test_segment_end, event_state)

                        if test_segment_start == event_start:
                            self.assertTrue(qualities_exact_event.quality_start == qualities_jittery_event.quality_start)
                        else:
                            self.assertTrue(qualities_exact_event.quality_start > qualities_jittery_event.quality_start)

                        if test_segment_end == event_end:
                            self.assertTrue(qualities_exact_event.quality_end == qualities_jittery_event.quality_end)
                        else:
                            self.assertTrue(qualities_exact_event.quality_end > qualities_jittery_event.quality_end)

    def test_quality_decreased_with_bad_call(self):
        """Tests that poorly determined call reduces all quality metrics."""
        pass

    def test_quality_start_end_symmetry(self):
        """Tests that start and end qualities are symmetric in symmetric situations."""
        pass

    def test_quality_single_point_event(self):
        """Tests that for a single-point event, all quality metrics are the same."""
        pass

    def test_quality_decreased_for_small_coherence_length(self):
        """Tests that a small coherence length (compared to the event length-scale) leads to lower qualities."""
        test_data_low_coherence_length = HMMSegmentationTestData(
            25, 75, 1, coherence_length=1.0, emission_depth=10.0)
        test_data_compatible_coherence_length = HMMSegmentationTestData(
            25, 75, 1, coherence_length=10.0, emission_depth=10.0)
        low_coherence_length_qualities = self._get_test_data_segment_quality(
            test_data_low_coherence_length, 25, 75, 1)
        compatible_coherence_length_qualities = self._get_test_data_segment_quality(
            test_data_compatible_coherence_length, 25, 75, 1)

        self.assertTrue(compatible_coherence_length_qualities.quality_some_called >
                        low_coherence_length_qualities.quality_some_called)
        self.assertTrue(compatible_coherence_length_qualities.quality_all_called >
                        low_coherence_length_qualities.quality_all_called)
        self.assertTrue(compatible_coherence_length_qualities.quality_start >
                        low_coherence_length_qualities.quality_start)
        self.assertTrue(compatible_coherence_length_qualities.quality_end >
                        low_coherence_length_qualities.quality_end)

    def test_quality_start_and_end_consistent_on_slope(self):
        """Tests that for an event located in a region with over-dispersion gradient, the breakpoint qualities
        are consistently higher/lower depending on the sign of the gradient."""
        pass

    def test_quality_some_and_all_called_decreased_with_noise(self):
        """Tests that noisy data decreases all quality metrics."""
        pass

    def test_quality_increased_with_depth(self):
        """Tests that higher emission depth implies higher qualities."""
        pass

    def test_quality_decreased_with_error_rate(self):
        """Tests that higher error rate implies lower qualities."""
        pass
