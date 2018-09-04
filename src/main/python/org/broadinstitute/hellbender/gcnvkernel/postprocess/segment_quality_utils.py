import logging
from typing import Dict, List

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
from scipy.misc import logsumexp

from ..utils.math import logsumexp_double_complement, logp_to_phred

_logger = logging.getLogger(__name__)


class HMMSegmentationQualityCalculator:
    """Calculates quality metrics for hidden state segments for a given HMM.

    Note:
        The initializer requires the emission and transition probabilities, as well as the forward
        and backward tables and the log posterior probability.
    """

    def __init__(self,
                 log_emission_tc: np.ndarray,
                 log_trans_tcc: np.ndarray,
                 alpha_tc: np.ndarray,
                 beta_tc: np.ndarray,
                 log_posterior_prob_tc: np.ndarray,
                 log_data_likelihood: float):
        """Initializer.

        Args:
            log_emission_tc: log copy-number emission matrix
            log_trans_tcc: log copy-number transition tensor
            alpha_tc: forward log likelihood matrix (from forward-backward algorithm)
            beta_tc: backward log likelihood matrix (from forward-backward algorithm)
            log_posterior_prob_tc: log copy-number posterior matrix (from forward-backward algorithm)
            log_data_likelihood: log data likelihood (from forward-backward algorithm)
        """
        assert isinstance(log_emission_tc, np.ndarray)
        assert log_emission_tc.ndim == 2
        self.num_sites, self.num_states = log_emission_tc.shape

        assert isinstance(log_trans_tcc, np.ndarray)
        assert log_trans_tcc.shape == (self.num_sites - 1, self.num_states, self.num_states)

        assert isinstance(alpha_tc, np.ndarray)
        assert alpha_tc.shape == (self.num_sites, self.num_states)

        assert isinstance(beta_tc, np.ndarray)
        assert beta_tc.shape == (self.num_sites, self.num_states)

        assert isinstance(log_posterior_prob_tc, np.ndarray)
        assert log_posterior_prob_tc.shape == (self.num_sites, self.num_states)

        self.log_emission_tc = log_emission_tc
        self.log_trans_tcc = log_trans_tcc
        self.alpha_tc = alpha_tc
        self.beta_tc = beta_tc
        self.log_posterior_prob_tc = log_posterior_prob_tc
        self.log_data_likelihood = log_data_likelihood

        self.all_states_set = set(range(self.num_states))
        self.all_states_list = list(range(self.num_states))
        self.leave_one_out_state_lists: Dict[int, List[int]] = {
            left_out_state: [state for state in range(self.num_states) if state != left_out_state]
            for left_out_state in range(self.num_states)}

    @staticmethod
    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_constrained_path_logp_theano_func() -> th.compile.function_module.Function:
        """Returns a theano function that calculates the log posterior probability of hidden state paths composed
        from a subset X of all hidden states S.

        More explicitly, this function calculates the logp of all paths constrained to the hidden-state set
        X for t_0 <= t <= t_N but unconstrained for before (t < t_0) and after (t > t_N):

        unconstrained           constrained           unconstrained
             t < t_0 |    t_0    t_1    ...    t_N   |  t > t_N
             c in S  |  c in X  c in X       c in X  |  c in S

        The inputs for the returned theano function are as follows:

            alpha_first_c: forward log likelihood (alpha) for t = t_0 and c in X
            beta_last_c: backward log likelihood (beta) for t = t_N and c in X
            log_emission_tc: log emission probabilities for t \in [t_0, ..., t_N] and c in X
            log_trans_tcc: log transition probabilities from `t` to `t+1` for t in [t_0, ..., t_N]
                and departure and destination states in X
            log_data_likelihood: log data likelihood of the unconstrained problem

        The output is a non-positive scalar value that signifies the desired probability in log space.

        Examples:

            If X = S (all hidden states), we expect logp = 0 (up to round-off error)

            if X = {a single hidden state}, then we expect the logp of a single path that takes on the same
                hidden-state for all positions [t_0, ..., t_N]

            In general, if X is a proper subset of S, we expect logp <= 0 (with logp = 0 iff the removed states
                are strictly forbidden by the prior and/or the transition matrix)

        Returns:
            a theano function
        """
        alpha_first_c = tt.vector('alpha_first_c')
        beta_last_c = tt.vector('beta_last_c')
        log_emission_tc = tt.matrix('log_emission_tc')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_data_likelihood = tt.scalar('log_data_likelihood')

        def update_alpha(c_log_emission_c: tt.vector,
                         c_log_trans_cc: tt.matrix,
                         p_alpha_c: tt.vector):
            return c_log_emission_c + pm.math.logsumexp(
                p_alpha_c.dimshuffle(0, 'x') + c_log_trans_cc, axis=0).dimshuffle(1)

        alpha_seg_iters, _ = th.scan(
            fn=update_alpha,
            sequences=[log_emission_tc, log_trans_tcc],
            outputs_info=[alpha_first_c])
        alpha_seg_end_c = alpha_seg_iters[-1, :]

        inputs = [alpha_first_c, beta_last_c, log_emission_tc, log_trans_tcc, log_data_likelihood]
        output = pm.math.logsumexp(alpha_seg_end_c + beta_last_c) - log_data_likelihood
        return th.function(inputs=inputs, outputs=output)

    # make a private static instance
    _constrained_path_logp_theano_func = _get_compiled_constrained_path_logp_theano_func.__func__()

    def get_log_constrained_posterior_prob(self,
                                           start_index: int, end_index: int,
                                           allowed_states: List[int]) -> float:
        """Calculates the constrained log posterior probability for contiguous set of sites in a Markov chain.
        At each site, only a subset of all states (as set by `allowed_states`) are allowed and the other states
        are strictly avoided.

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            allowed_states: the list of allowed states in the segment

        Returns:
            log constrained posterior probability (float)
        """
        assert start_index >= 0
        assert end_index < self.num_sites
        assert end_index >= start_index
        assert all(isinstance(item, int) and 0 <= item < self.num_states for item in allowed_states), \
            "The set of allowed states must be integers and in range [0, {0}]".format(self.num_states - 1)
        constrained_alpha_first_c = self.alpha_tc[start_index, allowed_states]
        constrained_beta_last_c = self.beta_tc[end_index, allowed_states]

        if end_index == start_index:  # single-site segment
            log_constrained_data_likelihood: float = logsumexp(constrained_alpha_first_c + constrained_beta_last_c)
            logp = log_constrained_data_likelihood - self.log_data_likelihood
        else:
            # calculate the required slices of the log emission and log transition representing
            # paths that only contain the allowed states
            constrained_log_emission_tc = \
                self.log_emission_tc[(start_index + 1):(end_index + 1), allowed_states]
            constrained_log_trans_tcc = \
                self.log_trans_tcc[start_index:end_index, allowed_states, :][:, :, allowed_states]
            logp = self._constrained_path_logp_theano_func(
                constrained_alpha_first_c, constrained_beta_last_c,
                constrained_log_emission_tc, constrained_log_trans_tcc, self.log_data_likelihood)

        return np.asscalar(logp)

    def get_segment_quality_some_called(self, start_index: int, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that one or more ("some") sites in a segment have
        the same hidden state ("call").

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert call_state in self.all_states_set

        all_other_states_logp = self.get_log_constrained_posterior_prob(
            start_index, end_index, self.leave_one_out_state_lists[call_state])

        return logp_to_phred(all_other_states_logp, complement=False)

    def get_segment_quality_all_called(self, start_index: int, end_index: int, call_state: int,
                                       quality_switch_threshold: float = 60.0) -> float:
        """Calculates the complementary phred-scaled posterior probability that "all" sites in a segment have
        the same hidden state ("call").

        Note:
            If all of the intervals in the segment overwhelmingly support the call state, the probability of
            deviations from the call state become very small, resulting in numerical instabilities.

            In such cases, we calculate the the quality assuming that the correlations in the posterior
            distribution are negligible, i.e.:

            log(1 - p(c_1 = c_2 = ... = call)) ~ log(1 - p(c_1 = call) x p(c_2 = call) x ...)
                                               = log(1 - (1 - p(c_1 != call)) x (1 - p(c_2 != call)) x ...)

            We calculate the latter expression using a robust numerical algorithm implemented in
            `gcnvkernel.utils.math.logsumexp_double_complement`.

            Since this calculation is relatively cheap, we always calculate the exact quality via the above
            scheme. If the uncorrelated phred-scale quality falls below `quality_switch_threshold`, we
            conclude that the segment is not high-quality and the correlated calculation is stable. Hence,
            we switch to the exact result.

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            call_state: segment call state index
            quality_switch_threshold: if the approximate quality is below this value, correlations will be
                taken into account.

        Returns:
            a phred-scaled probability
        """
        assert call_state in self.all_states_set

        if start_index == end_index:
            log_compl_prob = logsumexp(
                self.log_posterior_prob_tc[start_index, self.leave_one_out_state_lists[call_state]])
            return logp_to_phred(log_compl_prob, complement=False)
        else:
            # calculate the uncorrelated log complementary probability
            uncorrelated_log_compl_prob_array = np.asarray(
                [logsumexp(self.log_posterior_prob_tc[ti, self.leave_one_out_state_lists[call_state]])
                 for ti in range(start_index, end_index + 1)])
            uncorrelated_logp = logsumexp_double_complement(uncorrelated_log_compl_prob_array)
            uncorrelated_quality = logp_to_phred(uncorrelated_logp, complement=False)
            exact_quality = None
            if uncorrelated_quality < quality_switch_threshold:
                # fallback to exact quality calculation
                all_called_state_logp = self.get_log_constrained_posterior_prob(
                    start_index, end_index, [call_state])
                exact_quality = logp_to_phred(all_called_state_logp, complement=True)
            if exact_quality is not None and not np.isnan(exact_quality):
                return exact_quality
            else:
                return uncorrelated_quality

    def get_segment_quality_start(self, start_index: int, call_state: int) -> float:
        """Calculates the complement of phred-scaled posterior probability that a site marks the start of a segment.
        This is done by directly summing the probability of the following complementary paths in the log space:

            ...      [start_index-1]     [start_index]                          ...
                     (any state)     =>  (any state except for call_state)
                     call_state      =>  call_state

        Args:
            start_index: left breakpoint index of a segment
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert 0 <= start_index < self.num_sites

        all_other_states_list = self.leave_one_out_state_lists[call_state]
        if start_index == 0:
            logp = logsumexp(self.log_posterior_prob_tc[0, all_other_states_list])
        else:
            complementary_paths_unnorm_logp = [(self.alpha_tc[start_index - 1, prev_state] +
                                                self.log_trans_tcc[start_index - 1, prev_state, start_state] +
                                                self.log_emission_tc[start_index, start_state] +
                                                self.beta_tc[start_index, start_state])
                                               for prev_state in self.all_states_list
                                               for start_state in all_other_states_list]
            complementary_paths_unnorm_logp.append((self.alpha_tc[start_index - 1, call_state] +
                                                    self.log_trans_tcc[start_index - 1, call_state, call_state] +
                                                    self.log_emission_tc[start_index, call_state] +
                                                    self.beta_tc[start_index, call_state]))
            logp = logsumexp(np.asarray(complementary_paths_unnorm_logp)) - self.log_data_likelihood

        return logp_to_phred(logp)

    def get_segment_quality_end(self, end_index: int, call_state: int) -> float:
        """Calculates the complement of phred-scaled posterior probability that a site marks the end of a segment.
        This is done by directly summing the probability of the following complementary paths in the log space:

            ...      [end_index]                           [end_index+1]          ...
                     call_state                        =>  call_state
                     (any state except for call_state) =>  (any state)

        Args:
            end_index: right breakpoint index of a segment
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert 0 <= end_index < self.num_sites

        all_other_states_list = self.leave_one_out_state_lists[call_state]

        if end_index == self.num_sites - 1:
            logp = logsumexp(self.log_posterior_prob_tc[self.num_sites - 1, all_other_states_list])
        else:
            complementary_paths_unnorm_logp = [(self.alpha_tc[end_index, end_state] +
                                                self.log_trans_tcc[end_index, end_state, next_state] +
                                                self.log_emission_tc[end_index + 1, next_state] +
                                                self.beta_tc[end_index + 1, end_state])
                                               for end_state in all_other_states_list
                                               for next_state in self.all_states_list]
            complementary_paths_unnorm_logp.append((self.alpha_tc[end_index, call_state] +
                                                    self.log_trans_tcc[end_index, call_state, call_state] +
                                                    self.log_emission_tc[end_index + 1, call_state] +
                                                    self.beta_tc[end_index + 1, call_state]))
            logp = logsumexp(np.asarray(complementary_paths_unnorm_logp)) - self.log_data_likelihood

        return logp_to_phred(logp)