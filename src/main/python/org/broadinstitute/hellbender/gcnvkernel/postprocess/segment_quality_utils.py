from typing import Dict, List
import numpy as np
import theano as th
import pymc3 as pm
import theano.tensor as tt
from scipy.misc import logsumexp


class HMMSegmentationQualityCalculator:
    """Calculates quality metrics for hidden state segments for a given HMM.

    Note:
        The initializer requires the emission and transition probabilities, as well as the forward
        and backward tables and the log posterior probability.
    """

    # 10 / ln(10)
    INV_LN_10_TIMES_10 = 4.342944819032518

    # ln(1/2)
    LN_HALF = -0.6931471805599453

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
        self.all_states_set_except_for: Dict[int, List[int]] = {
            except_state: [state for state in range(self.num_states) if state != except_state]
            for except_state in range(self.num_states)}

    @staticmethod
    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_constrained_path_log_prob_theano_func():
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
    _constrained_path_log_prob_theano_func = _get_compiled_constrained_path_log_prob_theano_func.__func__()

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
            log_prob = log_constrained_data_likelihood - self.log_data_likelihood
        else:
            # calculate the required slices of the log emission and log transition representing
            # paths that only contain the allowed states
            constrained_log_emission_tc = \
                self.log_emission_tc[(start_index + 1):(end_index + 1), allowed_states]
            constrained_log_trans_tcc = \
                self.log_trans_tcc[start_index:end_index, allowed_states, :][:, :, allowed_states]
            log_prob = self._constrained_path_log_prob_theano_func(
                constrained_alpha_first_c, constrained_beta_last_c,
                constrained_log_emission_tc, constrained_log_trans_tcc, self.log_data_likelihood)

        return np.asscalar(log_prob)

    def get_segment_some_quality(self, start_index: int, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that one or more ("some") sites in a segment have
        the same hidden state ("call") normalized by the length of the segment.

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert call_state in self.all_states_set

        all_other_states_log_prob = self.get_log_constrained_posterior_prob(
            start_index, end_index, self.all_states_set_except_for[call_state])

        return self.log_prob_to_phred(all_other_states_log_prob, complement=False) / (end_index - start_index + 1)

    def get_segment_exact_quality(self, start_index: int, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that "all" sites in a segment have the same
        hidden state ("call").

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert call_state in self.all_states_set

        if start_index == end_index:
            return self.log_prob_to_phred(self.log_posterior_prob_tc[start_index, call_state], complement=True)
        else:
            all_called_state_log_prob = self.get_log_constrained_posterior_prob(start_index, end_index, [call_state])
            return self.log_prob_to_phred(all_called_state_log_prob, complement=True)

    def get_segment_start_quality(self, start_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that a site marks the start of a segment.

        Args:
            start_index: left breakpoint index of a segment
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert 0 <= start_index < self.num_sites
        if start_index == 0:
            log_prob = self.log_posterior_prob_tc[0, call_state]
        else:
            # calculate the probability of all paths that start from other states and end up with the called state
            all_other_states_list = self.all_states_set_except_for[call_state]
            prev_alpha_c = self.alpha_tc[start_index - 1, all_other_states_list]
            current_beta = self.beta_tc[start_index, call_state]
            current_log_emission = self.log_emission_tc[start_index, call_state]
            log_trans_c = self.log_trans_tcc[start_index - 1, all_other_states_list, call_state]
            log_breakpoint_likelihood = logsumexp(prev_alpha_c + log_trans_c + current_log_emission) + current_beta
            log_prob = log_breakpoint_likelihood - self.log_data_likelihood

        return self.log_prob_to_phred(log_prob, complement=True)

    def get_segment_start_quality_direct(self, start_index: int, call_state: int) -> float:
        """Calculates the complement of phred-scaled posterior probability that a site does _not_ mark the
        start of a segment. This is done by directly summing the probability of the following paths:

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

        all_other_states_list = self.all_states_set_except_for[call_state]
        if start_index == 0:
            log_prob = logsumexp(self.log_posterior_prob_tc[0, all_other_states_list])
        else:
            complementary_paths_unnorm_log_prob = [(self.alpha_tc[start_index - 1, prev_state] +
                                                    self.log_trans_tcc[start_index - 1, prev_state, start_state] +
                                                    self.log_emission_tc[start_index, start_state] +
                                                    self.beta_tc[start_index, start_state])
                                                   for prev_state in self.all_states_list
                                                   for start_state in all_other_states_list]
            complementary_paths_unnorm_log_prob.append((self.alpha_tc[start_index - 1, call_state] +
                                                        self.log_trans_tcc[start_index - 1, call_state, call_state] +
                                                        self.log_emission_tc[start_index, call_state] +
                                                        self.beta_tc[start_index, call_state]))
            log_prob = logsumexp(np.asarray(complementary_paths_unnorm_log_prob)) - self.log_data_likelihood

        return self.log_prob_to_phred(log_prob)

    def get_segment_end_quality(self, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that a site marks the end of a segment.

        Args:
            end_index: right breakpoint index of a segment
            call_state: segment call state index

        Returns:

        """
        assert 0 <= end_index < self.num_sites

        if end_index == self.num_sites - 1:
            log_prob = self.log_posterior_prob_tc[self.num_sites - 1, call_state]
        else:
            # calculate the probability of all paths that start from call state and end up with other states
            all_other_states_list = self.all_states_set_except_for[call_state]
            current_alpha = self.alpha_tc[end_index, call_state]
            next_beta_c = self.beta_tc[end_index + 1, all_other_states_list]
            next_log_emission_c = self.log_emission_tc[end_index + 1, all_other_states_list]
            log_trans_c = self.log_trans_tcc[end_index, call_state, all_other_states_list]
            log_breakpoint_likelihood = logsumexp(current_alpha + log_trans_c + next_log_emission_c + next_beta_c)
            log_prob = log_breakpoint_likelihood - self.log_data_likelihood

        return self.log_prob_to_phred(log_prob, complement=True)

    def get_segment_end_quality_direct(self, end_index: int, call_state: int) -> float:
        """Calculates the complement of phred-scaled posterior probability that a site does _not_ mark the
        end of a segment. This is done by directly summing the probability of the following paths:

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

        all_other_states_list = self.all_states_set_except_for[call_state]

        if end_index == self.num_sites - 1:
            log_prob = logsumexp(self.log_posterior_prob_tc[self.num_sites - 1, all_other_states_list])
        else:
            complementary_paths_unnorm_log_prob = [(self.alpha_tc[end_index, end_state] +
                                                    self.log_trans_tcc[end_index, end_state, next_state] +
                                                    self.log_emission_tc[end_index + 1, next_state] +
                                                    self.beta_tc[end_index + 1, end_state])
                                                   for end_state in all_other_states_list
                                                   for next_state in self.all_states_list]
            complementary_paths_unnorm_log_prob.append((self.alpha_tc[end_index, call_state] +
                                                        self.log_trans_tcc[end_index, call_state, call_state] +
                                                        self.log_emission_tc[end_index + 1, call_state] +
                                                        self.beta_tc[end_index + 1, call_state]))
            log_prob = logsumexp(np.asarray(complementary_paths_unnorm_log_prob)) - self.log_data_likelihood

        return self.log_prob_to_phred(log_prob)

    @staticmethod
    def log_prob_to_phred(log_prob: float, complement: bool = False) -> float:
        """Converts probabilities in natural log scale to phred scale.

        Args:
            log_prob: a probability in the natural log scale
            complement: invert the probability

        Returns:
            phred-scaled probability
        """
        final_log_prob = log_prob if not complement else HMMSegmentationQualityCalculator.log_prob_complement(log_prob)
        return -final_log_prob * HMMSegmentationQualityCalculator.INV_LN_10_TIMES_10

    @staticmethod
    def log_prob_complement(log_prob: float) -> float:
        """Calculates the complement of a probability in the natural log scale.

        Args:
            log_prob: a probability in the natural log scale

        Returns:
            complement of the the probability in the natural log scale
        """
        log_prob_zero_capped = min(0., log_prob)
        if log_prob_zero_capped >= HMMSegmentationQualityCalculator.LN_HALF:
            return np.log(-np.expm1(log_prob_zero_capped))
        else:
            return np.log1p(-np.exp(log_prob_zero_capped))
