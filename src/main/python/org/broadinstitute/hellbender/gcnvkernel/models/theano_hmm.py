from typing import Optional, List, Tuple, Union

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt

from . import commons
from .. import types


class TheanoForwardBackward:
    """Implementation of the forward-backward algorithm in theano."""
    def __init__(self,
                 log_posterior_probs_output_tc: Optional[types.TensorSharedVariable] = None,
                 resolve_nans: bool = False,
                 do_thermalization: bool = False,
                 do_admixing: bool = False,
                 include_update_size_output: bool = False,
                 include_alpha_beta_output: bool = False):
        """Initializes the forward-backward algorithm by compiling a theano function according to the
        boolean flags.

        Args:
            log_posterior_probs_output_tc: if not None, the new log posterior will be written to this shared tensor;
                otherwise, it will be returned as an ndarray
            resolve_nans: if True, expression such as inf - inf resulting in NaNs will be properly handled
            do_thermalization: if True, performs thermalization of HMM parameters
            do_admixing: if True, perform admixing of old and new hidden-state posterior probabilities
            include_update_size_output: if True, include update size in the returned values
            include_alpha_beta_output: include forward and backward tables in the return values
        """
        self.resolve_nans = resolve_nans
        self.log_posterior_probs_output_tc = log_posterior_probs_output_tc
        self.do_thermalization = do_thermalization
        self.do_admixing = do_admixing
        self.include_update_size_output = include_update_size_output
        self.include_alpha_beta_output = include_alpha_beta_output
        self._forward_backward_theano_func = self._get_compiled_forward_backward_theano_func()

    def perform_forward_backward(self,
                                 log_prior_c: np.ndarray,
                                 log_trans_tcc: np.ndarray,
                                 log_emission_tc: np.ndarray,
                                 prev_log_posterior_tc: Optional[np.ndarray] = None,
                                 admixing_rate: Optional[float] = None,
                                 temperature: Optional[float] = None) -> 'ForwardBackwardResult':
        """Runs the forward-backward algorithm.

        Notes:
            The inputs args must be compatible with the compiled theano function according to the
            class initializer flags.

        Args:
            log_prior_c: prior probability vector for the first node
            log_trans_tcc: transition probability matrices for each directed vertex
            log_emission_tc: emission probability vector for each node
            prev_log_posterior_tc: (optional) previous estimate of the log posterior
                (used if `self.do_admixing` is True)
            admixing_rate: (optional) a float in range [0, 1] denoting the amount of the new posterior probabilities
                to admix with the old posterior probabilities (higher = more of the new posterior)
            temperature: (optional) temperature (used if `self.do_thermalization` is True)

        Returns:
            an instance of `ForwardBackwardResult`
        """
        if self.do_admixing:
            assert prev_log_posterior_tc is not None,\
                "Posterior admixing is enabled but `prev_log_posterior_tc` is not specified."
            assert admixing_rate is not None,\
                "Posterior admixing is enabled but `admixing_rate` is not specified."
            assert 0.0 <= admixing_rate <= 1.0,\
                "Posterior admixing rate must be in the range [0, 1]."

        if self.do_thermalization:
            assert temperature is not None,\
                "Thermalization is enabled but `temperature` is not specified."
            assert temperature > 0,\
                "Temperature must be non-negative."

        if self.include_update_size_output:
            assert prev_log_posterior_tc is not None,\
                "Update size output is enabled but `prev_log_posterior_tc` is not specified."

        return self._decompose_theano_forward_backward_outputs(
            self._forward_backward_theano_func(*self._compose_theano_forward_backward_inputs(
                log_prior_c, log_trans_tcc, log_emission_tc,
                prev_log_posterior_tc, admixing_rate, temperature)))

    def _compose_theano_forward_backward_inputs(self,
                                                log_prior_c: np.ndarray,
                                                log_trans_tcc: np.ndarray,
                                                log_emission_tc: np.ndarray,
                                                prev_log_posterior_tc: Optional[np.ndarray],
                                                admixing_rate: Optional[float],
                                                temperature: Optional[float]) -> Tuple:
        inputs = (log_prior_c, log_trans_tcc, log_emission_tc)
        if self.do_admixing or self.include_update_size_output:
            inputs += (prev_log_posterior_tc,)
        if self.do_admixing:
            inputs += (admixing_rate,)
        if self.do_thermalization:
            inputs += (temperature,)
        return inputs

    def _decompose_theano_forward_backward_outputs(self, outputs: List[Union[np.ndarray, float]])\
            -> 'ForwardBackwardResult':
        result = ForwardBackwardResult()
        arg_idx = 0
        if self.log_posterior_probs_output_tc is None:
            result.log_posterior_probs_tc = outputs[arg_idx]
            arg_idx += 1
        result.log_data_likelihood = outputs[arg_idx]
        arg_idx += 1
        if self.include_update_size_output:
            result.update_norm_t = outputs[arg_idx]
            arg_idx += 1
        if self.include_alpha_beta_output:
            result.alpha_tc = outputs[arg_idx]
            arg_idx += 1
            result.beta_tc = outputs[arg_idx]
            arg_idx += 1
        return result

    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_forward_backward_theano_func(self) -> th.compile.function_module.Function:
        """Returns a compiled theano function that computes the posterior probabilities of hidden states using
        the forward-backward algorithm.

        Note:
            The input arguments and the output of the compiled theano function is determined by the initializer flags
            as follows:

            There are 3 basic input arguments:

                * log_prior_c (float vector),
                * log_trans_tcc (float tensor3),
                * log_emission_tc (float matrix)

            The rest of the input arguments must be concatenated to the basic arguments, in order, and as follows:

                If either `self.do_admixing == True` or `self.include_update_size_output == True`, then concatenate:
                    * prev_log_posterior_tc (float matrix)

                If `self.do_admixing` is True, then concatenate:
                    * admixing_rate (float scalar)

                if `self.do_thermalization` is True, then concatenate:
                    * temperature (float scalar)

            The outputs list is built, in order, and as follows:

                * If `log_posterior_probs_output_tc` (a shared tensor) is given to the class initializer, the computed
                  hidden state log posterior probabilities will be directly written to `log_posterior_probs_output_tc`.
                  Otherwise, it will be returned as the first entry of the outputs list.

                * The next entry is `log_data_likelihood` (float scalar).

                * If `self.include_update_size_output == True`, the next entry is `update_norm_t` (float vector).

                * If `self.include_alpha_beta_output == True`, the next two entries will be `alpha_tc` (float vector)
                  and `beta_tc` (float vector).

        Returns:
            A compiled theano function
        """
        # basic inputs
        log_prior_c = tt.vector('log_prior_c')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_emission_tc = tt.matrix('log_emission_tc')

        # optional inputs
        prev_log_posterior_tc = tt.matrix('prev_log_posterior_tc')
        admixing_rate = tt.scalar('admixing_rate')
        temperature = tt.scalar('temperature')

        if self.do_thermalization:
            processed_log_prior_c, processed_log_trans_tcc, processed_log_emission_tc =\
                self.get_symbolic_thermal_hmm_params(log_prior_c, log_trans_tcc, log_emission_tc, temperature)
        else:
            processed_log_prior_c, processed_log_trans_tcc, processed_log_emission_tc =\
                log_prior_c, log_trans_tcc, log_emission_tc

        # get symbolic forward-backward
        new_log_posterior_tc, log_data_likelihood_t, alpha_tc, beta_tc = self.get_symbolic_log_posterior(
            processed_log_prior_c, processed_log_trans_tcc, processed_log_emission_tc, self.resolve_nans)

        if self.do_admixing:
            processed_log_posterior_tc = commons.safe_logaddexp(
                new_log_posterior_tc + tt.log(admixing_rate),
                prev_log_posterior_tc + tt.log(1.0 - admixing_rate))
        else:
            processed_log_posterior_tc = new_log_posterior_tc

        log_data_likelihood = log_data_likelihood_t[-1]  # in theory, they are all the same

        # build the updates list
        if self.log_posterior_probs_output_tc is not None:
            updates = [(self.log_posterior_probs_output_tc, processed_log_posterior_tc)]
        else:
            updates = None

        # build the inputs list
        inputs = [log_prior_c, log_trans_tcc, log_emission_tc]
        if self.do_admixing or self.include_update_size_output:
            inputs += [prev_log_posterior_tc]
        if self.do_admixing:
            inputs += [admixing_rate]
        if self.do_thermalization:
            inputs += [temperature]

        # build the outputs list
        outputs = []
        if self.log_posterior_probs_output_tc is None:
            outputs += [processed_log_posterior_tc]

        outputs += [log_data_likelihood]

        if self.include_update_size_output:
            update_norm_t = commons.get_jensen_shannon_divergence(processed_log_posterior_tc, prev_log_posterior_tc)
            outputs += [update_norm_t]

        if self.include_alpha_beta_output:
            outputs += [alpha_tc, beta_tc]

        return th.function(inputs=inputs, outputs=outputs, updates=updates)

    @staticmethod
    def get_symbolic_log_posterior(log_prior_c: types.TheanoVector,
                                   log_trans_tcc: types.TheanoTensor3,
                                   log_emission_tc: types.TheanoMatrix,
                                   resolve_nans: bool):
        """Generates symbolic tensors representing hidden-state log posterior, log data likelihood,
        forward table (alpha), and backward table (beta).
        
        Returns:
            log_posterior_probs, log_data_likelihood
        """
        num_states = log_prior_c.shape[0]

        def calculate_next_alpha(c_log_trans_ab: types.TheanoMatrix,
                                 c_log_emission_b: types.TheanoVector,
                                 p_alpha_a: types.TheanoVector):
            """Calculates the next entry on the forward table, alpha_{t}, from alpha_{t-1}.

            Args:
                c_log_trans_ab: a 2d tensor with rows and columns corresponding to log transition probability
                    from the previous state at position t-1 and to the next state at position t, respectively
                c_log_emission_b: a 1d tensor representing the emission probability to each state at position t
                p_alpha_a: a 1d tensor representing alpha_{t-1}

            Returns:
                symbolic 1d tensor of alpha_{t}
            """
            mu_ba = tt.tile(p_alpha_a, (num_states, 1)) + c_log_trans_ab.T
            n_alpha_b = c_log_emission_b + pm.math.logsumexp(mu_ba, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(n_alpha_b), -np.inf, n_alpha_b)
            else:
                return n_alpha_b

        def calculate_prev_beta(n_log_trans_ab: types.TheanoMatrix,
                                n_log_emission_b: types.TheanoVector,
                                n_beta_b: types.TheanoVector):
            """Calculates the previous entry on the backward table, beta_{t-1}, from beta_{t}.

            Args:
                n_log_trans_ab: a 2d tensor with rows and columns corresponding to log transition probability
                    from the previous state at position t-1 and to the next state at position t, respectively
                n_log_emission_b: a 1d tensor representing the emission probability to each state at position t
                n_beta_b: a 1d tensor representing beta_{t}

            Returns:
                symbolic 1d tensor of beta_{t-1}
            """
            nu_ab = tt.tile(n_beta_b + n_log_emission_b, (num_states, 1)) + n_log_trans_ab
            p_beta_a = pm.math.logsumexp(nu_ab, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(p_beta_a), -np.inf, p_beta_a)
            else:
                return p_beta_a

        # first entry of the forward table
        first_alpha_c = log_prior_c + log_emission_tc[0, :]

        # the rest of the forward table
        rest_alpha_tc, alpha_updates = th.scan(
            fn=calculate_next_alpha,
            sequences=[log_trans_tcc, log_emission_tc[1:, :]],
            outputs_info=[first_alpha_c])

        # concatenate with the first alpha
        alpha_tc = tt.concatenate((first_alpha_c.dimshuffle('x', 0), rest_alpha_tc))

        # last entry of the backward table (zero for all states)
        last_beta_c = tt.zeros_like(log_prior_c)

        # the rest of the backward table
        rest_beta_tc, beta_updates = th.scan(
            fn=calculate_prev_beta,
            sequences=[log_trans_tcc, log_emission_tc[1:, :]],
            go_backwards=True,
            outputs_info=[last_beta_c])

        # concatenate with the last beta and reverse
        beta_tc = tt.concatenate((last_beta_c.dimshuffle('x', 0), rest_beta_tc))[::-1, :]

        # calculate normalized log posterior
        log_unnormalized_posterior_tc = alpha_tc + beta_tc
        log_data_likelihood_t = pm.math.logsumexp(log_unnormalized_posterior_tc, axis=1)
        log_posterior_probs_tc = log_unnormalized_posterior_tc - log_data_likelihood_t

        return log_posterior_probs_tc, log_data_likelihood_t.dimshuffle(0), alpha_tc, beta_tc

    @staticmethod
    def get_symbolic_thermal_hmm_params(log_prior_c: types.TheanoVector,
                                        log_trans_tcc: types.TheanoTensor3,
                                        log_emission_tc: types.TheanoMatrix,
                                        temperature: tt.scalar):
        inv_temperature = tt.inv(temperature)

        thermal_log_prior_c = inv_temperature * log_prior_c
        thermal_log_prior_c -= pm.math.logsumexp(thermal_log_prior_c)
        thermal_log_trans_tcc = inv_temperature * log_trans_tcc
        thermal_log_trans_tcc -= pm.math.logsumexp(thermal_log_trans_tcc, axis=-1)
        thermal_log_emission_tc = inv_temperature * log_emission_tc

        return thermal_log_prior_c, thermal_log_trans_tcc, thermal_log_emission_tc


class ForwardBackwardResult:
    """Stores the output of forward-backward algorithm."""
    def __init__(self,
                 log_posterior_probs_tc: Optional[np.ndarray] = None,
                 log_data_likelihood: Optional[float] = None,
                 alpha_tc: Optional[np.ndarray] = None,
                 beta_tc: Optional[np.ndarray] = None,
                 update_norm_t: Optional[np.ndarray] = None):
        self.log_posterior_probs_tc = log_posterior_probs_tc
        self.log_data_likelihood = log_data_likelihood
        self.alpha_tc = alpha_tc
        self.beta_tc = beta_tc
        self.update_norm_t = update_norm_t


class TheanoViterbi:
    """Implementation of the Viterbi algorithm in theano."""
    def __init__(self):
        self._viterbi_theano_func = self._get_compiled_viterbi_theano_func()

    def get_viterbi_path(self,
                         log_prior_c: np.ndarray,
                         log_trans_tcc: np.ndarray,
                         log_emission_tc: np.ndarray) -> List[int]:
        return self._viterbi_theano_func(log_prior_c, log_trans_tcc, log_emission_tc).tolist()

    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_viterbi_theano_func(self) -> th.compile.function_module.Function:
        """Returns a theano function that calculates the Viterbi path."""
        log_prior_c = tt.vector('log_prior_c')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_emission_tc = tt.matrix('log_emission_tc')

        return th.function(inputs=[log_prior_c, log_trans_tcc, log_emission_tc],
                           outputs=self._get_symbolic_viterbi_path(log_prior_c, log_trans_tcc, log_emission_tc))

    @staticmethod
    def _get_symbolic_viterbi_path(log_prior_c: types.TheanoVector,
                                   log_trans_tcc: types.TheanoTensor3,
                                   log_emission_tc: types.TheanoMatrix):
        """Generates a symbolic 1d integer tensor representing the most-likely chain of hidden states
        (Viterbi algorithm).

        Notes:
            In the following, `omega_tc` refers to the log data likelihood at position `t` for a max sum-product
            path ending with hidden state `c`. Also, `psi_tc` refers to the max sum-product backtracking
            table, i.e. `psi_tc` represents the best hidden state at position t-1 for a max sum-product path
            ending with hidden state `c` at position t.

        Returns:
            symbolic 1d integer tensor representing the most-likely chain of hidden states
        """

        def calculate_next_omega_psi(p_log_trans_ab, c_log_emission_b, p_omega_a):
            """Extends the max sum-product path by one position and calculates the log data likelihood of such paths
            for each final hidden state (`n_omega_b`), as well as the most-likely terminal state of the path at the
            previous position, assuming that it lands on state `b` after extension (`psi_b`).

            Args:
                p_log_trans_ab: log transition matrix from `a` to `b`
                c_log_emission_b: log emission probabilities at the current position for state `b`
                p_omega_a: previous log data likelihood of the max sum-product path ending with hidden state `a`

            Returns:
                next omega, next psi
            """
            tau_ab = p_log_trans_ab + p_omega_a.dimshuffle(0, 'x')
            max_tau_b, psi_b = tt.max_and_argmax(tau_ab, axis=0)
            n_omega_b = c_log_emission_b + max_tau_b
            return n_omega_b, psi_b

        def calculate_previous_best_state(c_psi_c, c_best_state):
            """Backtracks the max sum-product table to the previous position.

            Args:
                c_psi_c: backtracking table column at the current position
                c_best_state: best state at the current position

            Returns:
                most-likely state at the previous position
            """
            return c_psi_c[c_best_state]

        # log data likelihood for each hidden state at the first position
        omega_first_a = log_emission_tc[0, :] + log_prior_c

        # calculate the log data likelihood of the partial max sum-product paths (omega)
        # and the backtracking table (psi)
        omega_psi_list, _ = th.scan(
            fn=calculate_next_omega_psi,
            sequences=[log_trans_tcc, log_emission_tc[1:, :]],
            outputs_info=[omega_first_a, None])
        omega_tc = omega_psi_list[0]
        psi_tc = omega_psi_list[1]

        # the best terminal state
        last_best_state = tt.argmax(omega_tc[-1, :])

        # backtrack to obtain the previous states of the max sum-product path
        rest_best_states_t, _ = th.scan(
            fn=calculate_previous_best_state,
            sequences=[psi_tc],
            outputs_info=[last_best_state],
            go_backwards=True)

        # concatenate with the terminal state
        viterbi_path_t = tt.concatenate([tt.stack(last_best_state), rest_best_states_t])[::-1]

        return viterbi_path_t
