import numpy as np
import theano as th
import theano.tensor as tt
import pymc3 as pm
from typing import Optional
from .. import types
from . import commons


class TheanoForwardBackward:
    """Implementation of the forward-backward algorithm using `theano.scan`."""
    def __init__(self,
                 log_posterior_output: Optional[types.TensorSharedVariable],
                 admixing_rate: float,
                 include_alpha_beta_output: bool = False,
                 resolve_nans: bool = False):
        """Initializer.

        Args:
            log_posterior_output: if not None, the new log posterior will be written to this shared tensor;
                otherwise, it will be returned by `TheanoForwardBackward.perform_forward_backward`
            admixing_rate: a float in range [0, 1] denoting the amount of the new posterior to admix with the
                old posterior (higher = more of the new posterior)
            include_alpha_beta_output: include forward and backward tables in the return values
                of `TheanoForwardBackward.perform_forward_backward`
            resolve_nans: if True, expression such as inf - inf resulting in NaNs will be properly handled
        """
        self.admixing_rate = admixing_rate
        self.include_alpha_beta_output = include_alpha_beta_output
        self.resolve_nans = resolve_nans
        assert 0.0 < admixing_rate <= 1.0, "Admixing rate must be in range (0, 1]"
        self.log_posterior_output = log_posterior_output
        self._forward_backward_theano_func = self._get_compiled_forward_backward_theano_func()

    def perform_forward_backward(self,
                                 num_states: int,
                                 temperature: float,
                                 log_prior_c: np.ndarray,
                                 log_trans_tcc: np.ndarray,
                                 log_emission_tc: np.ndarray,
                                 prev_log_posterior_tc: np.ndarray):
        """Runs the forward-backward algorithm.

        Args:
            num_states: number of states for each node of the Markov chain
            temperature: a scale factor of the chain entropy
            log_prior_c: prior probability vector for the first node
            log_trans_tcc: transition probability matrices for each directed vertex
            log_emission_tc: emission probability vector for each node
            prev_log_posterior_tc: previous estimate of the log posterior (used for admixing)

        Returns:
            see `TheanoForwardBackward._get_compiled_forward_backward_theano_func`
        """
        return self._forward_backward_theano_func(
            num_states, temperature, log_prior_c, log_trans_tcc, log_emission_tc, prev_log_posterior_tc)

    # todo update docstring
    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_forward_backward_theano_func(self):
        """Returns a compiled theano function that perform forward-backward and either updates log posterior
        probabilities or returns it.

        Note:
            The returned theano function takes 6 inputs:

                num_states (integer scalar),
                temperature (float scalar),
                log_prior_c (float vector),
                og_trans_tcc (float tensor3),
                log_emission_tc (float matrix)
                prev_log_posterior_tc (float matrix)


            If a `log_posterior_output` shared tensor is given to the class initializer,
            the return tuple will be:

                update_norm_t, log_data_likelihood,
                (+ alpha_tc, beta_tc if self.include_alpha_beta_output == True)

            and the posterior will be directly written to `self.log_posterior_output`. Otherwise,
            return tuple will be:

                admixed_log_posterior_tc, update_norm_t, log_data_likelihood,
                (+ alpha_tc, beta_tc if self.include_alpha_beta_output == True)

        Returns:
            A compiled theano function
        """
        num_states = tt.iscalar('num_states')
        temperature = tt.scalar('temperature')
        log_prior_c = tt.vector('log_prior_c')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_emission_tc = tt.matrix('log_emission_tc')
        prev_log_posterior_tc = tt.matrix('prev_log_posterior_tc')

        new_log_posterior_tc, log_data_likelihood_t, alpha_tc, beta_tc = self._get_symbolic_log_posterior(
            num_states, temperature, log_prior_c, log_trans_tcc, log_emission_tc, self.resolve_nans)

        admixed_log_posterior_tc = commons.safe_logaddexp(
            new_log_posterior_tc + np.log(self.admixing_rate),
            prev_log_posterior_tc + np.log(1.0 - self.admixing_rate))

        log_data_likelihood = log_data_likelihood_t[-1]  # in theory, they are all the same
        update_norm_t = commons.get_jensen_shannon_divergence(admixed_log_posterior_tc, prev_log_posterior_tc)

        ext_output = [alpha_tc, beta_tc] if self.include_alpha_beta_output else []
        inputs = [num_states, temperature, log_prior_c, log_trans_tcc, log_emission_tc, prev_log_posterior_tc]
        if self.log_posterior_output is not None:
            return th.function(inputs=inputs,
                               outputs=[update_norm_t, log_data_likelihood] + ext_output,
                               updates=[(self.log_posterior_output, admixed_log_posterior_tc)])
        else:
            return th.function(inputs=inputs,
                               outputs=[admixed_log_posterior_tc, update_norm_t, log_data_likelihood] + ext_output)

    @staticmethod
    def _get_symbolic_log_posterior(num_states: tt.iscalar,
                                    temperature: tt.scalar,
                                    log_prior_c: types.TheanoVector,
                                    log_trans_tcc: types.TheanoTensor3,
                                    log_emission_tc: types.TheanoMatrix,
                                    resolve_nans: bool):
        """Generates symbolic tensors for log posterior and log data likelihood.
        
        Returns:
            log_posterior_probs, log_data_likelihood
        """

        def calculate_next_alpha(c_log_trans_mat: types.TheanoMatrix,
                                 c_log_emission_vec: types.TheanoVector,
                                 p_alpha_vec: types.TheanoVector):
            """Calculates the next entry on the forward table, alpha_{t}, from alpha_{t-1}.

            Args:
                c_log_trans_mat: a 2d tensor with rows and columns corresponding to log transition probability
                    from the previous state at position t-1 and to the next state at position t, respectively
                c_log_emission_vec: a 1d tensor representing the emission probability to each state at position t
                p_alpha_vec: a 1d tensor representing alpha_{t-1}

            Returns:
                symbolic 1d tensor of alpha_{t}
            """
            mu = tt.tile(p_alpha_vec, (num_states, 1)) + c_log_trans_mat.T
            n_alpha_vec = c_log_emission_vec + pm.math.logsumexp(mu, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(n_alpha_vec), -np.inf, n_alpha_vec)
            else:
                return n_alpha_vec

        def calculate_prev_beta(n_log_trans_mat: types.TheanoMatrix,
                                n_log_emission_vec: types.TheanoVector,
                                n_beta_vec: types.TheanoVector):
            """Calculates the previous entry on the backward table, beta_{t-1}, from beta_{t}.

            Args:
                n_log_trans_mat: a 2d tensor with rows and columns corresponding to log transition probability
                    from the previous state at position t-1 and to the next state at position t, respectively
                n_log_emission_vec: a 1d tensor representing the emission probability to each state at position t
                n_beta_vec: a 1d tensor representing beta_{t}

            Returns:
                symbolic 1d tensor of beta_{t-1}
            """
            nu = tt.tile(n_beta_vec + n_log_emission_vec, (num_states, 1)) + n_log_trans_mat
            p_beta_vec = pm.math.logsumexp(nu, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(p_beta_vec), -np.inf, p_beta_vec)
            else:
                return p_beta_vec

        # calculate thermal equivalent of various quantities
        # todo are normalizations (i.e. logsumexp subtraction) necessary?
        inv_temperature = tt.inv(temperature)
        thermal_log_prior_c = inv_temperature * log_prior_c
        thermal_log_prior_c -= pm.math.logsumexp(thermal_log_prior_c)
        thermal_log_trans_tcc = inv_temperature * log_trans_tcc
        thermal_log_trans_tcc -= pm.math.logsumexp(thermal_log_trans_tcc, axis=-1)
        thermal_log_emission_tc = inv_temperature * log_emission_tc

        # first entry of the forward table
        alpha_first = thermal_log_prior_c + thermal_log_emission_tc[0, :]

        # the rest of the forward table
        alpha_seq, alpha_updates = th.scan(
            fn=calculate_next_alpha,
            sequences=[thermal_log_trans_tcc, thermal_log_emission_tc[1:, :]],
            outputs_info=[alpha_first])

        # concatenate with the first alpha
        alpha_t = tt.concatenate((alpha_first.dimshuffle('x', 0), alpha_seq))

        # last entry of the backward table (zero for all states)
        beta_last = tt.zeros_like(log_prior_c)

        # the rest of the backward table
        beta_seq, beta_updates = th.scan(
            fn=calculate_prev_beta,
            sequences=[thermal_log_trans_tcc, thermal_log_emission_tc[1:, :]],
            go_backwards=True,
            outputs_info=[beta_last])

        # concatenate with the last beta and reverse
        beta_t = tt.concatenate((beta_last.dimshuffle('x', 0), beta_seq))[::-1, :]

        # calculate normalized log posterior
        log_unnormalized_posterior_t = alpha_t + beta_t
        log_data_likelihood_t = pm.math.logsumexp(log_unnormalized_posterior_t, axis=1)
        log_posterior_t = log_unnormalized_posterior_t - log_data_likelihood_t

        return log_posterior_t, log_data_likelihood_t.dimshuffle(0), alpha_t, beta_t
