import argparse
import inspect
import io
import logging
import time
from abc import abstractmethod
from typing import List, Callable, Optional, Set, Tuple, Any, Dict

import numpy as np
import pymc3 as pm
import theano as th
import theano.tensor as tt
import tqdm

from .. import types
from ..inference import fancy_optimizers
from ..inference.convergence_tracker import NoisyELBOConvergenceTracker
from ..inference.deterministic_annealing import ADVIDeterministicAnnealing
from ..inference.param_tracker import VariationalParameterTrackerConfig, VariationalParameterTracker
from ..io import io_commons
from ..models.fancy_model import GeneralizedContinuousModel

_logger = logging.getLogger(__name__)


class ConvergenceError(Exception):
    """Exception raised in case inference optimizer produces a NaN. """

    def __init__(self):
        self.message = "The optimization step for ELBO update returned a NaN"


class Sampler:
    """Base class for log emission posterior probability samplers to be used in the hybrid ADVI scheme."""
    def __init__(self, hybrid_inference_params: 'HybridInferenceParameters'):
        self.hybrid_inference_params = hybrid_inference_params

    @abstractmethod
    def update_approximation(self, approx: pm.approximations.MeanField) -> None:
        """Take a new mean-field approximation and update the sampler routine accordingly.

        Args:
            approx: an instance of PyMC3 mean-field posterior approximation

        Returns:
            None
        """
        raise NotImplementedError

    @abstractmethod
    def draw(self) -> np.ndarray:
        """Draw one sample (or average of several samples) from log emission posterior probability."""
        raise NotImplementedError

    @abstractmethod
    def reset(self):
        """Reset the internal state of the sampler (e.g. previously accumulated samples)."""
        raise NotImplementedError

    @abstractmethod
    def increment(self, update: np.ndarray):
        """Add an incremental update to the current estimate of the log emission posterior mean."""
        raise NotImplementedError

    @abstractmethod
    def get_latest_log_emission_posterior_mean_estimator(self) -> np.ndarray:
        """Returns the latest estimate of the log emission posterior mean."""
        raise NotImplementedError


class Caller:
    """Base class for callers, i.e. routines that update the posterior of discrete RVs, to be used in the
    hybrid ADVI scheme."""

    @abstractmethod
    def snapshot(self) -> None:
        """Takes a snapshot of the variables that change by `call` method. Taking a snapshot is useful if
        several calls are necessary to achieve convergence among discrete variables themselves.
        `finalize` may then admix the snapshot with the converged result."""
        raise NotImplementedError

    @abstractmethod
    def call(self) -> 'CallerUpdateSummary':
        """Update the posterior of discrete RVs and return a summary"""
        raise NotImplementedError

    @abstractmethod
    def finalize(self) -> None:
        """This method is called after exiting the call internal loop and before `update_auxiliary_vars` for
        finalizing the posteriors (e.g. admixing with the snapshot)."""

    @abstractmethod
    def update_auxiliary_vars(self) -> None:
        """Update auxiliary variables in workspaces, if any. This routine is called after one
        (or more) round(s) of invoking `Caller.call`."""
        raise NotImplemented


class CallerUpdateSummary:
    """Represents a summary of updates made to discrete RV posteriors by a `Caller`."""
    @abstractmethod
    def reduce_to_scalar(self) -> float:
        """A function that reduces arrays to scalars. It is used to summarize tensor-valued updates."""
        raise NotImplementedError

    @abstractmethod
    def __repr__(self):
        """Represents the caller update summary in a human readable format (used in logging)."""
        raise NotImplementedError


class LoggerTQDMAdapter(io.StringIO):
    """Output stream for `tqdm` which will output to logger module instead of `stderr`."""
    logger = None
    level = None
    buf = ''

    def __init__(self, logger, level=logging.INFO):
        super().__init__()
        self.logger = logger
        self.level = level

    def write(self, buf):
        self.buf = buf.strip('\r\n\t ')

    def flush(self):
        self.logger.log(self.level, self.buf)


class InferenceTask:
    """Base class of all inference tasks."""

    # Lay in a course for starbase one two warp nine point five...
    @abstractmethod
    def engage(self):
        """Initiate inference algorithm."""
        raise NotImplementedError("Core breach imminent!")

    @abstractmethod
    def disengage(self):
        """Wrap up the inference algorithm (clean up workspaces, etc.)"""
        raise NotImplementedError


class HybridInferenceTask(InferenceTask):
    """The hybrid inference framework is applicable to PGMs with the following general structure:

        +--------------+           +----------------+
        | discrete RVs + --------> + continuous RVs |
        +--------------+           +----------------+

    Note that discrete RVs do not have any continuous parents. The inference is approximately
    performed by factorizing the true posterior into an uncorrelated product of discrete RVs (DRVs)
    and continuous RVs (CRVs):

        p(CRVs, DRVs | observed) ~ q(CRVs) q(DRVs)

    q(CRVs) is updated via deterministic annealing mean-field ADVI.
    q(DRV) is updated by the provided "caller" and is out the scope of this class.

    Usage:
    ------

    Preliminaries. Let us decompose the log joint as follows:

        -log_P(CRVs, DRVs, observed) = F_c(CRVs, observed)
                                        + F_d(DRVs, observed)
                                            + F_cd(CRVs, DRVs, observed)

    The last term in the free energy (negative log joint) is the only term with cross terms between
    the discrete and continuous sectors.

    The user must supply the following components:

        (1) a pm.Model that yields the DRV-posterior-expectation of the free energy,

            F_c^{eff}(CRVs, observed) = E_{DRVs ~ q(DRVs)} [-log_P(CRVs, DRVs, observed)]
                                        = F_c(CRVs, observed)
                                            + E_{DRVs ~ q(DRVs)} [F_cd(CRVs, DRVs, observed)]
                                                + E_{DRVs ~ q(DRVs)} [F_d(DRVs, observed)]

            Note: the last term is fully determined by q(DRVs) and can be dropped while performing
            ADVI updates in the continuous sector.

        (2) a "sampler" that provides samples from the cross term, which we call "log emission",
            defined as:

            log_emission(DRVs) = E_{CRVs ~ q(CRVs)} [-F_cd (CRVs, DRV, observed)]

        (3) a "caller" that updates q(DRVs) given log_emission(DRV), i.e.:

            q(DRVs) \propto \exp[log_emission(DRVs) - F_d(DRVs, observed)]

            In practice, one does not need the complete joint posterior of DRVs: only sufficient
            statistics, to the extent required for calculating F_c^{eff} is needed. Calculating such
            sufficient statistics could be as simple as using the Bayes rule, or more complicated if
            the DRVs are strongly correlated.

    The general implementation motif is:

        (a) to store sufficient statistics from q(DRVs) as a shared theano tensor such that the the
            model can access it,
        (b) to store log_emission(DRVs) as a shared theano tensor (or ndarray) such that the caller
            can access it, and:
        (c) let the caller directly update the shared sufficient statistics.
    """

    # if the temperature is within the following tolerance of 1.0, it is assumed that annealing
    # has finished
    temperature_tolerance = 1e-6

    def __init__(self,
                 hybrid_inference_params: 'HybridInferenceParameters',
                 continuous_model: GeneralizedContinuousModel,
                 sampler: Optional[Sampler],
                 caller: Optional[Caller],
                 **kwargs):
        """Initializer.

        Args:
            hybrid_inference_params: inference configuration
            continuous_model: a PyMC3 model representing the continuous sector of the PGM
            sampler: log emission probability sampler
            caller: discrete RV posterior updater
            **kwargs: extra keywords

        Keyword Args:
            custom_optimizer: a custom stochastic optimizer to be used in place of the default optimizer (adamax
            elbo_normalization_factor: normalization factor of the full model ELBO (for logging)
            advi_task_name: name of the ADVI step (for logging)
            sampling_task_name: name of the sampling step (for logging)
            calling_task_name: name of the calling step (for logging)
        """
        assert hybrid_inference_params is not None
        self.hybrid_inference_params = hybrid_inference_params

        assert continuous_model is not None
        self.continuous_model = continuous_model
        self.continuous_model.verify_var_registry()

        if sampler is None:
            _logger.warning("No log emission sampler given; skipping the sampling step")
        self.sampler = sampler

        if caller is None:
            _logger.warning("No caller given; skipping the calling step")
        self.caller = caller

        if self.hybrid_inference_params.track_model_params:
            _logger.info("Instantiating the parameter tracker...")
            self.param_tracker = self._create_param_tracker()
        else:
            self.param_tracker = None

        _logger.info("Instantiating the convergence tracker...")
        self.advi_convergence_tracker = NoisyELBOConvergenceTracker(
            self.hybrid_inference_params.convergence_snr_averaging_window,
            self.hybrid_inference_params.convergence_snr_trigger_threshold,
            self.hybrid_inference_params.convergence_snr_countdown_window)

        _logger.info("Setting up DA-ADVI...")
        with self.continuous_model:
            if not hasattr(self, 'temperature'):
                initial_temperature = self.hybrid_inference_params.initial_temperature
                self.temperature: types.TensorSharedVariable = th.shared(
                    np.asarray([initial_temperature], dtype=types.floatX))
            initial_temperature = self.temperature.get_value()[0]
            if (np.abs(initial_temperature - 1.0) < self.temperature_tolerance or
                    hybrid_inference_params.disable_annealing):
                # no annealing
                temperature_update = None
            else:
                # linear annealing
                temperature_drop_per_iter = ((initial_temperature - 1.0) /
                                             self.hybrid_inference_params.num_thermal_advi_iters)
                temperature_update = [(self.temperature,
                                       tt.maximum(1.0, self.temperature - temperature_drop_per_iter))]

            self.continuous_model_advi = ADVIDeterministicAnnealing(
                random_seed=self.hybrid_inference_params.random_seed,
                temperature=self.temperature)
            self.continuous_model_approx: pm.MeanField = self.continuous_model_advi.approx

            if 'custom_optimizer' in kwargs.keys():
                opt = kwargs['custom_optimizer']
                assert issubclass(type(opt), fancy_optimizers.FancyStochasticOptimizer)
                self.fancy_opt = opt
            else:
                self.fancy_opt = fancy_optimizers.FancyAdamax(
                    learning_rate=hybrid_inference_params.learning_rate,
                    beta1=hybrid_inference_params.adamax_beta1,
                    beta2=hybrid_inference_params.adamax_beta2,
                    sample_specific_only=False)

            self.continuous_model_step_func = self.continuous_model_advi.objective.step_function(
                score=True,
                obj_optimizer=self.fancy_opt.get_optimizer(self.continuous_model, self.continuous_model_approx),
                total_grad_norm_constraint=self.hybrid_inference_params.total_grad_norm_constraint,
                obj_n_mc=self.hybrid_inference_params.obj_n_mc,
                more_updates=temperature_update)

        if self.sampler is not None:
            self.sampler.update_approximation(self.continuous_model_approx)

        if 'elbo_normalization_factor' in kwargs.keys():
            self.elbo_normalization_factor = kwargs['elbo_normalization_factor']
        else:
            self.elbo_normalization_factor = 1.0

        if 'advi_task_name' in kwargs.keys():
            self.advi_task_name = kwargs['advi_task_name']
        else:
            self.advi_task_name = "ADVI"

        if 'sampling_task_name' in kwargs.keys():
            self.sampling_task_name = kwargs['sampling_task_name']
        else:
            self.sampling_task_name = "sampling"

        if 'calling_task_name' in kwargs.keys():
            self.calling_task_name = kwargs['calling_task_name']
        else:
            self.calling_task_name = "calling_task_name"

        self._t0 = None
        self._t1 = None
        self.elbo_hist: List[float] = []
        self.rls_elbo_hist: List[float] = []
        self.snr_hist: List[float] = []
        self.i_epoch = 1
        self.i_advi = 1
        self.calling_hist: List[Tuple[int, bool, bool]] = []
        self.previous_sampling_rounds = 0
        self.latest_caller_update_summary: Optional[CallerUpdateSummary] = None
        self.tqdm_out = LoggerTQDMAdapter(_logger)

    @abstractmethod
    def disengage(self):
        raise NotImplementedError

    def engage(self):
        try:
            all_converged = False
            while self.i_epoch <= self.hybrid_inference_params.max_training_epochs:
                _logger.debug("Starting epoch {0}...".format(self.i_epoch))
                converged_continuous = self._update_continuous_posteriors()
                all_converged = converged_continuous
                if self.sampler is not None:
                    converged_sampling = self._update_log_emission_posterior_expectation()
                    all_converged = all_converged and converged_sampling
                if self.caller is not None:
                    converged_discrete = self._update_discrete_posteriors()
                    all_converged = all_converged and converged_discrete
                self.i_epoch += 1
                if all_converged and not self._premature_convergence():
                    break
                else:  # reset ADVI convergence tracker so that ADVI is run again
                    self.advi_convergence_tracker.reset_convergence_counter()
            if all_converged:
                _logger.info("Inference task completed successfully and convergence achieved.")
            else:
                _logger.warning("Inference task completed successfully but convergence not achieved.")
        except KeyboardInterrupt:
            pass

    def _log_start(self, task_name: str, i_epoch: int):
        self._t0 = time.time()
        _logger.debug("Starting {0} for epoch {1}...".format(task_name, i_epoch))

    def _log_stop(self, task_name: str, i_epoch: int):
        self._t1 = time.time()
        _logger.debug('The {0} for epoch {1} successfully finished in {2:.2f}s'.format(
            task_name, i_epoch, self._t1 - self._t0))

    def _log_interrupt(self, task_name: str, i_epoch: int):
        _logger.warning('The {0} for epoch {1} was interrupted'.format(task_name, i_epoch))

    def _premature_convergence(self):
        too_few_epochs = self.i_epoch < self.hybrid_inference_params.min_training_epochs
        temperature_still_high = np.abs(self.temperature.get_value()[0] - 1) > self.temperature_tolerance
        still_in_annealing = temperature_still_high and not self.hybrid_inference_params.disable_annealing
        return too_few_epochs or still_in_annealing

    def _create_param_tracker(self):
        assert all([param_name in self.continuous_model.vars or
                    param_name in self.continuous_model.deterministics
                    for param_name in self.hybrid_inference_params.param_tracker_config.var_names]),\
            "Some of the parameters chosen to be tracker are not present in the model"
        return VariationalParameterTracker(self.hybrid_inference_params.param_tracker_config)

    def _update_continuous_posteriors(self) -> bool:
        self._log_start(self.advi_task_name, self.i_epoch)
        max_advi_iters = self.hybrid_inference_params.max_advi_iter_subsequent_epochs if self.i_epoch > 1 \
            else self.hybrid_inference_params.max_advi_iter_first_epoch
        converged = False
        with tqdm.trange(max_advi_iters,
                         desc="({0}) starting...".format(self.advi_task_name),
                         file=self.tqdm_out) as progress_bar:
            try:
                for _ in progress_bar:
                    loss = self.continuous_model_step_func() / self.elbo_normalization_factor
                    if np.isnan(loss):
                        raise ConvergenceError

                    self.i_advi += 1

                    try:
                        self.advi_convergence_tracker(self.continuous_model_advi.approx, loss, self.i_advi)
                    except StopIteration:
                        if not self._premature_convergence():  # suppress signal if deemed premature
                            raise StopIteration
                        else:
                            self.advi_convergence_tracker.reset_convergence_counter()

                    snr = self.advi_convergence_tracker.snr
                    elbo_mean = self.advi_convergence_tracker.mean
                    elbo_variance = self.advi_convergence_tracker.variance
                    if snr is not None:
                        self.snr_hist.append(snr)
                    self.elbo_hist.append(-loss)
                    self.rls_elbo_hist.append(elbo_mean)
                    progress_bar.set_description("({0} epoch {1}) ELBO: {2}, SNR: {3}, T: {4:.2f}".format(
                        self.advi_task_name,
                        self.i_epoch,
                        "{0:.3f} +/- {1:.3f}".format(-elbo_mean, np.sqrt(elbo_variance))
                        if elbo_mean is not None and elbo_variance is not None else "N/A",
                        "{0:.1f}".format(snr) if snr is not None else "N/A",
                        self.temperature.get_value()[0]),
                        refresh=False)
                    if self.param_tracker is not None \
                            and self.i_advi % self.hybrid_inference_params.track_model_params_every == 0:
                        self.param_tracker(self.continuous_model_advi.approx, loss, self.i_advi)

            except StopIteration:
                converged = True
                progress_bar.refresh()
                progress_bar.close()
                self._log_stop(self.advi_task_name, self.i_epoch)

            except KeyboardInterrupt:
                progress_bar.refresh()
                progress_bar.close()
                self._log_interrupt(self.advi_task_name, self.i_epoch)
                raise KeyboardInterrupt

        return converged

    def _update_log_emission_posterior_expectation(self):
        self._log_start(self.sampling_task_name, self.i_epoch)
        if self.i_epoch == 1:
            self.sampler.reset()  # clear out log emission
        lag = min(self.previous_sampling_rounds, self.hybrid_inference_params.sampler_smoothing_window)
        converged = False
        median_rel_err = np.nan
        with tqdm.trange(self.hybrid_inference_params.log_emission_sampling_rounds,
                         desc="({0} epoch {1})".format(self.sampling_task_name, self.i_epoch),
                         file=self.tqdm_out) as progress_bar:
            try:
                for i_round in progress_bar:

                    # draw new log emission posterior samples
                    update_to_estimator = self.sampler.draw()

                    # update the estimator
                    latest_estimator = self.sampler.get_latest_log_emission_posterior_mean_estimator()
                    update_to_estimator = (update_to_estimator - latest_estimator) / (i_round + 1 + lag)
                    self.sampler.increment(update_to_estimator)

                    # relative update (and ensuring no NaNs are present, which can occur if latest_estimator
                    # has zero entries)
                    latest_estimator = self.sampler.get_latest_log_emission_posterior_mean_estimator()
                    rel_update = np.nan_to_num(np.abs(update_to_estimator / latest_estimator).flatten())
                    median_rel_err = np.median(rel_update)
                    std_rel_err = np.std(rel_update)

                    progress_bar.set_description("({0} epoch {1}) relative error: {2:2.4f} +/- {3:2.4f}".format(
                        self.sampling_task_name, self.i_epoch, median_rel_err, std_rel_err),
                        refresh=False)
                    if median_rel_err < self.hybrid_inference_params.log_emission_sampling_median_rel_error:
                        _logger.debug('{0} converged after {1} rounds with final '
                                      'median relative error {2:.3}.'.format(self.sampling_task_name, i_round + 1,
                                                                             median_rel_err))
                        raise StopIteration

            except StopIteration:
                converged = True
                progress_bar.refresh()
                progress_bar.close()
                self._log_stop(self.sampling_task_name, self.i_epoch)

            except KeyboardInterrupt:
                progress_bar.refresh()
                progress_bar.close()
                raise KeyboardInterrupt

            finally:
                self.previous_sampling_rounds = i_round + 1
                if not converged:
                    _logger.warning('{0} did not converge (median relative error '
                                    '= {1:.3}). Increase sampling rounds (current: {2}) if this behavior persists.'
                                    .format(self.sampling_task_name, median_rel_err,
                                            self.hybrid_inference_params.log_emission_sampling_rounds))

        return converged

    def _update_discrete_posteriors(self):
        self._log_start(self.calling_task_name, self.i_epoch)
        first_call_converged = False  # if convergence is achieved on the first call (stronger)
        iters_converged = False  # if internal loop is converged (weaker, does not imply global convergence)
        with tqdm.trange(self.hybrid_inference_params.max_calling_iters,
                         desc="({0} epoch {1})".format(self.calling_task_name, self.i_epoch),
                         file=self.tqdm_out) as progress_bar:
            try:
                # take a snapshot of the posteriors before running the internal convergence loop
                self.caller.snapshot()
                for i_calling_iter in progress_bar:
                    caller_summary = self.caller.call()
                    self.latest_caller_update_summary = caller_summary
                    progress_bar.set_description("({0} epoch {1}) {2}".format(
                        self.calling_task_name, self.i_epoch, repr(caller_summary)), refresh=False)
                    caller_update_size_scalar = caller_summary.reduce_to_scalar()
                    if caller_update_size_scalar < self.hybrid_inference_params.caller_update_convergence_threshold:
                        iters_converged = True
                        if i_calling_iter == 0:
                            first_call_converged = True
                        raise StopIteration

            except StopIteration:
                progress_bar.refresh()
                progress_bar.close()
                self._log_stop(self.calling_task_name, self.i_epoch)

            except KeyboardInterrupt:
                progress_bar.refresh()
                progress_bar.close()
                self._log_interrupt(self.calling_task_name, self.i_epoch)
                raise KeyboardInterrupt

            finally:
                self.caller.finalize()
                self.caller.update_auxiliary_vars()
                self.calling_hist.append((self.i_advi, iters_converged, first_call_converged))
                # if there is a self-consistency loop and not converged ...
                if not iters_converged and self.hybrid_inference_params.max_calling_iters > 1:
                    _logger.warning('{0} did not converge. Increase maximum calling rounds (current: {1}) '
                                    'if this behavior persists.'.format(self.calling_task_name,
                                                                        self.hybrid_inference_params.max_calling_iters))

        return iters_converged

    def save_elbo_history(self, output_file):
        io_commons.write_ndarray_to_tsv(output_file, np.asarray(self.elbo_hist), write_shape_info=False)


class HybridInferenceParameters:
    """Hybrid ADVI inference parameters."""
    def __init__(self,
                 learning_rate: float = 0.05,
                 adamax_beta1: float = 0.9,
                 adamax_beta2: float = 0.99,
                 obj_n_mc: int = 1,
                 random_seed: int = 1984,
                 total_grad_norm_constraint: Optional[float] = None,
                 log_emission_samples_per_round: int = 50,
                 log_emission_sampling_median_rel_error: float = 5e-3,
                 log_emission_sampling_rounds: int = 10,
                 max_advi_iter_first_epoch: int = 1000,
                 max_advi_iter_subsequent_epochs: int = 100,
                 min_training_epochs: int = 5,
                 max_training_epochs: int = 50,
                 initial_temperature: float = 2.0,
                 num_thermal_advi_iters: int = 500,
                 track_model_params: bool = False,
                 track_model_params_every: int = 10,
                 param_tracker_config: Optional['VariationalParameterTrackerConfig'] = None,
                 convergence_snr_averaging_window: int = 500,
                 convergence_snr_trigger_threshold: float = 0.1,
                 convergence_snr_countdown_window: int = 10,
                 max_calling_iters: int = 10,
                 caller_update_convergence_threshold: float = 1e-3,
                 caller_internal_admixing_rate: float = 0.75,
                 caller_external_admixing_rate: float = 0.75,
                 sampler_smoothing_window: int = 0,
                 caller_summary_statistics_reducer: Callable[[np.ndarray], float] = np.mean,
                 disable_sampler: bool = False,
                 disable_caller: bool = False,
                 disable_annealing: bool = False):
        """See `expose_args` for the description of arguments."""
        self.learning_rate = learning_rate
        self.adamax_beta1 = adamax_beta1
        self.adamax_beta2 = adamax_beta2
        self.obj_n_mc = obj_n_mc
        self.random_seed = random_seed
        self.total_grad_norm_constraint = total_grad_norm_constraint
        self.log_emission_samples_per_round = log_emission_samples_per_round
        self.log_emission_sampling_median_rel_error = log_emission_sampling_median_rel_error
        self.log_emission_sampling_rounds = log_emission_sampling_rounds
        self.max_advi_iter_first_epoch = max_advi_iter_first_epoch
        self.max_advi_iter_subsequent_epochs = max_advi_iter_subsequent_epochs
        self.min_training_epochs = min_training_epochs
        self.max_training_epochs = max_training_epochs
        self.initial_temperature = initial_temperature
        self.num_thermal_advi_iters = num_thermal_advi_iters
        self.track_model_params = track_model_params
        self.track_model_params_every = track_model_params_every
        self.param_tracker_config = param_tracker_config
        self.convergence_snr_averaging_window = convergence_snr_averaging_window
        self.convergence_snr_trigger_threshold = convergence_snr_trigger_threshold
        self.convergence_snr_countdown_window = convergence_snr_countdown_window
        self.max_calling_iters = max_calling_iters
        self.caller_update_convergence_threshold = caller_update_convergence_threshold
        self.caller_internal_admixing_rate = caller_internal_admixing_rate
        self.caller_external_admixing_rate = caller_external_admixing_rate
        self.sampler_smoothing_window = sampler_smoothing_window
        self.caller_summary_statistics_reducer = caller_summary_statistics_reducer
        self.disable_sampler = disable_sampler
        self.disable_caller = disable_caller
        self.disable_annealing = disable_annealing

        self._assert_params()

    def _assert_params(self):
        assert self.learning_rate > 0
        assert self.adamax_beta1 > 0
        assert self.adamax_beta2 > 0
        assert self.obj_n_mc > 0
        assert self.log_emission_samples_per_round > 0
        assert 0.0 < self.log_emission_sampling_median_rel_error < 1.0
        assert self.log_emission_sampling_rounds > 0
        assert self.max_advi_iter_first_epoch >= 0
        assert self.max_advi_iter_subsequent_epochs >= 0
        assert self.min_training_epochs > 0
        assert self.max_training_epochs > 0
        assert self.max_training_epochs >= self.min_training_epochs
        total_max_advi_iters = (self.max_advi_iter_first_epoch +
                                (self.max_training_epochs - 1) * self.max_advi_iter_subsequent_epochs)
        assert self.num_thermal_advi_iters <= total_max_advi_iters
        assert self.initial_temperature >= 1.0
        assert self.disable_annealing or self.num_thermal_advi_iters > 0
        assert self.track_model_params_every > 0
        assert self.convergence_snr_averaging_window > 0
        assert self.convergence_snr_trigger_threshold > 0
        assert self.max_calling_iters > 0
        assert self.caller_update_convergence_threshold > 0
        assert self.caller_internal_admixing_rate > 0
        assert self.caller_external_admixing_rate > 0
        assert self.sampler_smoothing_window >= 0

        if self.track_model_params:
            assert self.param_tracker_config is not None

        if self.disable_annealing and self.initial_temperature > 1.0:
            _logger.warning("The initial temperature ({0}) is above 1.0 and annealing is disabled. This run "
                            "makes inferences in a thermal state. Ignore this warning if "
                            "this setting is intended.".format(self.initial_temperature))

    @staticmethod
    def expose_args(args: argparse.ArgumentParser, override_default: Dict[str, Any] = None, hide: Set[str] = None):
        """Exposes arguments of `__init__` to a given instance of `ArgumentParser`.

        Args:
            args: an instance of `ArgumentParser`
            override_default: a dictionary containing arguments the default values of which
                are to be overridden before passing to the argument parser
            hide: a set of arguments not to expose

        Returns:
            None
        """
        group = args.add_argument_group(title="Inference parameters")
        if override_default is None:
            override_default = dict()
        if hide is None:
            hide = set()

        initializer_params = inspect.signature(HybridInferenceParameters.__init__).parameters
        valid_args = {"--" + arg for arg in initializer_params.keys()}
        for hidden_arg in hide:
            assert hidden_arg in valid_args, \
                "Initializer argument to be hidden {0} is not a valid initializer arguments; possible " \
                "choices are: {1}".format(hidden_arg, valid_args)
        for override_default_arg in override_default.keys():
            assert override_default_arg in valid_args, \
                "Initializer argument of which the default is to be overridden {0} is not a valid; possible " \
                "choices are: {1}".format(override_default_arg, valid_args)

        def process_and_maybe_add(arg, **kwargs):
            full_arg = "--" + arg
            if full_arg in hide:
                return
            if full_arg in override_default:
                kwargs['default'] = override_default[full_arg]
            else:
                kwargs['default'] = initializer_params[arg].default
            group.add_argument(full_arg, **kwargs)

        def str_to_bool(value: str):
            if value.lower() in ('yes', 'true', 't', 'y', '1'):
                return True
            elif value.lower() in ('no', 'false', 'f', 'n', '0'):
                return False
            else:
                raise argparse.ArgumentTypeError('Boolean value expected.')

        process_and_maybe_add("learning_rate",
                              type=float,
                              help="Adamax optimizer learning rate")

        process_and_maybe_add("adamax_beta1",
                              type=float,
                              help="Adamax first moment estimator forgetting factor")

        process_and_maybe_add("adamax_beta2",
                              type=float,
                              help="Adamax second moment estimator forgetting factor")

        process_and_maybe_add("random_seed",
                              type=int,
                              help="Random seed for the inference task")

        process_and_maybe_add("log_emission_samples_per_round",
                              type=int,
                              help="Number of log emission posterior samples per sampling round")

        process_and_maybe_add("log_emission_sampling_median_rel_error",
                              type=float,
                              help="Maximum tolerated median relative error in log emission posterior sampling")

        process_and_maybe_add("log_emission_sampling_rounds",
                              type=int,
                              help="Maximum log emission posterior sampling rounds")

        process_and_maybe_add("max_advi_iter_first_epoch",
                              type=int,
                              help="Maximum ADVI iterations in the first epoch (before sampling and calling)")

        process_and_maybe_add("max_advi_iter_subsequent_epochs",
                              type=int,
                              help="Maximum ADVI iterations after the first epoch")

        process_and_maybe_add("min_training_epochs",
                              type=int,
                              help="Minimum number of training epochs before evaluating convergence")

        process_and_maybe_add("max_training_epochs",
                              type=int,
                              help="Maximum number of training epochs before stopping")

        process_and_maybe_add("initial_temperature",
                              type=float,
                              help="Initial temperature for deterministic annealing (must be >= 1.0)")

        process_and_maybe_add("num_thermal_advi_iters",
                              type=int,
                              help="Annealing duration (in the units of ADVI iterations)")

        process_and_maybe_add("convergence_snr_averaging_window",
                              type=int,
                              help="Averaging window for calculating training SNR for evaluating convergence")

        process_and_maybe_add("convergence_snr_trigger_threshold",
                              type=float,
                              help="The SNR threshold to be reached for triggering convergence")

        process_and_maybe_add("convergence_snr_countdown_window",
                              type=int,
                              help="The number of ADVI iterations during which the SNR is required to stay below the "
                                   "set threshold for convergence")

        process_and_maybe_add("max_calling_iters",
                              type=int,
                              help="Maximum number of calling internal self-consistency iterations")

        process_and_maybe_add("caller_update_convergence_threshold",
                              type=float,
                              help="Maximum tolerated calling update size for convergence")

        process_and_maybe_add("caller_internal_admixing_rate",
                              type=float,
                              help="Admixing ratio of new and old caller posteriors (between 0 and 1; higher means "
                                   "using more of the new posterior) in internal calling loops")

        process_and_maybe_add("caller_external_admixing_rate",
                              type=float,
                              help="Admixing ratio of new and old caller posteriors (between 0 and 1; higher means "
                                   "using more of the new posterior) after internal convergence")

        process_and_maybe_add("disable_sampler",
                              type=str_to_bool,
                              help="(advanced) Disable sampler")

        process_and_maybe_add("disable_caller",
                              type=str_to_bool,
                              help="(advanced) Disable caller")

        process_and_maybe_add("disable_annealing",
                              type=str_to_bool,
                              help="(advanced) Keep the temperature pinned at its initial value and disable annealing")

    @staticmethod
    def from_args_dict(args_dict: Dict):
        """Initialize an instance of `HybridInferenceParameters` from a dictionary of arguments.

        Args:
            args_dict: a dictionary of arguments; the keys must match argument names in
                `HybridInferenceParameters.__init__`

        Returns:
            an instance of `HybridInferenceParameters`
        """
        relevant_keys = set(inspect.getfullargspec(HybridInferenceParameters.__init__).args)
        relevant_kwargs = {k: v for k, v in args_dict.items() if k in relevant_keys}
        return HybridInferenceParameters(**relevant_kwargs)
