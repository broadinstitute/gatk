import numpy as np
from pymc3.variational.callbacks import Callback

from ..utils.rls import NonStationaryLinearRegression


class NoisyELBOConvergenceTracker(Callback):
    """Convergence stopping criterion based on the linear trend of the noisy ELBO observations."""

    MIN_WINDOW_SIZE = 10

    def __init__(self,
                 window: int = 100,
                 snr_stop_trigger_threshold: float = 0.5,
                 stop_countdown_window: int = 10):
        """Initializer.

        Args:
            window: window size for performing linear regression
            snr_stop_trigger_threshold: signal-to-noise (SNR) ratio threshold for triggering countdown to stop
            stop_countdown_window: once the trigger is pulled, SNR must remain below the given threshold for at
                least `stop_countdown_window` subsequent ELBO observations to raise StopIteration. the countdown
                will be reset if at any point the snr goes about `snr_stop_trigger_threshold`.
        """
        self.window = window
        self.snr_stop_trigger_threshold = snr_stop_trigger_threshold
        self.stop_countdown_window = stop_countdown_window

        assert self.window > self.MIN_WINDOW_SIZE, \
            "ELBO linear regression window size is too small (minimum is {0})".format(self.MIN_WINDOW_SIZE)
        assert self.snr_stop_trigger_threshold > 0, "Bad SNR stop trigger threshold (must be positive)"
        assert self.stop_countdown_window >= 1, "Bad SNR-under-threshold countdown window (must be >= 1)"

        self._lin_reg = NonStationaryLinearRegression(window=self.window)
        self._n_obs: int = 0
        self._n_obs_snr_under_threshold: int = 0

        self.egpi: float = None  # effective gain per iteration
        self.snr: float = None  # signal-to-noise ratio
        self.variance: float = None  # variance of elbo in the window
        self.drift: float = None  # absolute elbo change in the window
        self.mean: float = None  # average elbo in the window

    def __call__(self, approx, loss, i):
        """Add new iteration.

        Args:
            approx: approximation (ignored)
            loss: current value of the loss function
            i: current iteration number (ignored)

        Returns:
            None
        """
        self._lin_reg.add_observation(loss)
        self._n_obs += 1
        self.egpi = self._lin_reg.get_slope()
        self.mean = self._lin_reg.get_intercept()
        self.variance = self._lin_reg.get_variance()
        if self.egpi is not None and self.variance is not None:
            self.egpi *= -1
            self.variance = np.abs(self.variance)
            self.drift = np.abs(self.egpi) * self.window
            self.snr = self.drift / np.sqrt(2 * self.variance)
            if self.snr < self.snr_stop_trigger_threshold and i > self.window:
                self._n_obs_snr_under_threshold += 1
            else:  # reset countdown
                self._n_obs_snr_under_threshold = 0
            if self._n_obs_snr_under_threshold >= self.stop_countdown_window:
                raise StopIteration("Convergence criterion satisfied: SNR remained below {0} for "
                                    "{1} iterations.".format(self.snr_stop_trigger_threshold,
                                                             self.stop_countdown_window))

    def reset_convergence_counter(self):
        self._n_obs_snr_under_threshold = 0
