import numpy as np
from collections import deque
from typing import Deque, Tuple, Optional


class NonStationaryLinearRegression:
    """This class performs maximum-likelihood linear regression for sequentially observed data
    on an equally spaced grid. The data can be non-stationary, and a window size needs
    to be provided that determined the forgetting factor. This is a non-recursive implementation of
    the recursive least squares (RLS) algorithm.
    """

    def __init__(self, window: int = 50):
        assert window >= 1, "Window size must be >= 1"
        self._n_obs: int = 0
        self._window: int = window
        self._obs_buffer: Deque[float] = deque([], window)
        self._slope: Optional[float] = None
        self._intercept: Optional[float] = None
        self._variance: Optional[float] = None
        self._weights: Optional[np.ndarray] = None
        self._summed_weights: Optional[float] = None
        self._x: Optional[np.ndarray] = None
        self._kern: Optional[np.ndarray] = None

    @staticmethod
    def _get_lambda(window: int) -> float:
        """Convert averaging window length to forgetting factor."""
        return np.exp(-np.log(2.0) / window)

    @staticmethod
    def _generate_auxiliary_vars(lam: float, n: int) -> Tuple[np.ndarray, float, np.ndarray, np.ndarray]:
        assert n >= 2, "At least two observations are required"
        weights = np.asarray([lam ** (n - i) for i in range(1, n + 1)])
        summed_weights: float = np.sum(weights)
        x = np.arange(1, n + 1)
        kern: np.ndarray = np.linalg.inv(np.asarray(
            [[np.sum(weights), np.sum(weights * x)],
             [np.sum(weights * x), np.sum(weights * x * x)]]))
        return weights, summed_weights, x, kern

    def add_observation(self, y: float):
        self._n_obs += 1
        self._obs_buffer.append(y)
        if self._n_obs < 2:
            return
        elif 2 <= self._n_obs <= self._window:
            self._weights, self._summed_weights, self._x, self._kern = self._generate_auxiliary_vars(
                self._get_lambda(self._n_obs), self._n_obs)
        self._update_regression()

    def _update_regression(self):
        y = np.asarray(self._obs_buffer)
        wy = np.sum(self._weights * y)
        wxy = np.sum(self._weights * self._x * y)
        sol = np.dot(self._kern, np.asarray([wy, wxy]).reshape((2, 1)))
        self._intercept = sol[0, 0]
        self._slope = sol[1, 0]
        dev = y - self._intercept - self._slope * self._x
        self._variance = np.sum(self._weights * (dev ** 2)) / self._summed_weights

    def get_slope(self) -> Optional[float]:
        """Get the latest estimate of the linear regression slope.

        Returns:
            float value of the slope if estimate is available
            None if the estimate is not available yet
        """
        return self._slope

    def get_intercept(self) -> Optional[float]:
        """Get the latest estimate of the linear regression intercept.

        Returns:
            float value of the intercept if estimate is available
            None if the estimate is not available yet
        """
        return self._intercept

    def get_variance(self) -> Optional[float]:
        """Get the latest estimate of the linear regression variance.

        Returns:
            float value of the variance if estimate is available
            None if the estimate is not available yet
        """
        return self._variance
