import pytensor


def change_test_value_flag(mode: str):
    """Returns a decorator that sets PyTensor's `compute_test_value` flag for the decorated callable.

    PyTensor removed test values, along with the `compute_test_value` configuration flag, in 3.0.0.
    Setting the flag there raises a KeyError, so the decorator is a no-op with those versions.

    Args:
        mode: value for the `compute_test_value` flag, e.g. "off" or "ignore"

    Returns:
        a decorator
    """
    if hasattr(pytensor.config, "compute_test_value"):
        return pytensor.config.change_flags(compute_test_value=mode)

    def identity_decorator(func):
        return func

    return identity_decorator


try:
    from pytensor.tensor.shape import unbroadcast
except ImportError:
    def unbroadcast(x, *axes):
        """Returns `x` unchanged.

        PyTensor removed broadcastable flags, along with `unbroadcast`, in 3.0.0: a dimension of
        static shape 1 is no longer implicitly broadcastable there, so nothing has to be undone.
        """
        return x


try:
    from pytensor.sparse import SparseConstant
except ImportError:
    # PyTensor renamed SparseConstant to Constant in 3.0.0.
    from pytensor.sparse import Constant as SparseConstant


try:
    from pytensor import reduce as scan_reduce
except ImportError:
    def scan_reduce(fn, sequences, outputs_info, non_sequences=None, go_backwards=False,
                    mode=None, name=None):
        """Constructs a `Scan` `Op` that functions like `reduce`, returning only the last step.

        PyTensor removed the scan views (`reduce`, `foldl`, `foldr`, `map`) in 3.0.0; this is the
        implementation they had, expressed with `pytensor.scan`.
        """
        outputs, updates = pytensor.scan(fn=fn,
                                         sequences=sequences,
                                         outputs_info=outputs_info,
                                         non_sequences=non_sequences,
                                         go_backwards=go_backwards,
                                         truncate_gradient=-1,
                                         mode=mode,
                                         name=name)
        if isinstance(outputs, (list, tuple)):
            return [output[-1] for output in outputs], updates
        return outputs[-1], updates
