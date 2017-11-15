# let theano share memory workspace on large tensors with numpy
borrow_numpy = True

# used for stabilizing log calculations
log_eps = 1e-12

# if a normalized PMF violates total probability by the following threshold, it will
# be normalized and a warning will be emitted
prob_sum_tol = 1e-10
