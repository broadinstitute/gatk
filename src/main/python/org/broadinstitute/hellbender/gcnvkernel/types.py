import numpy as np
import pytensor
import pytensor.tensor as pt

# the following dtype will be used for all float numpy ndarrays
floatX = pytensor.config.floatX

# big uint dtype (used for aggregated counts)
big_uint = np.uint64

# medium uint dtype (used for single interval counts)
med_uint = np.uint32

# small uint dtype (used for copy number, contig ploidy, gc bin index)
small_uint = np.uint16

# all int dtypes
int_dtypes = [np.int8, np.int16, np.int32, np.int64]
uint_dtypes = [np.uint8, np.uint16, np.uint32, np.uint64]

# pytensor tensor types
PytensorVector = pt.TensorType(floatX, (False,))
PytensorMatrix = pt.TensorType(floatX, (False, False))
PytensorTensor3 = pt.TensorType(floatX, (False, False, False))
TensorSharedVariable = pytensor.tensor.sharedvar.TensorSharedVariable
