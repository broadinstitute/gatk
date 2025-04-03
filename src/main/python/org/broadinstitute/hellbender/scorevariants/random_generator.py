import numpy as np

# Implementation of random number generator to get the same result as GATK CNNScoreVariants reservoir down-sampler

# Reservoir Downsampler: Selects n reads out of a stream whose size is not known in advance, with
# every read in the stream having an equal chance of being selected for inclusion.

# An implementation of "Algorithm R" from the paper "Random Sampling with a Reservoir" (Jeffrey Scott Vitter, 1985)

class Random:
    multiplier = np.int64(0)
    addend = np.int64(0)
    mask = np.int64(0)

    def __init__(self, newSeed):
        self.setSeed(newSeed)

    def setSeed(self, newSeed):
        self.seed = np.int64((newSeed ^ Random.multiplier) & Random.mask)

    def next(self, bits):
        self.seed = np.int64((self.seed * Random.multiplier + Random.addend) & Random.mask)
        return np.int32(self.seed >> (48 - bits))

    def nextInt(self):
        return np.int32(self.next(32))

    def nextInt(self, n):
        if (n & -n) == n:  # i.e., n is a power of 2
            return np.int32(n * self.next(31) >> 31)

        val = 0
        while True:
            bits = self.next(31)
            val = np.int32(bits % n)
            if np.int32(bits - val + (n-1)) >= 0:
                break
        return val
