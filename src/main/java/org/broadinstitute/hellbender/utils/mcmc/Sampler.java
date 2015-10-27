package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

/**
 * Interface for generating random samples of a quantity,
 * given an {@link AbstractParameterizedState} and a {@link DataCollection}.
 * @param <U>   type of quantity to be sampled
 * @param <S>   type of AbstractParameterizedState
 * @param <T>   type of DataCollection
 */
public interface Sampler<U, S extends AbstractParameterizedState, T extends DataCollection> {
    /**
     * Returns a random sample of a quantity that is dependent on an {@link AbstractParameterizedState} and a
     * {@link DataCollection}.
     * @param rng               RandomGenerator to use in generating random sample
     * @param state             AbstractParameterizedState to use in generating random sample
     * @param dataCollection    DataCollection to use in generating random sample
     * @return                  random sample of quantity
     */
    U sample(final RandomGenerator rng, final S state, final T dataCollection);
}