package org.broadinstitute.hellbender.utils.mcmc;


import org.apache.commons.math3.random.RandomGenerator;

/**
 * Interface for generating random samples of a {@link Parameter} value,
 * given an {@link ParameterizedState} and a {@link DataCollection}.
 * @param <U>   type of parameter value to be sampled
 * @param <V>   type of enumerated parameters, see {@link ParameterEnum}
 * @param <S>   type of {@link ParameterizedState}
 * @param <T>   type of {@link DataCollection}
 */
@FunctionalInterface
public interface ParameterSampler<U, V extends Enum<V> & ParameterEnum, S extends ParameterizedState<V>, T extends DataCollection> {
    /**
     * Returns a random sample of a value that is dependent on an {@link ParameterizedState} and a
     * {@link DataCollection}.
     * @param rng               RandomGenerator to use in generating random sample
     * @param state             ParameterizedState to use in generating random sample
     * @param dataCollection    DataCollection to use in generating random sample
     * @return                  random sample of value
     */
    U sample(final RandomGenerator rng, final S state, final T dataCollection);
}