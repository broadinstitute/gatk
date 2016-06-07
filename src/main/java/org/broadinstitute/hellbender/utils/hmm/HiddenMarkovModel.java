package org.broadinstitute.hellbender.utils.hmm;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Hidden Markov Model interface with hidden-state transition probabilities and data
 * emission that might depends on the position.
 *
 * <p>
 * This class models a process that emits a fixed-length data sequence.
 * </p>
 * <p>
 * At each data point (or position) the process is in one of a discrete set
 * of hidden states which are shared by all positions.
 * </p>
 * <p>
 * The emission probability of a datum given the hidden state may depend on the current hidden state and also
 * the position in the sequence.
 * </p>
 * <p>
 * The transition probabilities between the current hidden state and the next may depend on these and also the current
 * and next position.
 * </p>
 * <p>
 * This extended finite HMM model is equivalent to the classical non-position dependent HMM where
 * the set of common hidden-state set is unrolled across positions.
 * </p>
 *
 * @param <D> is the observed data component.
 * @param <S> is the state component.
 * @param <T> represent the observation position type.
 */
public interface HiddenMarkovModel<D, T, S> {

    int RANDOM_SEED_FOR_CHAIN_GENERATION = 13;

    /**
     * Returns all hidden states.
     * <p>
     *     Contract:
     *     <ul>
     *         <li>This method must always return the same state list every time it is called; this cannot be null.</li>
     *         <li>There must be at least one hidden state.</li>
     *     </ul>
     * </p>
     * This method must always return the same state set with at least one element.
     *
     * @return never {@code null} nor an empty list.
     */
    List<S> hiddenStates();

    /**
     * Returns the prior probability of an state.
     * <p>
     *      Contract:
     *     <ul>
     *         <li>The return value cannot be greater than 0 (1 in linear scale).</li>
     *         <li>For a given position, the sum of probabilities returned
     *         for every hidden state must sum to 1.</li>
     *         <li>The return value cannot be a {@link Double#NaN}</li>
     *     </ul>
     * </p>
     *
     * @param state the query state.
     * @param position the query position.
     * @return a value between {@link Double#NEGATIVE_INFINITY} to 0 inclusive.
     * @throws IllegalArgumentException if {@code state} is not recognized by the model.
     */
    double logPriorProbability(final S state, final T position);

    /**
     * Returns the transition probability from one state to another.
     *
     * <p>
     *     Contract:
     *     <ul>
     *         <li>The returned value cannot be greater than 0 (or 1 in linear scale).</li>
     *         <li>Given a current position and hidden state and next position, the sum of probabilities
     *             returned for every next hidden state must sum to 1.</li>
     *         <li>The returned value cannot be a {@link Double#NaN}.</li>
     *     </ul>
     * </p>
     *
     * @param currentState the destination state.
     * @param currentPosition the source time before the transition.
     * @param nextState the source state.
     * @param nextPosition the destination time.
     * @return a value between {@link Double#NEGATIVE_INFINITY} to 0 inclusive.
     * @throws IllegalStateException if any of the input positions or states is not recognized by the model.
     */
    double logTransitionProbability(final S currentState, final T currentPosition, final S nextState, final T nextPosition);

    /**
     * Returns the emission probability of the data given hidden state and position in the sequence.
     *
     * <p>
     *     Contract:
     *     <ul>
     *         <li>The return value cannot be greater than 0 (or 1 in linear scale).</li>
     *         <li>The return value cannot be a {@link Double#NaN}.</li>
     *     </ul>
     * </p>
     *
     * @param data the observed data value.
     * @param state the hidden state at observation.
     * @param position the observation position.
     * @return a value between {@link Double#NEGATIVE_INFINITY} to 0 inclusive.
     * @throws IllegalArgumentException if either {@code state} or {@code position} is not recognized by the model.
     */
    double logEmissionProbability(final D data, final S state, final T position);

    default List<S> generateHiddenStateChain(final List<T> positions) {
        final RandomGenerator rg = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED_FOR_CHAIN_GENERATION));
        final List<S> hiddenStates = hiddenStates();
        final List<S> result = new ArrayList<>(positions.size());

        final S initialState = GATKProtectedMathUtils.randomSelect(hiddenStates, s -> Math.exp(logPriorProbability(s, positions.get(0))), rg);
        result.add(initialState);

        IntStream.range(1, positions.size()).forEach(n ->
            result.add(GATKProtectedMathUtils.randomSelect(hiddenStates,
                    s -> Math.exp(logTransitionProbability(result.get(n-1), positions.get(n - 1), s, positions.get(n))), rg)));
        return result;
    }
}
