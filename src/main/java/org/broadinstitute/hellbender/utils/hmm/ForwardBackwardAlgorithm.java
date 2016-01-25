package org.broadinstitute.hellbender.utils.hmm;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Performs the Forward-backward algorithm for
 * {@link HiddenMarkovModel}
 * on a sequence of data points and positions.
 *
 * <p>This done by calling {@link #apply} that returns an
 * {@link ForwardBackwardAlgorithm.Result} instance which
 * allow for subsequent queries on the posterior probability of the
 * hidden states at arbitrary points of the input sequence</p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ForwardBackwardAlgorithm {

    /**
     * The forward-backward algorithm result query interface.
     *
     * @param <D> the observed data type.
     * @param <T> the observation time/position type.
     * @param <S> the hidden state type.
     */
    public interface Result<D, T, S> {

        /**
         * Returns the data passed to the forward-backward algorithm.
         * @return never {@code null}.
         */
        List<D> data();

        /**
         * Returns the list of positions passed to the forward-backward algorithm
         * @return never {@code null}.
         */
        List<T> positions();

        /**
         * Returns the model passed to the foward-backward algorithm.
         * @return never {@code null}.
         */
        HiddenMarkovModel<D, T, S> model();

        /**
         * Returns the forward probability of the hidden state at certain position.
         * @param state the query hidden state.
         * @param position the query position.
         * @return a valid probability in log scale (from -Inf to 0).
         *
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logForwardProbability(final S state, final T position);

        /**
         * Returns the backward probability of the hidden state at certain position.
         * @param state the query hidden state.
         * @param position the query position.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logBackwardProbability(final S state, final T position);

        /**
         * Returns the posterior probability of the hidden state at certain position.
         * @param state the query hidden state.
         * @param position the query position.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logProbability(final S state, final T position);

        /**
         * Returns the posterior probability of a full sequence of hidden states given the original data
         * and position sequences.
         * @param states the list of states.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if {@code states}, is {@code null}, contains some value not recognized
         * as a hidden state by the model, or is not of the same length as the original data.
         */
        double logProbability(final List<S> states);

        /**
         * Returns the likelihood of the original data given the model.
         *
         * @return a valid probability in log scale (from -Inf to 0)
         */
        double logDataLikelihood();
    }

    /**
     * Runs the forward-backward algorithm on a data-positions list pair given a type-compatible model.
     * <p>
     *     Calling code can then use the returned {@link Result} instance to perform various queries based on
     *     the calculations performed by the algorithm. Refer to {@link Result} documentation for details.
     * </p>
     * <p>
     *     This method invocation may be computationally intensive with large data-sets. The complexity of the
     *     algorithm is {@code O(L*(N^2))} where {@code L is the length of the data/position sequence and
     *     {@code N} is the number of hidden states. The memory size of the result is of the order of {@code O(L*N)}.
     * </p>
     *
     * @param data the observed data sequence.
     * @param positions the observation time/position points.
     * @param model the HMM model.
     * @param <D> the observed data type.
     * @param <T> the observation time/position type.
     * @param <S> the hidden state type.
     * @return never {@code null}.
     * @throws IllegalArgumentException if any of the arguments is {@code null}, {@code data} and {@code positions}
     *   have different length, or if the {@code model} does not recognized any of the values in {@code data} or {@code}
     *   positions.
     */
    public static <D, T, S> Result<D, T, S> apply(final List<D> data, final List<T> positions,
                                 final HiddenMarkovModel<D, T, S> model) {
        Utils.nonNull(data, "the input data sequence cannot be null.");
        Utils.nonNull(positions, "the input position sequence cannot be null.");
        Utils.nonNull(model, "the input model cannot be null");

        final List<D> dataList = Collections.unmodifiableList(new ArrayList<>(data));
        final List<T> positionList = Collections.unmodifiableList(new ArrayList<>(positions));
        if (dataList.size() != positionList.size()) {
            throw new IllegalArgumentException("the data sequence and position sequence must have the same number of elements");
        }

        final double[][] forwardProbabilities = calculateLogForwardProbabilities(model, dataList, positionList);
        final double[][] backwardProbabilities = calculateLogBackwardProbabilities(model, dataList, positionList);

        return new ArrayResult<>(dataList, positionList, model, forwardProbabilities, backwardProbabilities);
    }

    /**
     * Calculates the forward probabilities (the forward phase of the algorithm). These are in log scale.
     * <p>
     *     The returned array is uni-dimensional. Forward probabilities are arranged by position, thus the
     *     forward probability for the ith hidden state and jth position is stored in the
     *     index: i + j * number-of-hidden-states.
     * </p>
     * @param model the HMM model.
     * @param data the observed data sequence.
     * @param positions the observation time/position points.
     * @param <D> the observed data type.
     * @param <T> the observation time/position type.
     * @param <S> the hidden state type.
     * @return never {@code null}, the calling code can modify this array at will.
     * @throws IllegalArgumentException if the {@code model} does not recognize any of the values in {@code data}
     *   or {@code positions}.
     */
    private static <D, T, S> double[][] calculateLogForwardProbabilities(
            final HiddenMarkovModel<D, T, S> model,
            final List<D> data,
            final List<T> positions) {

        final List<S> states = model.hiddenStates();
        final int numStates = states.size();
        final int length = data.size();
        final double[][] result = new double[length][numStates];
        // Empty data? then we just return the empty array:
        if (length == 0) {
            return result;
        }
        // We first initialize the probabilities at the first position.
        final T position0 = positions.get(0);
        final D datum0 = data.get(0);
        for (int stateIndex = 0; stateIndex < states.size(); stateIndex++) {
            final S state = states.get(stateIndex);
            result[0][stateIndex] =
                    model.logPriorProbability(state, position0)
                + model.logEmissionProbability(datum0, state, position0);
        }

        // Then we do the rest t_1, t_2, ... and so on.
        // Array re-used to hold the elements of log-sum-exp operations:
        final double[] logSumBuffer = new double[states.size()];
        for (int thisPositionIndex = 1; thisPositionIndex < data.size(); thisPositionIndex++) {
            final int previousPositionIndex = thisPositionIndex - 1;
            final T previousPosition = positions.get(previousPositionIndex);
            final T thisPosition = positions.get(thisPositionIndex);
            for (int thisStateIndex = 0; thisStateIndex < numStates; thisStateIndex++) {
                final S thisState = states.get(thisStateIndex);
                for (int previousStateIndex = 0; previousStateIndex < numStates; previousStateIndex++) {
                    logSumBuffer[previousStateIndex] =
                            result[previousPositionIndex][previousStateIndex]
                                    + model.logTransitionProbability(states.get(previousStateIndex),
                                        previousPosition, thisState, thisPosition);
                }
                result[thisPositionIndex][thisStateIndex] = GATKProtectedMathUtils.naturalLogSumNaturalLog(logSumBuffer)
                        + model.logEmissionProbability(data.get(thisPositionIndex), thisState, thisPosition);
            }
        }
        return result;
    }

    /**
     * Calculates the backward probabilities (the backward phase of the algorithm). These are in log scale.
     * <p>
     *     The returned array is uni-dimensional. Back probabilities are arranged by position, thus the
     *     back probability for the ith hidden state and jth position is stored in the
     *     index: i + j * number-of-hidden-states.
     * </p>
     * @param model the HMM model.
     * @param dataList the observed data sequence.
     * @param positionList the observation time/position points.
     * @param <D> the observed data type.
     * @param <T> the observation time/position type.
     * @param <S> the hidden state type.
     * @return never {@code null}, the calling code can modify this array at will.
     * @throws IllegalArgumentException if the {@code model} does not recognize any of the values in {@code dataList}
     *   or {@code positionList}.
     */
    private static <D, T, S> double[][] calculateLogBackwardProbabilities(final HiddenMarkovModel<D, T, S> model,
                                                                        final List<D> dataList,
                                                                        final List<T> positionList) {

        final List<S> states = model.hiddenStates();
        final int numStates = states.size();
        final int length = dataList.size();

        // result contains is implicitly initialized to all 0.
        final double[][] result = new double[length][numStates];

        // Empty data? then we just return the empty array:
        if (length == 0) {
            return result;
        }

        // result already contain the correct value for the
        // last position in the input data/position sequence.
        // = 0 (i.e. log(1))
        // thus we proceed directly to t_L - 1.

        // "small" buffer array reused to do the log-sum-exp-log trick:
        final double[] logSumBuffer = new double[states.size()];

        for (int thisPositionIndex = length - 2; thisPositionIndex >= 0; --thisPositionIndex) {
            final int nextPositionIndex = thisPositionIndex + 1;
            final T thisPosition = positionList.get(thisPositionIndex);
            final T nextPosition = positionList.get(nextPositionIndex);
            for (int thisStateIndex = 0; thisStateIndex < numStates; thisStateIndex++) {
                final S thisState = states.get(thisStateIndex);
                for (int nextStateIndex = 0; nextStateIndex < numStates; nextStateIndex++) {
                    logSumBuffer[nextStateIndex] = result[nextPositionIndex][nextStateIndex]
                                    + model.logTransitionProbability(thisState, thisPosition,
                                                states.get(nextStateIndex), nextPosition)
                                    + model.logEmissionProbability(dataList.get(nextPositionIndex),
                                                states.get(nextStateIndex), nextPosition);
                }
                result[thisPositionIndex][thisStateIndex] = GATKProtectedMathUtils.naturalLogSumNaturalLog(logSumBuffer);
            }
        }
        return result;
    }

    /**
     * Implementation of the interface {@link Result} returned by the {@link #apply} method.
     * @param <D> the observed data type.
     * @param <T> the observation time/position type.
     * @param <S> the hidden state type.
     */
    private static class ArrayResult<D, T, S> implements Result<D, T, S> {

        private final List<D> data;

        private final List<T> positions;

        private final HiddenMarkovModel<D, T, S> model;

        private final Object2IntMap<T> positionIndex;
        private final Object2IntMap<S> stateIndex;

        private final int numStates;

        private final double[][] logForwardProbabilities;

        private final double[][] logBackwardProbabilities;
        private final double logDataLikelihood;

        private ArrayResult(final List<D> data, final List<T> positions,
                            final HiddenMarkovModel<D, T, S> model,
                            final double[][] logForwardProbabilities,
                            final double[][] logBackwardProbabilities) {
            this.data = Collections.unmodifiableList(new ArrayList<>(data));
            this.positions = Collections.unmodifiableList(new ArrayList<>(positions));
            this.model = model;
            this.positionIndex = composeIndexMap(this.positions);
            this.stateIndex = composeIndexMap(model.hiddenStates());
            this.numStates = stateIndex.size();
            this.logBackwardProbabilities = logBackwardProbabilities;
            this.logForwardProbabilities = logForwardProbabilities;
            this.logDataLikelihood = calculateLogDataLikelihood(logForwardProbabilities, logBackwardProbabilities);
        }

        /**
         * Calculates the data likelihood in log scale.
         * <p>
         *     This value can be obtained by adding up the posterior probabilities at any position; their sum is supposed
         *     to be the same across at each position.
         * </p>
         * <p>
         *     This implementation uses always uses the first position.
         * </p>
         * @param logForwardProbabilities the log forward probabilities array.
         * @param logBackwardProbabilities the log backward probabilities array.
         * @return a valid probability in log scale (between -Inf and 0 inclusive).
         */
        private static double calculateLogDataLikelihood(final double[][] logForwardProbabilities,
                                                         final double[][] logBackwardProbabilities) {
            if (logForwardProbabilities.length == 0) {
                return 0;
            } else {
                return GATKProtectedMathUtils.naturalLogSumNaturalLog(IntStream.range(0, logForwardProbabilities[0].length)
                        .mapToDouble(i -> logBackwardProbabilities[0][i] + logForwardProbabilities[0][i]).toArray());
            }
        }

        /**
         * Composes a object ot index map given an object list.
         * @param list the list to index.
         * @param <E> the element type.
         * @return never {@code null}.
         */
        private <E> Object2IntMap<E> composeIndexMap(final List<E> list) {
            return IntStream.range(0, list.size())
                    .collect(
                            () -> new Object2IntOpenHashMap<>(list.size()),
                            (map, i) -> map.put(list.get(i), i),
                            (map1, map2) -> map2.object2IntEntrySet().forEach(
                                    e -> map1.put(e.getKey(), e.getIntValue())
                            ));
        }

        @Override
        public List<D> data() {
            return data;
        }

        @Override
        public List<T> positions() {
            return positions;
        }

        @Override
        public HiddenMarkovModel<D, T, S> model() {
            return model;
        }

        @Override
        public double logForwardProbability(final S state, final T position) {
            final int stateIndex = this.stateIndex.getOrDefault(state, -1);
            final int positionIndex = this.positionIndex.getOrDefault(position, -1);
            if (stateIndex == -1) {
                throw new IllegalArgumentException("the input state is not recognized by the model");
            } else if (positionIndex == -1) {
                throw new IllegalArgumentException("unknown input position");
            } else {
                return logForwardProbabilities[positionIndex][stateIndex];
            }
        }

        @Override
        public double logBackwardProbability(final S state, final T position) {
            final int stateIndex = this.stateIndex.getOrDefault(state, -1);
            final int positionIndex = this.positionIndex.getOrDefault(position, -1);
            if (stateIndex == -1) {
                throw new IllegalArgumentException("the input state is not recognized by the model");
            } else if (positionIndex == -1) {
                throw new IllegalArgumentException("unknown input position");
            } else {
                return logBackwardProbabilities[positionIndex][stateIndex];
            }
        }

        @Override
        public double logProbability(final S state, final T position) {
            final int stateIndex = this.stateIndex.getOrDefault(state, -1);
            final int positionIndex = this.positionIndex.getOrDefault(position, -1);
            if (stateIndex == -1) {
                throw new IllegalArgumentException("the input state is not recognized by the model");
            } else if (positionIndex == -1) {
                throw new IllegalArgumentException("unknown input position");
            } else {
                return logBackwardProbabilities[positionIndex][stateIndex]
                        + logForwardProbabilities[positionIndex][stateIndex] - logDataLikelihood;
            }
        }

        @Override
        public double logProbability(final List<S> states) {
            Utils.nonNull(states, "the input states sequence cannot be null");
            if (states.size() != positions.size()) {
                throw new IllegalArgumentException("the input states sequence does not have the same length as the data sequence");
            } else if (states.size() == 0) {
                return 0;
            } else {
                final int lastIndex = states.size() - 1;
                double result = logForwardProbability(states.get(0), positions.get(0));
                result += logBackwardProbability(states.get(lastIndex), positions.get(lastIndex));
                for (int i = 1; i <= lastIndex; i++) {
                    result += model.logTransitionProbability(states.get(i - 1), positions.get(i - 1), states.get(i), positions.get(i));
                    result += model.logEmissionProbability(data.get(i), states.get(i), positions.get(i));
                }
                result -= logDataLikelihood;
                return result;
            }
        }

        @Override
        public double logDataLikelihood() {
            return logDataLikelihood;
        }
    }

}
