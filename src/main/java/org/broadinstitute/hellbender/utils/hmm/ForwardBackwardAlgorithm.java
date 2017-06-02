package org.broadinstitute.hellbender.utils.hmm;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.lang.math.IntRange;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Performs the Forward-backward algorithm for
 * {@link HMM}
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
         * Returns the model passed to the forward-backward algorithm.
         * @return never {@code null}.
         */
        HMM<D, T, S> model();


        /**
         * Returns the forward probability of the hidden state at certain position.
         * @param position the query position index.
         * @param state the query hidden state.
         * @return a valid probability in log scale (from -Inf to 0).
         *
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logForwardProbability(final int position, final S state);

        /**
         * Returns the forward probability of the hidden state at certain position.
         * @param position the query position.
         * @param state the query hidden state.
         * @return a valid probability in log scale (from -Inf to 0).
         *
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logForwardProbability(final T position, final S state);

        /**
         * Returns the backward probability of the hidden state at certain position.
         * @param position the query position index.
         * @param state the query hidden state.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logBackwardProbability(final int position, final S state);

        /**
         * Returns the backward probability of the hidden state at certain position.
         * @param position the query position.
         * @param state the query hidden state.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logBackwardProbability(final T position, final S state);

        /**
         * Returns the posterior probability of the hidden state at certain position.
         * @param position the query position index.
         * @param state the query hidden state.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logProbability(final int position, final S state);

        /**
         * Returns the posterior probability of the hidden state at certain position.
         * @param position the query position.
         * @param state the query hidden state.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if either the {@code state} or the {@code position} are not
         *  recognized by the model or were not passed to the forward-backward algorithm.
         */
        double logProbability(final T position, final S state);

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
         * Returns the posterior probability that the hidden state takes on a constant value along
         * a position interval.
         * @param from the starting position index for the query interval (inclusive).
         * @param to the stop position index for the query interval (exclusive).
         * @param state the query constant hidden state.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if {@code state}, is not a valid state in the underlying model or
         * the {@code from} and {@code to} do not represent a valid position index range.
         */
        default double logProbability(final int from, final int to, final S state) {
            ParamUtils.inRange(to, from, positions().size(), "the 'to' index must be between 'from' and the length of the data/position sequence");
            return logProbability(from, Collections.nCopies(to - from, state));
        }

        /**
         * Returns the posterior probability of a sequence of hidden states starting at a particular position.
         * @param from the starting position index, so that the first element in {@code states} corresponds to the
         *             state at this position.
         * @param states the list of states.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if {@code states}, is {@code null}, contains some value not recognized
         * as a hidden state by the model, or the {@code states} sequence is too long and would go beyond the last
         * position.
         */
        double logProbability(final int from, final List<S> states);

        /**
         * Returns the posterior probability of a sequence of hidden states starting at a particular position.
         * @param from the starting position, so that the first element in {@code states} corresponds to the
         *             state at this position.
         * @param states the list of states.
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if {@code states}, is {@code null}, contains some value not recognized
         * as a hidden state by the model, or the {@code states} sequence is too long and would go beyond the last
         * position.
         */
        double logProbability(final T from, final List<S> states);

        /**
         * Return the posterior probability of a sequence of hidden state constraints.
         *
         * @param from the first position index, that corresponds to the first constraint in {@code stateConstraints}
         * @param stateConstraints the query constrain sequence.
         * @return a log scaled probability between 0 and -Inf.
         * @throws IllegalArgumentException if {@code stateConstraints} is {@code null}, it contains a {@code null}
         *  set or any of the states in any of the sets in {@code stateConstraints} is not recognized by the model or
         *  {@code from} is not a valid position in the original data/position sequence.
         */
        double logConstrainedProbability(final int from, final List<Set<S>> stateConstraints);

        /**
         * Return the posterior probability of a sequence of hidden state constraints.
         *
         * @param from the first target, that corresponds to the first constraint in {@code stateConstraints}
         * @param stateConstraints the query constrain sequence.
         * @return a log scaled probability between 0 and -Inf.
         * @throws IllegalArgumentException if {@code stateConstraints} is {@code null}, it contains a {@code null}
         *  set or any of the states in any of the sets in {@code stateConstraints} is not recognized by the model or
         *  {@code from} is not a valid position in the original data/position sequence.
         */
        double logConstrainedProbability(final T from, final List<Set<S>> stateConstraints);

        /**
         * Return the posterior joint probability of going from one state to another in a specified interval
         * @param from the starting position (inclusive)
         * @param to the stop position (inclusive)
         * @param fromState the initial state
         * @param toState the final state
         * @return a valid probability in log scale (from -Inf to 0)
         * @throws IllegalArgumentException if {@code fromState} or {@code toState} are not valid states in the
         * underlying model or the {@code from} and {@code to} do not represent a valid position index range.
         */
        default double logJointProbability(final int from, final int to, final S fromState, final S toState) {
            ParamUtils.inRange(from, 0, positions().size() - 1, "The 'from' index must be between 0 and the length" +
                    " of the data/position sequence - 1");
            ParamUtils.inRange(to, from + 1, positions().size() - 1, "the 'to' index must be between 'from' + 1 and" +
                    " the length of the data/position sequence - 1");
            final List<Set<S>> paths = new ArrayList<>(to - from + 1);
            paths.add(Collections.singleton(fromState));
            if (to > from + 1) { /* intermediate states */
                paths.addAll(Collections.nCopies(to - from - 1, new HashSet<>(model().hiddenStates())));
            }
            paths.add(Collections.singleton(toState));
            return logConstrainedProbability(from, paths);
        }

        /**
         * Returns the likelihood of the original data given the model.
         *
         * @return a valid probability in log scale (from -Inf to 0)
         */
        double logDataLikelihood();

        /**
         * Returns the likelihood of the original data given the model as evaluated at a particular
         * position.
         *
         * <p>
         * Although in theory one does not need to indicate a position to evaluate the data-likelihood
         * in practice, lack of float point precision may make those small differences important (e.g.
         * some computations with small probabilities may result in a probability just over 1.0 which makes
         * no sense).
         * </p>
         * <p>
         * As a rule of thumb, to avoid such problems you should request the data-likelihood at the target
         * from which you evaluate the backward probability.
         * </p>
         *
         * @param position the position.
         * @return a valid log scaled probability from -Inf to 0.
         * @throws IllegalArgumentException if {@code target} is unknown to the model.
         */
        double logDataLikelihood(final int position);

        /**
         * Returns the likelihood of the original data given the model as evaluated at a particular
         * position.
         *
         * <p>
         * Although in theory one does not need to indicate a position to evaluate the data-likelihood
         * in practice, lack of float point precision may make those small differences important (e.g.
         * some computations with small probabilities may result in a probability just over 1.0 which makes
         * no sense).
         * </p>
         * <p>
         * As a rule of thumb, to avoid such problems you should request the data-likelihood at the target
         * from which you evaluate the backward probability.
         * </p>
         *
         * @param target the position.
         * @return a valid log scaled probability from -Inf to 0.
         * @throws IllegalArgumentException if {@code target} is unknown to the model.
         */
        double logDataLikelihood(final T target);

        /**
         * Returns the posterior probability of the hidden chain, defined as:
         *
         *   \log \pi_i E[z_{i,0}] + \log T_{t, t+1}^{i, j} E[z_{i,t} z_{j,t+1}]
         *
         * This quantity is the logDataLikelihood excluding the emission probability
         *
         * @return a log scaled probability from -Inf to 0.
         */
        double logChainPosteriorProbability();

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
                                 final HMM<D, T, S> model) {
        Utils.nonNull(data, "the input data sequence cannot be null.");
        Utils.nonNull(positions, "the input position sequence cannot be null.");
        Utils.nonNull(model, "the input model cannot be null");

        final List<D> dataList = Collections.unmodifiableList(new ArrayList<>(data));
        final List<T> positionList = Collections.unmodifiableList(new ArrayList<>(positions));
        Utils.validateArg(dataList.size()== positionList.size(), "the data sequence and position sequence must have the same number of elements");

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
            final HMM<D, T, S> model,
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
                result[thisPositionIndex][thisStateIndex] = GATKProtectedMathUtils.logSumExp(logSumBuffer)
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
    private static <D, T, S> double[][] calculateLogBackwardProbabilities(final HMM<D, T, S> model,
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

        // "small" buffer array reused to do the log-sum-exp trick:
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
                result[thisPositionIndex][thisStateIndex] = GATKProtectedMathUtils.logSumExp(logSumBuffer);
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
    private static class ArrayResult<D, T, S> implements Result<D, T, S>, Serializable {

        private static final long serialVersionUID = -8556604447304292642L;

        private final List<D> data;

        private final List<T> positions;
        private final IntRange positionIndexRange;

        private final HMM<D, T, S> model;

        private final Object2IntMap<T> positionIndex;
        private final Object2IntMap<S> stateIndex;

        private final double[][] logForwardProbabilities;

        private final double[][] logBackwardProbabilities;
        private final double[] logDataLikelihood;

        private ArrayResult(final List<D> data, final List<T> positions,
                            final HMM<D, T, S> model,
                            final double[][] logForwardProbabilities,
                            final double[][] logBackwardProbabilities) {
            this.data = Collections.unmodifiableList(new ArrayList<>(data));
            this.positions = Collections.unmodifiableList(new ArrayList<>(positions));
            positionIndexRange = new IntRange(0, positions.size() - 1);
            this.model = model;
            positionIndex = composeIndexMap(this.positions);
            stateIndex = composeIndexMap(model.hiddenStates());
            this.logBackwardProbabilities = logBackwardProbabilities;
            this.logForwardProbabilities = logForwardProbabilities;
            logDataLikelihood = calculateLogDataLikelihood(logForwardProbabilities, logBackwardProbabilities);
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
        private static double[] calculateLogDataLikelihood(final double[][] logForwardProbabilities,
                                                           final double[][] logBackwardProbabilities) {
            return IntStream.range(0, logForwardProbabilities.length)
                    .mapToObj(i ->
                            IntStream.range(0, logForwardProbabilities[i].length)
                                    .mapToDouble(j -> logBackwardProbabilities[i][j]
                                                    + logForwardProbabilities[i][j])
                                    .toArray())
                    .mapToDouble(GATKProtectedMathUtils::logSumExp)
                    .toArray();
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
        public HMM<D, T, S> model() {
            return model;
        }

        @Override
        public double logForwardProbability(final int positionIndex, final S state) {
            ParamUtils.inRange(positionIndexRange, positionIndex, "position index");
            final int stateIndex = validStateIndex(state);
            return logForwardProbabilities[positionIndex][stateIndex];
        }

        @Override
        public double logForwardProbability(final T position, final S state) {
            final int positionIndex = validPositionIndex(position);
            final int stateIndex = validStateIndex(state);
            return logForwardProbabilities[positionIndex][stateIndex];
        }

        @Override
        public double logBackwardProbability(final int positionIndex, S state) {
            ParamUtils.inRange(positionIndexRange, positionIndex, "position index");
            final int stateIndex = validStateIndex(state);
            return logBackwardProbabilities[positionIndex][stateIndex];
        }

        @Override
        public double logBackwardProbability(final T position, final S state) {
            final int stateIndex = validStateIndex(state);
            final int positionIndex = validPositionIndex(position);
            return logBackwardProbabilities[positionIndex][stateIndex];
        }

        @Override
        public double logProbability(final int positionIndex, final S state) {
            final int stateIndex = validStateIndex(state);
            ParamUtils.inRange(positionIndexRange, positionIndex, "position index");
            return logBackwardProbabilities[positionIndex][stateIndex]
                    + logForwardProbabilities[positionIndex][stateIndex]
                    - logDataLikelihood[positionIndex];
        }

        @Override
        public double logProbability(final T position, final S state) {
            final int stateIndex = this.stateIndex.getOrDefault(state, -1);
            final int positionIndex = this.positionIndex.getOrDefault(position, -1);
            Utils.validateArg(stateIndex != -1, "the input state is not recognized by the model");
            Utils.validateArg(positionIndex != -1, "unknown input position");
            return logBackwardProbabilities[positionIndex][stateIndex]
                        + logForwardProbabilities[positionIndex][stateIndex] - logDataLikelihood[positionIndex];
        }

        @Override
        public double logProbability(final List<S> states) {
            Utils.nonNull(states);
            Utils.validateArg(states.size() == data.size(), "the input states sequence does not have the same length as the data sequence");
            return states.isEmpty() ? 0 : logProbability(0, states);
        }

        @Override
        public double logProbability(final int startIndex, final List<S> states) {
            ParamUtils.inRange(positionIndexRange, startIndex, "position index");
            Utils.nonNull(states, "the input states sequence cannot be null");
            if (states.isEmpty()) {
                return 0;
            } else {
                final int statesLength = states.size();
                final int lastIndex = statesLength + startIndex - 1;
                Utils.validateArg(lastIndex < data.size(), "the input state sequence is too long");
                double result = logForwardProbability(positions.get(startIndex), states.get(0));
                for (int statesOffset = 1, dataOffset = startIndex + 1; dataOffset <= lastIndex; statesOffset++, dataOffset++) {
                    result += model.logTransitionProbability(states.get(statesOffset - 1), positions.get(dataOffset - 1), states.get(statesOffset), positions.get(dataOffset));
                    result += model.logEmissionProbability(data.get(dataOffset), states.get(statesOffset), positions.get(dataOffset));
                }
                result += logBackwardProbability(positions.get(lastIndex), states.get(statesLength - 1));
                result -= logDataLikelihood[lastIndex];
                return result;
            }
        }

        @Override
        public double logProbability(final T position, final List<S> states) {
            return logProbability(validPositionIndex(position), states);
        }

        @Override
        public double logConstrainedProbability(final int startIndex, final List<Set<S>> stateConstraints) {
            ParamUtils.inRange(positionIndexRange, startIndex, "position index");
            Utils.nonNull(stateConstraints, "the input state constraints sequence cannot be null");
            if (stateConstraints.isEmpty()) {
                return 0;
            } else {
                final int length = stateConstraints.size();
                final int lastIndex = length + startIndex - 1;
                Utils.validateArg(lastIndex < data.size(), "the input state sequence is too long");

                // calculate the likelihoods of each state at the first position using the forward-probabilities.
                List<S>  currentStates = new ArrayList<>(Utils.nonNull(stateConstraints.get(0)));
                double[] currentLikelihoods = currentStates.stream()
                        .mapToInt(stateIndex::getInt)
                        .mapToDouble(i -> logForwardProbabilities[startIndex][i])
                        .toArray();
                // We move forward across contiguous positions updating the current state likelihoods
                // with the previous ones honoring transition and emission probabilities:
                for (int statesOffset = 1, dataOffset = startIndex + 1; dataOffset <= lastIndex; statesOffset++, dataOffset++) {
                    final double[] previousLikelihoods = currentLikelihoods;
                    final List<S> previousStates = currentStates;
                    final T previousPosition = positions.get(dataOffset - 1);
                    final T thisPosition = positions.get(dataOffset);
                    final D thisData = data.get(dataOffset);

                    currentStates = new ArrayList<>(Utils.nonNull(stateConstraints.get(statesOffset)));
                    // The update likelihoods are calculated per current state
                    // by adding the likelihoods of previous states multiply by the transition
                    // probabilities from those state to the current one and then multiplying it
                    // by the emission probability of the current datum.
                    currentLikelihoods = currentStates.stream()
                            .mapToDouble(thisState ->
                                    GATKProtectedMathUtils.logSumExp(IntStream.range(0, previousStates.size())
                                            .mapToDouble(previousStateIndex -> {
                                                final S previousState = previousStates.get(previousStateIndex);
                                                return previousLikelihoods[previousStateIndex]
                                                        + model.logTransitionProbability(previousState, previousPosition,
                                                        thisState, thisPosition);
                                            }).toArray())
                                            + model.logEmissionProbability(thisData, thisState, thisPosition)
                            ).toArray();
                }

                // finally we add the backward-probabilities at the last position.
                final List<S> lastStates = currentStates;
                for (int i = 0; i < currentLikelihoods.length; i++) {
                    currentLikelihoods[i] += logBackwardProbabilities[lastIndex][stateIndex.getInt(lastStates.get(i))];
                }
                return GATKProtectedMathUtils.logSumExp(currentLikelihoods) - logDataLikelihood[lastIndex];
            }
        }

        @Override
        public double logConstrainedProbability(final T position, final List<Set<S>> stateConstraints) {
            return logConstrainedProbability(validPositionIndex(position), stateConstraints);
        }

        @Override
        public double logDataLikelihood() {
            return logDataLikelihood.length == 0 ? 0 : logDataLikelihood[0];
        }

        @Override
        public double logDataLikelihood(final int positionIndex) {
            ParamUtils.inRange(positionIndexRange, positionIndex, "position index");
            return logDataLikelihood[positionIndex];
        }

        @Override
        public double logDataLikelihood(final T position) {
            return logDataLikelihood[validPositionIndex(position)];
        }

        @Override
        public double logChainPosteriorProbability() {
            final List<S> states = model.hiddenStates();
            if (positions.isEmpty() || states.isEmpty()) {
                return 0;
            }
            /* calculate the contribution of emission nodes */
            final double logEmissionPosteriorExpectation = IntStream.range(0, positions.size())
                    .mapToDouble(positionIndex -> IntStream.range(0, states.size())
                            .mapToDouble(stateIndex -> FastMath.exp(logProbability(positionIndex, states.get(stateIndex))) *
                                    model.logEmissionProbability(data.get(positionIndex), states.get(stateIndex),
                                            positions.get(positionIndex)))
                            .sum())
                    .sum();
            return logDataLikelihood() - logEmissionPosteriorExpectation;
        }

        /**
         * Translates a model hidden state to its index.
         * @param state the input state object.
         * @return 0 or greater.
         * @throws IllegalArgumentException if the input state is not recognized by the model.
         */
        private int validStateIndex(final S state) {
            final int stateIndex = this.stateIndex.getOrDefault(state, -1);
            Utils.validateArg(stateIndex != -1, "the input state is not recognized by the model");
            return stateIndex;
        }

        /**
         * Translates a sequence position into its index.
         * @param position the input state object.
         * @return 0 or greater.
         * @throws IllegalArgumentException if the input position is not part of the original FWBW algorithm input.
         */
        private int validPositionIndex(final T position) {
            final int positionIndex = this.positionIndex.getOrDefault(position, -1);
            Utils.validateArg(positionIndex != -1, "the input position is not recognized by the model");
            return positionIndex;
        }
    }

}
