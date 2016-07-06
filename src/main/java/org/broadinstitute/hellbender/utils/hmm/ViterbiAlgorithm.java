package org.broadinstitute.hellbender.utils.hmm;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Implements the Viterbi Algorithm.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ViterbiAlgorithm {

    /**
     * Calculates the most likely the hidden state sequence that explains the observed data
     * given the corresponding sequence of observation positions.
     *
     * @param data the observed data.
     * @param positions the observation positions.
     * @param model the HMM model.
     * @param <D> observed data type.
     * @param <T> time data type.
     * @param <S> hidden state data-type.
     * @return never {@code null}, with only valid states as returned by <code>model.
     * {@link HiddenMarkovModel#hiddenStates() hiddenStates()}</code>.
     *         its length will be the sames as the length of the input {@code data}.
     * @throws IllegalArgumentException if any of the following is true:
     * <ul>
     *     <li>any {@code data}, {@code positions} or {@code model} is {@code null},</li>
     *     <li>{@code data}, {@code positions} contains elements that are not expected by the implementation
     *     of {@code model},</li>
     *     <li>{@code data} and {@code positions} have different lengths.</li>
     * </ul>
     */
    public static <D, T, S> List<S> apply(final List<D> data, final List<T> positions,
                                                          final HiddenMarkovModel<D, T, S> model) {
        checkApplyArguments(data, positions, model);

        if (data.isEmpty()) {
            return new ArrayList<>(0);
        }

        @SuppressWarnings("unchecked")
        final S[] states = (S[]) model.hiddenStates().stream().toArray(Object[]::new);
        final int length = data.size();
        final int numStates = states.length;

        // We alternate between the two arrays as we move along time.
        // These two arrays will contain the best paths of two contiguous observation positions.
        // t and t + 1.
        // The array index for time x is determined by the lowest bit in x.
        @SuppressWarnings({"rawtypes", "unchecked"})
        final Path<S>[][] bestPaths = (Path<S>[][]) new Path[][] {
            initializeBestPaths(data, positions, model, states), // best paths at time 0.
            (Path<S>[]) new Path[numStates]                      // empty best paths for time 1.
        };

        for (int thisPositionIndex = 1; thisPositionIndex < length; thisPositionIndex++) {
            final int previousPositionIndex = thisPositionIndex - 1;
            final Path<S>[] bestCurrentPaths = bestPaths[1 & thisPositionIndex]; // time thisPositionIndex best path; where the new best paths go.
            final Path<S>[] bestPreviousPaths = bestPaths[1 & previousPositionIndex]; // time previousPositionIndex best paths; previous best paths.
            final T previousPosition = positions.get(previousPositionIndex);
            final D thisDatum = data.get(thisPositionIndex);
            final T thisPosition = positions.get(thisPositionIndex);
            for (int thisStateIndex = 0; thisStateIndex < numStates; thisStateIndex++) {
                final S thisState = states[thisStateIndex];
                // Initialize best-previousPositionIndex-state search setting the best so
                // far to 0th indexed state:
                int bestPreviousStateIndex = 0;
                double bestPreviousStateLogProb = bestPreviousPaths[0].logProb
                        + model.logTransitionProbability(states[0], previousPosition, thisState, thisPosition);
                // Then we check on the 1th state, the 2nd state and so forth:
                for (int previousStateIndex = 1; previousStateIndex < numStates; previousStateIndex++) {
                    final double candidatePreviousStateLogProb = bestPreviousPaths[previousStateIndex].logProb
                            + model.logTransitionProbability(states[previousStateIndex], previousPosition, thisState, thisPosition);
                    if (candidatePreviousStateLogProb > bestPreviousStateLogProb) {
                        bestPreviousStateLogProb = candidatePreviousStateLogProb;
                        bestPreviousStateIndex = previousStateIndex;
                    }
                }
                bestCurrentPaths[thisStateIndex] = bestPreviousPaths[bestPreviousStateIndex].extend(thisState,
                        bestPreviousStateLogProb + model.logEmissionProbability(thisDatum, thisState, thisPosition));
            }
        }

        return composeBestStateSequence(length, bestPaths[(length - 1) & 1]);
    }

    private static <S> List<S> composeBestStateSequence(final int length, final Path<S>[] bestPaths) {

        // Get the best path amongst the best path that finish in every state.
        Path<S> bestPath = Stream.of(bestPaths)
                .sorted(Comparator.comparing(p -> -p.logProb))
                .findFirst()
                .get();

        // Fill out the array backwards.
        @SuppressWarnings("unchecked")
        final S[] result = (S[]) new Object[length];
        for (int i = result.length - 1; i >= 0; --i) {
            result[i] = bestPath.state;
            bestPath = bestPath.prefix;
        }

        // Return the proper data-type: a modifiable List<S>.
        return Stream.of(result).collect(Collectors.toList());
    }


    @SuppressWarnings("unchecked")
    private static <D, T, S> Path<S>[] initializeBestPaths(
            final List<D> data,
            final List<T> positions,
            final HiddenMarkovModel<D, T, S> model, S[] states) {
        final D data0 = data.get(0);
        final T position0 = positions.get(0);

        // This array's position i, contains the most likely paths up to the ending in the ith state at the
        // current time t.
        return (Path<S>[]) Stream.of(states)
                .map(s -> Path.makeFirst(s, model.logPriorProbability(s, position0)
                        + model.logEmissionProbability(data0, s, position0)))
                .toArray(Path[]::new);
    }

    private static <D, T, S> void checkApplyArguments(List<D> data, List<T> times, HiddenMarkovModel<D, T, S> model) {
        Utils.nonNull(data);
        Utils.nonNull(times);
        Utils.nonNull(model);
        if (data.size() != times.size()) {
            throw new IllegalArgumentException("the data and time input sequences must have the same length");
        }
    }

    /**
     * Internal linked representation for the best paths.
     * @param <S> the hidden-state type.
     */
    private static class Path<S> {

        /**
         * Prefix of this path (as long as this path - 1 state).
         */
        private final Path<S> prefix;

        /**
         * Last state of the path.
         */
        private final S state;

        /**
         * Probability of the path.
         */
        private final double logProb;

        @SuppressWarnings("unchecked")
        private static <S> Path<S> makeFirst(final S state, final double logProb) {
            return new Path<>(null, state, logProb);
        }

        private Path<S> extend(final S state, final double logProb) {
            return new Path<>(this, state, logProb);
        }

        private Path(final Path<S> prefix, final S state, final double logProb) {
            this.state = state;
            this.prefix = prefix;
            this.logProb = logProb;
        }
    }
}
