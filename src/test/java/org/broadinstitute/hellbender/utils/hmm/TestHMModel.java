package org.broadinstitute.hellbender.utils.hmm;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

/**
 * Simple HMM for tests only.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class TestHMModel implements HiddenMarkovModel<TestHMModel.Datum, Integer, TestHMModel.State> {

    private final double[] priors;

    private final double[][] transitionProbabilities;

    private final double[][] emissionProbabilities;


    private final static double LN_OF_10 = Math.log(10);

    static TestHMModel fromPhredProbabilities(final double ... phred) {

        final double[] priors = Arrays.copyOfRange(phred, 0, State.values().length);
        final double[][] transition = reshape(
                Arrays.copyOfRange(phred, State.values().length, State.values().length * (State.values().length + 1)), State.values().length);
        final double[][] emission = reshape(
                Arrays.copyOfRange(phred, (State.values().length + 1) * State.values().length, phred.length), State.values().length);

        final DoubleUnaryOperator phredToLog = d -> d * -.1 * LN_OF_10;
        final double[] logPriorsRaw = DoubleStream.of(priors).map(phredToLog).toArray();
        final double logPriorsSum = GATKProtectedMathUtils.naturalLogSumExp(logPriorsRaw);
        final double[] logPriors = DoubleStream.of(logPriorsRaw).map(d -> d - logPriorsSum).toArray();

        final double[][] logEmissionProbs = Stream.of(emission)
                .map(x -> { final double[] result = DoubleStream.of(x).map(phredToLog).toArray();
                            final double sum = GATKProtectedMathUtils.naturalLogSumExp(result);
                            return DoubleStream.of(result).map(d -> d - sum).toArray(); })
                .toArray(double[][]::new);
        final double[][] logTransitionProbs = Stream.of(transition)
                .map(x -> { final double[] result = DoubleStream.of(x).map(phredToLog).toArray();
                            final double sum = GATKProtectedMathUtils.naturalLogSumExp(result);
                            return DoubleStream.of(result).map(d -> d - sum).toArray(); })
                .toArray(double[][]::new);
        return new TestHMModel(logPriors, logEmissionProbs, logTransitionProbs);

    }

    private static double[][] reshape(final double[] values, final int rows) {
        final int cols = values.length / rows;
        final double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = values[j + rows * i];
            }
        }
        return result;
    }

    /**
     * Creates a new model.
     *
     * @param priors the state priors.
     * @param emission the data emission.
     * @param transition transition probabilities
     */
    private TestHMModel(final double[] priors, final double[][] emission, final double[][] transition) {
        this.priors = priors;
        this.emissionProbabilities = emission;
        this.transitionProbabilities = transition;
    }

    @Override
    public List<State> hiddenStates() {
        return Collections.unmodifiableList(Arrays.asList(State.A, State.B, State.C));
    }

    @Override
    public double logPriorProbability(final State state, final Integer position) {
        Utils.nonNull(state);
        return priors[state.ordinal()];
    }

    @Override
    public double logTransitionProbability(final State currentState, final Integer currentPosition, final State nextState, final Integer nextPosition) {
        Utils.nonNull(currentState);
        Utils.nonNull(nextState);
        return transitionProbabilities[currentState.ordinal()][nextState.ordinal()];
    }

    @Override
    public double logEmissionProbability(final Datum datum, final State state, final Integer position) {
        Utils.nonNull(state);
        Utils.nonNull(datum);
        return emissionProbabilities[state.ordinal()][datum.ordinal()];
    }

    enum State {
        A, B, C;
    }

    enum Datum {
        X, Y, Z;
    }

    public Pair<List<State>, List<Datum>> generate(final List<Integer> positions, final Random rdn) {
        Utils.nonNull(positions);
        Utils.nonNull(rdn);
        final List<Datum> data = new ArrayList<>(positions.size());
        final List<State> states = new ArrayList<>(positions.size());
        final Integer firstPosition = positions.get(0);
        final State firstState = generateFirstState(firstPosition, rdn);
        states.add(firstState);
        data.add(generateDatum(firstState, firstPosition, rdn));
        State previousState = firstState;
        Integer previousPosition = firstPosition;
        for (int i = 1; i < positions.size(); i++) {
            final Integer position = positions.get(i);
            final State nextState = generateNextState(previousState, previousPosition, position, rdn);
            data.add(generateDatum(nextState, position, rdn));
            states.add(previousState = nextState);
            previousPosition = position;
        }
        return new Pair<>(states, data);
    }

    static String toHMMInstallRString() {
        return "library(HMM)";
    }

    static String toHMMModelDeclarationRString(final List<TestHMModel> models) {
        return String.format("list(\n\t%s\n)", models.stream().map(TestHMModel::toRExpressionString)
              .collect(Collectors.joining(",\n\t")));
    };

    /**
     * Composes the R code to create an equivalent model
     * @return never {@code null}.
     */
    public String toRExpressionString() {
        final String statesString = String.format("c(%s)",
                Stream.of(State.values())
                        .map(s -> "\"" + s + "\"").collect(Collectors.joining(", ")));

        final String symbolsString = String.format("c(%s)",
                Stream.of(Datum.values())
                        .map(s -> "\"" + s + "\"").collect(Collectors.joining(", ")));

        final String startProbsString = String.format("c(%s)", DoubleStream.of(priors)
                .map(d -> Math.exp(d))
                .mapToObj(d -> String.format("%.4g", d))
                .collect(Collectors.joining(", ")));

        final String transProbsString = String.format("matrix(c(%s), nrow=%d, byrow=T)",
                DoubleStream.of(Stream.of(this.transitionProbabilities)
                        .reduce(ArrayUtils::addAll).get())
                        .map(d -> Math.exp(d))
                        .mapToObj(d -> String.format("%.4g", d))
                        .collect(Collectors.joining(", "))
                , State.values().length);

        final String emissionProbsString = String.format("matrix(c(%s), nrow=%d, byrow=T)",
                DoubleStream.of(Stream.of(this.emissionProbabilities)
                        .reduce(ArrayUtils::addAll).get())
                        .map(d -> Math.exp(d))
                        .mapToObj(d -> String.format("%.4g", d))
                        .collect(Collectors.joining(", "))
                , State.values().length);

        final String result = String.format("initHMM(States=%s, Symbols=%s, startProbs=%s, transProbs=%s, emissionProbs=%s)",
                statesString, symbolsString, startProbsString, transProbsString, emissionProbsString);
        return result;
    }

    private State generateNextState(final State previousState, final Integer previousPosition, final Integer nextPosition, final Random rdn) {
        final double[] logTransitionProbs = new double[State.values().length];
        for (int i = 0; i < logTransitionProbs.length; i++) {
            logTransitionProbs[i] = logTransitionProbability(previousState, previousPosition, State.values()[i], nextPosition);
        }
        final double total = Math.exp(GATKProtectedMathUtils.naturalLogSumExp(logTransitionProbs));
        final double uniform = rdn.nextDouble() * total;
        double accumulative = 0;
        for (int i = 0; i < logTransitionProbs.length; i++) {
            accumulative += Math.exp(logTransitionProbs[i]);
            if (uniform < accumulative) {
                return State.values()[i];
            }
        }
        throw new IllegalStateException("may be some of the trans probs are negative");
    }

    private Datum generateDatum(final State state, final Integer position, final Random rdn) {
        final double[] logEmissionProbs = new double[State.values().length];
        for (int i = 0; i < logEmissionProbs.length; i++) {
            logEmissionProbs[i] = logEmissionProbability(Datum.values()[i], state, position);
        }
        final double total = Math.exp(GATKProtectedMathUtils.naturalLogSumExp(logEmissionProbs));

        final double uniform = rdn.nextDouble() * total;
        double accumulate = 0;
        for (int i = 0; i < Datum.values().length; i++) {
            accumulate += Math.exp(logEmissionProbs[i]);
            if (uniform < accumulate) {
                return Datum.values()[i];
            }
        }
        throw new IllegalStateException("may be some of the emission probs are negative?");
    }

    private State generateFirstState(final Integer firstPosition, final Random rdn) {
        final double uniform = rdn.nextDouble() * Math.exp(GATKProtectedMathUtils.naturalLogSumExp(priors));
        double accumulative = 0;
        for (int i = 0; i < priors.length; i++) {
            accumulative += Math.exp(logPriorProbability(State.values()[i], firstPosition));
            if (uniform < accumulative) {
                return State.values()[i];
            }
        }
        throw new IllegalStateException("bug in testing code: seems that the test model priors sum to more than 1: " + Arrays.toString(priors));
    }
}
