package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Segmentation and parameter learning corresponding to {@link ClusteringGenomicHMM}
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ClusteringGenomicHMMSegmenter<DATA, HIDDEN> {
    protected final Logger logger = LogManager.getLogger(ClusteringGenomicHMMSegmenter.class);

    private List<HIDDEN> hiddenStateValues;
    private double memoryLength;
    private boolean parametersHaveBeenLearned = false;

    private static final int MIN_EM_ITERATIONS = 3;
    private static final int MAX_EM_ITERATIONS = 15;
    protected final List<DATA> data;

    protected final List<SimpleInterval> positions;
    private final double[] distances;   // distances[n] is the n to n+1 distance

    protected static final double NEGLIGIBLE_POSTERIOR_FOR_M_STEP = 0.001;

    private static final double MINIMUM_MEMORY_LENGTH = 1;
    private static final double MAXIMUM_MEMORY_LENGTH = 1e10;

    protected static final double RELATIVE_CONVERGENCE_THRESHOLD = 0.01;

    protected static final double RELATIVE_TOLERANCE_FOR_OPTIMIZATION = 0.001;
    protected static final double ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION = 0.001;
    protected static final int MAX_EVALUATIONS_FOR_OPTIMIZATION = 100;

    /**
     * Initialize the segmenter with everything given i.e. without default values
     */
    public ClusteringGenomicHMMSegmenter(final List<SimpleInterval> positions,
                                         final List<DATA> data,
                                         final List<HIDDEN> hiddenStateValues,
                                         final double memoryLength) {
        this.data = Utils.nonEmpty(data);
        this.positions = Utils.nonEmpty(positions);
        Utils.validateArg(data.size() == positions.size(), "The number of data must equal the number of positions.");
        distances = calculateDistances(positions);
        this.hiddenStateValues = Utils.nonEmpty(hiddenStateValues);
        this.memoryLength = memoryLength;
    }

    private static double[] calculateDistances(final List<SimpleInterval> positions) {
        return IntStream.range(0, positions.size() - 1)
                .mapToDouble(n -> ClusteringGenomicHMM.calculateDistance(positions.get(n), positions.get(n + 1)))
                .toArray();
    }


    /**
     * given current values of memory length, weights, hidden state values, and any other parameters that a child class
     * may have, generate the model.  This is needed for running the Viterbi and forward-backward algorithms.
     */
    protected abstract ClusteringGenomicHMM<DATA, HIDDEN> makeModel();

    public void makeSureParametersHaveBeenLearned() {
        if (!parametersHaveBeenLearned) {
            learn();
        }
    }

    public List<Pair<SimpleInterval, HIDDEN>> findSegments() {
        makeSureParametersHaveBeenLearned();
        return findSegmentsDuringLearning();
    }

    private List<Pair<SimpleInterval, HIDDEN>> findSegmentsDuringLearning() {
        final ClusteringGenomicHMM<DATA, HIDDEN> model = makeModel();
        List<Integer> states = ViterbiAlgorithm.apply(data, positions, model);
        List<Pair<SimpleInterval, HIDDEN>> result = new ArrayList<>();

        int beginningOfCurrentSegment = 0;
        int currentState = states.get(0);
        String currentContig = positions.get(0).getContig();
        for (int n = 0; n <= positions.size(); n++) {
            //if contig or state has switched, make previous segment and start new one
            if (n == positions.size() || currentState != states.get(n) || !currentContig.equals(positions.get(n).getContig())) {
                final int previousSegmentStart = positions.get(beginningOfCurrentSegment).getStart();
                final int previousSegmentEnd = positions.get(n-1).getEnd();
                final HIDDEN previousSegmentHiddenStateValue = hiddenStateValues.get(currentState);
                final SimpleInterval interval = new SimpleInterval(currentContig, previousSegmentStart, previousSegmentEnd);
                result.add(new ImmutablePair<>(interval, previousSegmentHiddenStateValue));
                if (n < positions.size()) {
                    currentState = states.get(n);
                    currentContig = positions.get(n).getContig();
                    beginningOfCurrentSegment = n;
                }
            }
        }
        logger.info(String.format("Number of segments: %d", result.size()));
        return result;
    }

    private void learn() {
        int iteration = 0;
        boolean converged = false;
        List<Pair<SimpleInterval, HIDDEN>> oldSegments = new ArrayList<>();
        while (!converged && iteration++ < MAX_EM_ITERATIONS) {
            logger.info(String.format("Beginning iteration %d of learning.", iteration));
            logger.info(String.format("Current memory length: %f bases.", memoryLength));

            final List<HIDDEN> oldHiddenStateValues = new ArrayList<>(hiddenStateValues);
            performEMIteration();
            final List<Pair<SimpleInterval, HIDDEN>> currentSegments = findSegmentsDuringLearning();
            converged =
                    iteration > MIN_EM_ITERATIONS &&
                    oldSegments.size() == currentSegments.size() &&
                    hiddenStateValuesHaveConverged(oldHiddenStateValues);
            oldSegments = currentSegments;
            if (converged) {
                logger.info("Segmentation converged.");
            }
        }
        parametersHaveBeenLearned = true;
    }

    protected abstract boolean hiddenStateValuesHaveConverged(final List<HIDDEN> oldHiddenStateValues);

    // update the model and the concentration parameter with a single EM step
    private void performEMIteration() {
        final ExpectationStep expectationStep = new ExpectationStep();
        relearnHiddenStateValues(expectationStep);
        reportStates("After M step");
        relearnMemoryLength(expectationStep);
        relearnAdditionalParameters(expectationStep);
        pruneStatesByCount(expectationStep);
    }

    private void pruneStatesByCount(final ExpectationStep eStep) {
        //TODO: magic constants
        final List<Integer> indicesToKeep = IntStream.range(0, numStates())
                .filter(state -> eStep.totalCountsForState(state) > 0)
                .filter(state -> IntStream.range(0, numPositions()).mapToDouble(p -> eStep.pStateAtPosition(state,p)).max().getAsDouble() > 0.25)
                .boxed().collect(Collectors.toList());
        hiddenStateValues = indicesToKeep.stream().map(hiddenStateValues::get).collect(Collectors.toList());
    }

    private void reportStates(final String timing) {
        logger.info(timing + ", there are " + numStates() + " hidden states.");
        final StringBuilder message = new StringBuilder("States are: ");
        for (int n = 0; n < numStates(); n++) {
            message.append(String.format("%s", getState(n).toString()) + ((n < numStates() - 1) ? "; " : "."));
        }
        logger.info(message);
    }

    protected abstract void relearnAdditionalParameters(final ExpectationStep eStep);


    private void relearnMemoryLength(final ExpectationStep eStep) {
        final Function<Double, Double> objective = D -> IntStream.range(0, distances.length)
                .mapToDouble(n -> eStep.pForget(n)*Math.log(1 - Math.exp(-distances[n]/D)) - (1 - eStep.pForget(n))*(distances[n]/D))
                .sum();
        memoryLength = OptimizationUtils.argmax(objective, MINIMUM_MEMORY_LENGTH, MAXIMUM_MEMORY_LENGTH, memoryLength,
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);
    }

    protected abstract void relearnHiddenStateValues(final ExpectationStep eStep);

    /**
     * Stores the results of the expectation (E) step in which we run the forward-backward algorithm
     * to obtain posterior probabilities of 1) each hidden state at each het position; 2) each possible pair
     * of hidden states at consecutive het positions; 3) hidden state memory being lost at each het position;
     * and the expected number of 4) transitions to each state with memory loss, totalled over all het positions
     */
    protected final class ExpectationStep {
        private final int K = numStates();
        private final int N = positions.size();

        private final double[][] pStateByPosition = new double[K][N];

        // probability that memory was lost in the n to n+1 transition
        private final double[] pForget = new double[N-1];

        // transition pseudocounts for each hidden state

        public ExpectationStep() {
            final ForwardBackwardAlgorithm.Result<DATA, SimpleInterval, Integer> fbResult =
                    ForwardBackwardAlgorithm.apply(data, positions, makeModel());

            IntStream.range(0, K).forEach(state ->
                    pStateByPosition[state] = IntStream.range(0, N).mapToDouble(n -> Math.exp(fbResult.logProbability(n, state))).toArray());

            for (int n = 0; n < N-1; n++) {
                for (int from = 0; from < K; from++) {
                    for (int to = 0; to < K; to++) {
                        // probability that from -> to transition occurred going from position n to n + 1
                        final double pTransition = Math.exp(fbResult.logProbability(n, Arrays.asList(from, to)));
                        if (to != from ) {
                            pForget[n] += pTransition;
                        } else {
                            final double priorPForget = 1 - Math.exp(-distances[n] / memoryLength);
                            // Bayes' Rule gives the probability that the state was forgotten given that toState == fromState
                            final double pForgetAndTransition = pTransition * (priorPForget / K / ((1-priorPForget) + priorPForget / K));
                            pForget[n] += pForgetAndTransition;
                        }
                    }
                }
            }
        }

        public double pForget(final int position) { return pForget[position]; }
        public double pStateAtPosition(final int state, final int position) { return pStateByPosition[state][position]; }
        public double totalCountsForState(final int state) { return MathUtils.sum(pStateByPosition[state]); }
        public int numPositions() { return N; }
    }

    protected ExpectationStep getExpectationStep() { return new ExpectationStep(); }

    public int numStates() { return hiddenStateValues.size(); }
    public int numPositions() { return positions.size(); }
    public SimpleInterval getPosition(final int n) { return positions.get(n); }
    public DATA getDatum(final int n) { return data.get(n); }
    public double getMemoryLength() { return memoryLength; }
    protected HIDDEN getState(final int n) { return hiddenStateValues.get(n); }
    protected void setState(final int n, final HIDDEN value) { hiddenStateValues.set(n, value); }
    protected List<HIDDEN> getStates() { return Collections.unmodifiableList(hiddenStateValues); }
}
