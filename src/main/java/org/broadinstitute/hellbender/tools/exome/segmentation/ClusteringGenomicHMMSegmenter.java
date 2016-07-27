package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.util.*;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * Segmentation and parameter learning corresponding to {@link ClusteringGenomicHMM}
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ClusteringGenomicHMMSegmenter<T> {
    protected final Logger logger = LogManager.getLogger(ClusteringGenomicHMM.class);

    private double concentration;
    protected double[] weights; //one per hidden state
    protected double[] hiddenStateValues;
    protected double memoryLength;
    private boolean parametersHaveBeenLearned = false;

    protected static final int NEUTRAL_VALUE_INDEX = 0;

    private static final int MAX_EM_ITERATIONS = 50;
    protected final List<T> data;
    protected final List<SimpleInterval> positions;
    private final double[] distances;   // distances[n] is the n to n+1 distance

    protected static double NEGLIGIBLE_POSTERIOR_FOR_M_STEP = 0.001;

    private static final double MINIMUM_MEMORY_LENGTH = 1;
    private static final double MAXIMUM_MEMORY_LENGTH = 1e10;
    private static final double DEFAULT_MEMORY_LENGTH = 5e7;

    // parameters for pruning unused hidden states
    final double DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_SPURIOUS = 0.01;  // if a hidden state is this close to another state, it might be false
    final double MAX_WEIGHT_CONSIDERED_FOR_PRUNING = 0.02;  // never prune a state with greater weight than this
    final double AUTOMATICALLY_PRUNED_WEIGHT = 5e-4;    // a weight so low it is always pruned

    // (unnormalized) vague gamma prior on concentration
    private static final Function<Double, Double> PRIOR_ON_CONCENTRATION = alpha -> alpha*Math.exp(-alpha);

    // a concentration parameter smaller than this would suggest less than one hidden state
    private static final double MINIMUM_CONCENTRATION = 1e-4;

    // finite cutoff for numerical integrals.  Concentration higher than this is basically impossible
    private static final double MAXIMUM_CONCENTRATION = 5;
    private static final int MAX_INTEGRATION_EVALUATIONS = 1000;
    private static final UnivariateIntegrator UNIVARIATE_INTEGRATOR = new SimpsonIntegrator(1e-3, 1e-3, 5, 20);

    private static final double CONVERGENCE_THRESHOLD = 0.01;
    private static final double MEMORY_LENGTH_CONVERGENCE_THRESHOLD = 1e4;

    protected static final double RELATIVE_TOLERANCE_FOR_OPTIMIZATION = 0.01;
    protected static final double ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION = 0.01;
    protected static final int MAX_EVALUATIONS_FOR_OPTIMIZATION = 100;

    /**
     * Initialize the segmenter with its data and panel of normals, giving equal weight to a set of evenly-spaced
     * hidden minor allele fraction values.
     *
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     * @param data              The data attached to this segmenter
     */
    public ClusteringGenomicHMMSegmenter(final int initialNumStates, final List<SimpleInterval> positions, final List<T> data) {
        this.data = Utils.nonNull(data);
        concentration = 1;
        weights = Collections.nCopies(initialNumStates, 1.0/initialNumStates).stream().mapToDouble(x->x).toArray();   //uniform
        memoryLength = DEFAULT_MEMORY_LENGTH;
        this.positions = positions;
        distances = IntStream.range(0, positions.size() - 1)
                .mapToDouble(n -> ClusteringGenomicHMM.calculateDistance(positions.get(n), positions.get(n + 1)))
                .toArray();
        initializeHiddenStateValues(initialNumStates);
        initializeAdditionalParameters();
    }

    protected abstract void initializeHiddenStateValues(final int K);

    // override this if any parameters other than weights, memory length, and hidden state values must be initialized
    protected abstract void initializeAdditionalParameters();

    /**
     * given current values of memory length, weights, hidden state values, and any other parameters that a child class
     * may have, generate the model.  This is needed for running the Viterbi and forward-backward algorithms.
     */
    protected abstract ClusteringGenomicHMM<T> makeModel();

    public List<ModeledSegment> findSegments() {
        if (!parametersHaveBeenLearned) {
            learn();
        }
        final ClusteringGenomicHMM<T> model = makeModel();
        List<Integer> states = ViterbiAlgorithm.apply(data, positions, model);
        List<ModeledSegment> result = new ArrayList<>();

        int beginningOfCurrentSegment = 0;
        int currentState = states.get(0);
        String currentContig = positions.get(0).getContig();
        for (int n = 0; n <= positions.size(); n++) {
            //if contig or state has switched, make previous segment and start new one
            if (n == positions.size() || currentState != states.get(n) || !currentContig.equals(positions.get(n).getContig())) {
                final int previousSegmentStart = positions.get(beginningOfCurrentSegment).getStart();
                final int previousSegmentEnd = positions.get(n-1).getEnd();
                final int previousSegmentSnpCount = n - beginningOfCurrentSegment;
                final double previousSegmentHiddenStateValue = hiddenStateValues[currentState];
                final SimpleInterval interval = new SimpleInterval(currentContig, previousSegmentStart, previousSegmentEnd);
                result.add(new ModeledSegment(interval, previousSegmentSnpCount, previousSegmentHiddenStateValue));
                if (n < positions.size()) {
                    currentState = states.get(n);
                    currentContig = positions.get(n).getContig();
                    beginningOfCurrentSegment = n;
                }
            }
        }
        return result;
    }

    private void learn() {
        int iteration = 0;
        boolean converged = false;
        while (!converged && iteration++ < MAX_EM_ITERATIONS) {
            logger.info(String.format("Beginning iteration %d of learning.", iteration));
            logger.info(String.format("Current values of hidden states: %s.", Arrays.toString(hiddenStateValues)));
            logger.info(String.format("Current weights of hidden states: %s.", Arrays.toString(weights)));
            logger.info(String.format("Current memory length: %f bases.", memoryLength));

            final double oldMemoryLength = memoryLength;
            final double[] oldWeights = weights.clone();
            final double[] oldHiddenStateValues = hiddenStateValues.clone();
            performEMIteration();
            converged = oldWeights.length == weights.length &&
                    Math.abs(oldMemoryLength - memoryLength) < MEMORY_LENGTH_CONVERGENCE_THRESHOLD &&
                    GATKProtectedMathUtils.maxDifference(oldWeights, weights) < CONVERGENCE_THRESHOLD &&
                    GATKProtectedMathUtils.maxDifference(oldHiddenStateValues, hiddenStateValues) < CONVERGENCE_THRESHOLD;
        }
        parametersHaveBeenLearned = true;
    }

    // update the model and the concentration parameter with a single EM step
    private void performEMIteration() {
        logger.info("Performing E step via the forward-backward algorithm.");
        final ExpectationStep expectationStep = new ExpectationStep();
        logger.info("Performing M step optimizations of parameters.");
        logger.info("Relearning hidden state values.");
        relearnHiddenStateValues(expectationStep);
        logger.info("Relearning hidden state weights.");
        relearnWeights(expectationStep);
        logger.info("Relearning memory length.");
        relearnMemoryLength(expectationStep);
        logger.info("Attempting larger change in memory length");
        attemptBigChangeInMemoryLength();
        logger.info("Relearning additional parameters.");
        relearnAdditionalParameters(expectationStep);
        logger.info("Attempting to reduce the number of hidden states used.");
        pruneUnusedComponents();
        logger.info("Number of states after pruning: " + hiddenStateValues.length);
        logger.info("Relearning concentration hyperparameter.");
        relearnConcentration();
    }

    protected abstract void relearnAdditionalParameters(final ExpectationStep eStep);

    private void relearnWeights(final ExpectationStep expectationStep) {
        final Dirichlet priorOnWeights = Dirichlet.symmetricDirichlet(weights.length, concentration);
        weights =  new Dirichlet(priorOnWeights, expectationStep.transitionCounts()).effectiveMultinomialWeights();
    }

    // filter out components that have low weight and are too close to another component -- these will
    // die out eventually in EM, but very slowly, so we hasten their demise for quicker convergence
    private void pruneUnusedComponents() {
        final int K = weights.length;

        final Set<Integer> componentsToPrune = new TreeSet<>();
        for (final int state : IntStream.range(0, K).filter(n -> n != NEUTRAL_VALUE_INDEX).toArray()) {
            final int closestOtherState = MathUtils.maxElementIndex(IntStream.range(0, K)
                    .mapToDouble(n -> n == state ? Double.NEGATIVE_INFINITY : -Math.abs(hiddenStateValues[n] - hiddenStateValues[state]))
                    .toArray());
            final boolean hasLowWeight = weights[state] < MAX_WEIGHT_CONSIDERED_FOR_PRUNING;
            final boolean isCloseToNeighbor = Math.abs(hiddenStateValues[state] - hiddenStateValues[closestOtherState]) < DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_SPURIOUS;
            final boolean hasLessWeightThanNeighbor = weights[state] < weights[closestOtherState];
            if (weights[state] < AUTOMATICALLY_PRUNED_WEIGHT || (hasLowWeight && isCloseToNeighbor && hasLessWeightThanNeighbor)) {
                componentsToPrune.add(state);
            }
        }

        weights = IntStream.range(0, K)
                .filter(n -> !componentsToPrune.contains(n)).mapToDouble(n -> weights[n]).toArray();
        hiddenStateValues = IntStream.range(0, K)
                .filter(n -> !componentsToPrune.contains(n)).mapToDouble(n -> hiddenStateValues[n]).toArray();
    }

    /**
     * Compute the effective value of the Dirichlet concentration parameter, which defines the prior on weights in
     * subsequent iterations.  This value is the expectation of the concentration with respect to it mean-field
     * variational Bayes posterior distribution.  This mean-field comprises the prior on concentration, the
     * concentration-dependent Dirichlet distribution normalization constant, and the Dirichlet likelihood of the
     * effective weights.
     */
    private void relearnConcentration() {
        final double geometricMeanOfEffectiveWeights = Math.exp(Arrays.stream(weights).map(Math::log).average().getAsDouble());

        final int K = weights.length;
        final Function<Double, Double> distribution = alpha -> PRIOR_ON_CONCENTRATION.apply(alpha)
                * Math.pow(geometricMeanOfEffectiveWeights, alpha)  //likelihood
                * Math.exp(Gamma.logGamma(alpha) - K * Gamma.logGamma(alpha/K));    //normalization constant

        concentration = UNIVARIATE_INTEGRATOR.integrate(MAX_INTEGRATION_EVALUATIONS, alpha ->  alpha * distribution.apply(alpha), MINIMUM_CONCENTRATION, MAXIMUM_CONCENTRATION)
                / UNIVARIATE_INTEGRATOR.integrate(MAX_INTEGRATION_EVALUATIONS, distribution::apply, MINIMUM_CONCENTRATION, MAXIMUM_CONCENTRATION);
    }

    private void relearnMemoryLength(final ExpectationStep eStep) {
        final Function<Double, Double> objective = D -> IntStream.range(0, distances.length)
                .mapToDouble(n -> eStep.pForget(n)*Math.log(1 - Math.exp(-distances[n]/D)) - (1 - eStep.pForget(n))*(distances[n]/D))
                .sum();
        memoryLength = OptimizationUtils.argmax(objective, MINIMUM_MEMORY_LENGTH, MAXIMUM_MEMORY_LENGTH, memoryLength,
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);
    }

    // memoryLength may move slowly because the E and M steps are entangled.  We try to get around this by maximizing the
    // exact model log-likelihood
    private void attemptBigChangeInMemoryLength() {
        final double memoryLengthMultiplier = Math.exp(new Random().nextGaussian()/2);
        final double currentLogLikelihood = ForwardBackwardAlgorithm.apply(data, positions, makeModel()).logDataLikelihood();
        final double oldMemoryLength = memoryLength;
        memoryLength = memoryLength * memoryLengthMultiplier;
        final double proposalLogLikelihood = ForwardBackwardAlgorithm.apply(data, positions, makeModel()).logDataLikelihood();
        if (proposalLogLikelihood < currentLogLikelihood) {
            memoryLength = oldMemoryLength;
        }
    }

    private void relearnHiddenStateValues(final ExpectationStep eStep) {
        final ClusteringGenomicHMM<T> model = makeModel();
        // by convention, state = 0 represents the neutral value (minor allele fraction = 1/2 or copy ratio = 1)
        // which we always retain and do not wish to adjust via MLE.  Thus we start at state 1
        for (final int state : IntStream.range(0, hiddenStateValues.length).filter(n -> n != NEUTRAL_VALUE_INDEX).toArray()) {
            final Function<Double, Double> objective = f -> IntStream.range(0, data.size())
                    .filter(n -> eStep.pStateAtPosition(state, n) > NEGLIGIBLE_POSTERIOR_FOR_M_STEP)
                    .mapToDouble(n -> eStep.pStateAtPosition(state, n) * model.logEmissionProbability(data.get(n), f))
                    .sum();
            hiddenStateValues[state] = OptimizationUtils.singleNewtonArgmaxUpdate(objective, minHiddenStateValue(),
                    maxHiddenStateValue(), hiddenStateValues[state]);
        }
    }

    protected abstract double minHiddenStateValue();
    protected abstract double maxHiddenStateValue();

    /**
     * Stores the results of the expectation (E) step in which we run the forward-backward algorithm
     * to obtain posterior probabilities of 1) each hidden state at each het position; 2) each possible pair
     * of hidden states at consecutive het positions; 3) hidden state memory being lost at each het position;
     * and the expected number of 4) transitions to each state with memory loss, totalled over all het positions
     */
    protected final class ExpectationStep {
        private final int K = weights.length;
        private final int N = positions.size();

        private final double[][] pStateByPosition = new double[K][N];

        // probability that memory was lost in the n to n+1 transition
        private final double[] pForget = new double[N-1];

        // transition pseudocounts for each hidden state
        private final double[] transitionCountsByState = new double[K];

        public ExpectationStep() {
            final ForwardBackwardAlgorithm.Result<T, SimpleInterval, Integer> fbResult =
                    ForwardBackwardAlgorithm.apply(data, positions, makeModel());

            IntStream.range(0, K).forEach(state ->
                    pStateByPosition[state] = IntStream.range(0, N).mapToDouble(n -> Math.exp(fbResult.logProbability(n, state))).toArray());

            for (int n = 0; n < N-1; n++) {
                for (int from = 0; from < K; from++) {
                    for (int to = 0; to < K; to++) {
                        // probability that from -> to transition occurred going from position n to n + 1
                        final double pTransition = Math.exp(fbResult.logProbability(n, Arrays.asList(from, to)));
                        if (to != from ) {
                            transitionCountsByState[to] += pTransition;
                            pForget[n] += pTransition;
                        } else {
                            final double priorPForget = 1 - Math.exp(-distances[n] / memoryLength);
                            // Bayes' Rule gives the probability that the state was forgotten given that toState == fromState
                            final double pForgetAndTransition = pTransition * (priorPForget * weights[to] / ((1-priorPForget) + priorPForget * weights[to]));
                            transitionCountsByState[to] += pForgetAndTransition;
                            pForget[n] += pForgetAndTransition;
                        }
                    }
                }
            }
        }

        public double pForget(final int position) { return pForget[position]; }
        public double pStateAtPosition(final int state, final int position) { return pStateByPosition[state][position]; }
        public double[] transitionCounts() { return transitionCountsByState; }
    }
}
