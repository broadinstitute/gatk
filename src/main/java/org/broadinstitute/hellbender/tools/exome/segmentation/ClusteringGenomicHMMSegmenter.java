package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.util.*;
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

    private double concentration;
    private List<Double> weights; //one per hidden state
    private List<HIDDEN> hiddenStateValues;
    private double memoryLength;
    private boolean parametersHaveBeenLearned = false;

    private static final int MAX_EM_ITERATIONS = 50;
    protected final List<DATA> data;

    protected final List<SimpleInterval> positions;
    private final double[] distances;   // distances[n] is the n to n+1 distance

    protected static double NEGLIGIBLE_POSTERIOR_FOR_M_STEP = 0.001;

    private static final double MINIMUM_MEMORY_LENGTH = 1;
    private static final double MAXIMUM_MEMORY_LENGTH = 1e10;

    // parameters for pruning unused hidden states
    final double DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_SPURIOUS = 0.02;  // if a hidden state is this close to another state, it might be false
    final double DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_DEFINITELY_SPURIOUS = 0.01;  // if a hidden state is this close to another state, one is assumed a clone
    final double MAX_WEIGHT_CONSIDERED_FOR_PRUNING = 0.04;  // never prune a state with greater weight than this
    final double AUTOMATICALLY_PRUNED_WEIGHT = 5e-4;    // a weight so low it is always pruned

    // (unnormalized) vague gamma prior on concentration
    private static final Function<Double, Double> PRIOR_ON_CONCENTRATION = alpha -> alpha*Math.exp(-alpha);

    // a concentration parameter smaller than this would suggest less than one hidden state
    private static final double MINIMUM_CONCENTRATION = 1e-4;

    // finite cutoff for numerical integrals.  Concentration higher than this is basically impossible
    private static final double MAXIMUM_CONCENTRATION = 5;
    private static final int MAX_INTEGRATION_EVALUATIONS = 1000;
    private static final UnivariateIntegrator UNIVARIATE_INTEGRATOR = new SimpsonIntegrator(1e-3, 1e-3, 5, 20);

    protected static final double CONVERGENCE_THRESHOLD = 0.01;
    private static final double MEMORY_LENGTH_CONVERGENCE_THRESHOLD = 1e4;

    protected static final double RELATIVE_TOLERANCE_FOR_OPTIMIZATION = 0.01;
    protected static final double ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION = 0.01;
    protected static final int MAX_EVALUATIONS_FOR_OPTIMIZATION = 100;

    //RNG used in attemptBigChangeInMemoryLength method
    private static final int RANDOM_SEED = 239;
    private final Random random = new Random(RANDOM_SEED);


    /**
     * Initialize the segmenter with everything given i.e. without default values
     */
    public ClusteringGenomicHMMSegmenter(final List<SimpleInterval> positions,
                                         final List<DATA> data,
                                         final List<HIDDEN> hiddenStateValues,
                                         final List<Double> weights,
                                         final double concentration,
                                         final double memoryLength) {
        this.data = Utils.nonEmpty(data);
        this.positions = Utils.nonEmpty(positions);
        Utils.validateArg(data.size() == positions.size(), "The number of data must equal the number of positions.");
        distances = calculateDistances(positions);
        this.concentration = concentration;
        this.hiddenStateValues = Utils.nonEmpty(hiddenStateValues);
        this.weights = Utils.nonEmpty(weights);
        Utils.validateArg(hiddenStateValues.size() == weights.size(), "The number of hidden states must equal the number of weights.");
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
        return result;
    }

    private void learn() {
        int iteration = 0;
        boolean converged = false;
        while (!converged && iteration++ < MAX_EM_ITERATIONS) {
            logger.info(String.format("Beginning iteration %d of learning.", iteration));
            logger.info(String.format("Current memory length: %f bases.", memoryLength));

            final double oldMemoryLength = memoryLength;
            final List<Double> oldWeights = new ArrayList<>(weights);
            final List<HIDDEN> oldHiddenStateValues = new ArrayList<>(hiddenStateValues);
            performEMIteration();
            converged = oldWeights.size() == numStates() &&
                    Math.abs(oldMemoryLength - memoryLength) < MEMORY_LENGTH_CONVERGENCE_THRESHOLD &&
                    GATKProtectedMathUtils.maxDifference(oldWeights, weights) < CONVERGENCE_THRESHOLD &&
                    hiddenStateValuesHaveConverged(oldHiddenStateValues);
        }
        parametersHaveBeenLearned = true;
    }

    protected abstract boolean hiddenStateValuesHaveConverged(final List<HIDDEN> oldHiddenStateValues);

    // update the model and the concentration parameter with a single EM step
    private void performEMIteration() {
        final ExpectationStep expectationStep = new ExpectationStep();
        relearnHiddenStateValues(expectationStep);
        relearnWeights(expectationStep);
        reportStatesAndWeights("After M step");
        relearnMemoryLength(expectationStep);
        attemptBigChangeInMemoryLength();
        relearnAdditionalParameters(expectationStep);
        pruneUnusedComponents();
        reportStatesAndWeights("After pruning");
        relearnConcentration();
    }

    private void reportStatesAndWeights(final String timing) {
        logger.info(timing + ", there are " + numStates() + " hidden states.");
        final StringBuilder message = new StringBuilder("(state, weight) pairs are: ");
        for (int n = 0; n < numStates(); n++) {
            message.append(String.format("(%s, %f)", getState(n).toString(), weights.get(n)) + ((n < numStates() - 1) ? "; " : "."));
        }
        logger.info(message);
    }

    protected abstract void relearnAdditionalParameters(final ExpectationStep eStep);

    private void relearnWeights(final ExpectationStep expectationStep) {
        final double symmetricPriorWeight = concentration / numStates();
        final double[] transitionCounts = expectationStep.transitionCounts();
        final double[] posteriorDirichletParameters = Arrays.stream(transitionCounts).map(x -> x + symmetricPriorWeight).toArray();
        weights = Doubles.asList(new Dirichlet(posteriorDirichletParameters).effectiveMultinomialWeights());
    }

    protected abstract void pruneUnusedComponents();

    protected void removeStates(Collection<Integer> componentsToPrune) {
        weights = IntStream.range(0, numStates())
                .filter(n -> !componentsToPrune.contains(n)).mapToObj(weights::get).collect(Collectors.toList());
        hiddenStateValues = IntStream.range(0, numStates())
                .filter(n -> !componentsToPrune.contains(n)).mapToObj(hiddenStateValues::get).collect(Collectors.toList());
    }

    /**
     * Compute the effective value of the Dirichlet concentration parameter, which defines the prior on weights in
     * subsequent iterations.  This value is the expectation of the concentration with respect to it mean-field
     * variational Bayes posterior distribution.  This mean-field comprises the prior on concentration, the
     * concentration-dependent Dirichlet distribution normalization constant, and the Dirichlet likelihood of the
     * effective weights.
     */
    private void relearnConcentration() {
        final double geometricMeanOfEffectiveWeights = Math.exp(weights.stream().mapToDouble(Math::log).average().getAsDouble());

        final int K = numStates();
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
        final double memoryLengthMultiplier = Math.exp(random.nextGaussian()/2);
        final double currentLogLikelihood = ForwardBackwardAlgorithm.apply(data, positions, makeModel()).logDataLikelihood();
        final double oldMemoryLength = memoryLength;
        memoryLength = memoryLength * memoryLengthMultiplier;
        final double proposalLogLikelihood = ForwardBackwardAlgorithm.apply(data, positions, makeModel()).logDataLikelihood();
        if (proposalLogLikelihood < currentLogLikelihood) {
            memoryLength = oldMemoryLength;
        }
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
        private final double[] transitionCountsByState = new double[K];

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
                            transitionCountsByState[to] += pTransition;
                            pForget[n] += pTransition;
                        } else {
                            final double priorPForget = 1 - Math.exp(-distances[n] / memoryLength);
                            // Bayes' Rule gives the probability that the state was forgotten given that toState == fromState
                            final double pForgetAndTransition = pTransition * (priorPForget * getWeight(to) / ((1-priorPForget) + priorPForget * getWeight(to)));
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
        public int numPositions() { return N; }
    }

    protected ExpectationStep getExpectationStep() { return new ExpectationStep(); }

    public int numStates() { return hiddenStateValues.size(); }
    public int numPositions() { return positions.size(); }
    public SimpleInterval getPosition(final int n) { return positions.get(n); }
    public DATA getDatum(final int n) { return data.get(n); }
    public double getConcentration() { return concentration; }
    public double getMemoryLength() { return memoryLength; }
    protected double getWeight(final int n) { return weights.get(n); }
    protected HIDDEN getState(final int n) { return hiddenStateValues.get(n); }
    protected void setState(final int n, final HIDDEN value) { hiddenStateValues.set(n, value); }
    protected List<Double> getWeights() { return Collections.unmodifiableList(weights); }
    protected List<HIDDEN> getStates() { return Collections.unmodifiableList(hiddenStateValues); }
}
