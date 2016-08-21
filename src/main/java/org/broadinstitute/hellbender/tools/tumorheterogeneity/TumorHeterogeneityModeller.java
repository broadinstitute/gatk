package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.AllelicCNV;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.ModelSampler;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedModel;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedModel.EnsembleBuilder;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Performs Markov-Chain Monte Carlo inference for {@link TumorHeterogeneity} and stores the generated samples.
 * Uses affine-invariant ensemble sampling (Goodman & Weare 2010, implemented in {@link ParameterizedModel.EnsembleBuilder})
 * to deconvolve a mixture of subclones with copy-number variation from the result of {@link AllelicCNV}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class TumorHeterogeneityModeller {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityModeller.class);

    private static final int NUM_SAMPLES_PER_LOG_ENTRY = 500;
    private static final int MAX_NUM_PROPOSALS_INITIAL_WALKER_BALL = 25;
    private static final double SCALE_PARAMETER = 2.;

    private final EnsembleBuilder<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> builder;
    private final ParameterizedModel<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> model;
    private final TumorHeterogeneityData data;
    private final int numWalkers;

    private final List<Double> concentrationSamples = new ArrayList<>();
    private final List<Double> copyRatioNormalizationSamples = new ArrayList<>();
    private final List<Double> copyRatioNoiseConstantSamples = new ArrayList<>();
    private final List<Double> ploidySamples = new ArrayList<>();
    private final List<PopulationMixture> populationMixtureSamples = new ArrayList<>();

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final int numPopulations,
                                      final int numWalkers,
                                      final double initialWalkerBallSize,
                                      final RandomGenerator rng) {
        this(data, TumorHeterogeneityState.initializeNormalState(data.priors(), data.numSegments(), numPopulations), numWalkers, initialWalkerBallSize, rng);
    }

    public TumorHeterogeneityModeller(final TumorHeterogeneityData data,
                                      final TumorHeterogeneityState initialState,
                                      final int numWalkers,
                                      final double initialWalkerBallSize,
                                      final RandomGenerator rng) {
        Utils.nonNull(data);
        Utils.nonNull(initialState);
        Utils.validateArg(numWalkers >= 2, "Number of walkers must be greater than or equal to two.");
        Utils.validateArg(initialWalkerBallSize > 0., "Initial walker-ball size must be positive.");
        Utils.nonNull(rng);

        this.data = data;
        this.numWalkers = numWalkers;

        //define log-target function
        final Function<TumorHeterogeneityState, Double> logTargetTumorHeterogeneity = state ->
                TumorHeterogeneityUtils.calculateLogJacobianFactor(state, data) + TumorHeterogeneityUtils.calculateLogPosterior(state, data);

        //enumerate copy-number product states
        final int numVariantPopulations = initialState.populationMixture().numVariantPopulations();
        final List<PloidyState> ploidyStates = data.priors().ploidyStatePrior().ploidyStates();
        final Set<Integer> totalCopyNumberStates = ploidyStates.stream().map(PloidyState::total).collect(Collectors.toSet());
        final List<List<Integer>> totalCopyNumberProductStates =
                new ArrayList<>(Sets.cartesianProduct(Collections.nCopies(numVariantPopulations, totalCopyNumberStates)));
        //enumerate ploidy-state product states for each total-copy-number state
        final Map<Integer, Set<PloidyState>> ploidyStateSetsMap = new HashMap<>();
        for (final int totalCopyNumber : totalCopyNumberStates) {
            final Set<PloidyState> ploidyStateSet = ploidyStates.stream().filter(ps -> ps.total() == totalCopyNumber).collect(Collectors.toSet());
            ploidyStateSetsMap.put(totalCopyNumber, ploidyStateSet);
        }
        //define walker transformation
        final Function<WalkerPosition, TumorHeterogeneityState> transformWalkerPositionToState = walkerPosition ->
                TumorHeterogeneityUtils.transformWalkerPositionToState(walkerPosition, data, totalCopyNumberProductStates, ploidyStateSetsMap);

        //initialize walker positions in a ball around initialState
        final List<WalkerPosition> initialWalkerPositions = initializeWalkerBall(rng, initialState, initialWalkerBallSize, logTargetTumorHeterogeneity, transformWalkerPositionToState);

        builder = new EnsembleBuilder<>(SCALE_PARAMETER, initialWalkerPositions, data, transformWalkerPositionToState, logTargetTumorHeterogeneity);
        model = builder.build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link TumorHeterogeneityState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    public void fitMCMC(final int numSamples, final int numBurnIn) {
        Utils.validateArg(numSamples > 0, "Total number of samples must be positive.");
        Utils.validateArg(0 <= numBurnIn && numBurnIn < numSamples,
                "Number of burn-in samples to discard must be non-negative and strictly less than total number of samples.");
        //run MCMC
        final ModelSampler<TumorHeterogeneityParameter, TumorHeterogeneityState, TumorHeterogeneityData> modelSampler
                = new ModelSampler<>(numWalkers * numSamples, model);
        modelSampler.setNumSamplesPerLogEntry(NUM_SAMPLES_PER_LOG_ENTRY);
        modelSampler.runMCMC();
        //update posterior samples
        concentrationSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.CONCENTRATION,
                Double.class, numWalkers * numBurnIn));
        copyRatioNormalizationSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.COPY_RATIO_NORMALIZATION,
                Double.class, numWalkers * numBurnIn));
        copyRatioNoiseConstantSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.COPY_RATIO_NOISE_CONSTANT,
                Double.class, numWalkers * numBurnIn));
        ploidySamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.PLOIDY,
                Double.class, numWalkers * numBurnIn));
        //collapse populations
        populationMixtureSamples.addAll(modelSampler.getSamples(TumorHeterogeneityParameter.POPULATION_MIXTURE,
                PopulationMixture.class, numWalkers * numBurnIn).stream()
                .map(pm -> pm.collapseNormalPopulations(data.priors().normalPloidyState()))
                .collect(Collectors.toList()));
        logger.info("Final acceptance rate: " + builder.calculateAcceptanceRate());
        logger.info("Maximum log posterior: " + builder.getMaxLogTarget());
    }

    public List<Double> getConcentrationSamples() {
        return Collections.unmodifiableList(concentrationSamples);
    }

    public List<Double> getCopyRatioNormalizationSamples() {
        return Collections.unmodifiableList(copyRatioNormalizationSamples);
    }
    
    public List<Double> getCopyRatioNoiseConstantSamples() {
        return Collections.unmodifiableList(copyRatioNoiseConstantSamples);
    }

    public List<Double> getPloidySamples() {
        return Collections.unmodifiableList(ploidySamples);
    }

    public List<PopulationMixture.PopulationFractions> getPopulationFractionsSamples() {
        return Collections.unmodifiableList(populationMixtureSamples.stream()
                .map(PopulationMixture::populationFractions).collect(Collectors.toList()));
    }

    public List<PopulationMixture.VariantProfileCollection> getVariantProfileCollectionSamples() {
        return Collections.unmodifiableList(populationMixtureSamples.stream()
                .map(PopulationMixture::variantProfileCollection).collect(Collectors.toList()));
    }

    public TumorHeterogeneityData getData() {
        return data;
    }

    /**
     * Returns the maximum a posteriori state that was sampled by the {@link ParameterizedModel.EnsembleBuilder} over
     * the entire sampling run.  Note that this state may not be contained in the internally held samples
     * if it was sampled during burn-in.  Normal populations are collapsed.
     */
    public TumorHeterogeneityState getPosteriorMode() {
        final TumorHeterogeneityState posteriorMode = builder.getMaxLogTargetState();
        return new TumorHeterogeneityState(
                posteriorMode.concentration(),
                posteriorMode.copyRatioNormalization(),
                posteriorMode.copyRatioNoiseConstant(),
                posteriorMode.initialPloidy(),
                posteriorMode.ploidy(),
                posteriorMode.populationMixture().collapseNormalPopulations(data.priors().normalPloidyState()));
    }

    /**
     * Given bin sizes in purity and ploidy, returns a list of the indices of the samples falling in the
     * purity-ploidy bin centered on the given state.
     */
    public List<Integer> collectIndicesOfSamplesInBin(final TumorHeterogeneityState binCenter,
                                                      final double purityBinSize,
                                                      final double ploidyBinSize) {
        Utils.validateArg(0. <= purityBinSize && purityBinSize <= 1., "Invalid purity bin size.");
        Utils.validateArg(0. <= ploidyBinSize && ploidyBinSize <= data.priors().ploidyStatePrior().maxCopyNumber(), "Invalid ploidy bin size.");
        final int numSamples = getConcentrationSamples().size();
        if (numSamples == 0) {
            throw new IllegalStateException("Cannot output modeller result before samples have been generated.");
        }

        //determine purity bins centered on given state
        final int numVariantPopulations = binCenter.populationMixture().numVariantPopulations();
        final List<Double> purityBinCenters = binCenter.populationMixture().populationFractions().subList(0, numVariantPopulations);
        final List<Double> purityBinMins = purityBinCenters.stream().map(c -> Math.max(0., c - purityBinSize / 2)).collect(Collectors.toList());
        final List<Double> purityBinMaxs = purityBinCenters.stream().map(c -> Math.min(1., c + purityBinSize / 2)).collect(Collectors.toList());
        IntStream.range(0, numVariantPopulations).forEach(i -> logger.info("Population fraction " + i + " bin: [" + purityBinMins.get(i) + ", " + purityBinMaxs.get(i) + ")"));

        //determine ploidy bin centered on given state
        final double ploidyBinCenter = binCenter.ploidy();
        final double ploidyBinMin = Math.max(TumorHeterogeneityUtils.PLOIDY_MIN, ploidyBinCenter - ploidyBinSize / 2);
        final double ploidyBinMax = Math.min(data.priors().ploidyStatePrior().maxCopyNumber(), ploidyBinCenter + ploidyBinSize / 2);
        logger.info("Ploidy bin: [" + ploidyBinMin + ", " + ploidyBinMax + ")");

        //collect indices of samples falling into purity-ploidy bin
        final List<PopulationMixture.PopulationFractions> populationFractionsSamples = getPopulationFractionsSamples();

        final List<Integer> sampleIndices = IntStream.range(0, numSamples).boxed()
                .filter(i -> IntStream.range(0, numVariantPopulations).allMatch(pi -> purityBinMins.get(pi) <= populationFractionsSamples.get(i).get(pi)
                        && populationFractionsSamples.get(i).get(pi) < purityBinMaxs.get(pi))
                        && ploidyBinMin <= getPloidySamples().get(i)
                        && getPloidySamples().get(i) < ploidyBinMax)
                .collect(Collectors.toList());

        logger.info("Number of samples in bin: " + sampleIndices.size());
        return sampleIndices;
    }

    /**
     * Randomly initialize walker positions in a ball around a given initial {@link TumorHeterogeneityState}.
     */
    private List<WalkerPosition> initializeWalkerBall(final RandomGenerator rng,
                                                      final TumorHeterogeneityState initialState,
                                                      final double initialWalkerBallSize,
                                                      final Function<TumorHeterogeneityState, Double> logTargetTumorHeterogeneity,
                                                      final Function<WalkerPosition, TumorHeterogeneityState> transformWalkerPositionToState) {
        //number of walker dimensions = number of global parameters + (numPopulations - 1) simplex parameters
        final int numVariantPopulations = initialState.populationMixture().numVariantPopulations();
        final int numDimensions = TumorHeterogeneityUtils.NUM_GLOBAL_PARAMETERS + numVariantPopulations;
        final NormalDistribution ballGaussian = new NormalDistribution(rng, 0., initialWalkerBallSize);
        final WalkerPosition walkerPositionOfInitialState = TumorHeterogeneityUtils.transformStateToWalkerPosition(initialState, data);
        final List<WalkerPosition> initialWalkerPositions = new ArrayList<>(numWalkers);
        for (int walkerIndex = 0; walkerIndex < numWalkers; walkerIndex++) {
            boolean acceptedProposedPosition = false;
            WalkerPosition initialWalkerPosition = walkerPositionOfInitialState;
            for (int proposalIndex = 0; proposalIndex < MAX_NUM_PROPOSALS_INITIAL_WALKER_BALL; proposalIndex++) {
                final WalkerPosition proposedWalkerPosition = new WalkerPosition(
                            IntStream.range(0, numDimensions).boxed()
                                    .map(dimensionIndex -> walkerPositionOfInitialState.get(dimensionIndex) + ballGaussian.sample())
                                    .collect(Collectors.toList()));
                final TumorHeterogeneityState proposedState = transformWalkerPositionToState.apply(proposedWalkerPosition);
                proposedState.values().forEach(p -> logger.debug("Proposed " + p.getName().name() + ": " + p.getValue()));
                //only accept the position if its transformed state is within parameter bounds
                if (Double.isFinite(logTargetTumorHeterogeneity.apply(proposedState))) {
                    initialWalkerPosition = proposedWalkerPosition;
                    acceptedProposedPosition = true;
                    break;
                }
            }
            if (!acceptedProposedPosition) {
                logger.debug("Unable to initialize walker within parameter bounds, using position of initial state. " +
                        "Reduce walker-ball size.");
            }
            initialWalkerPositions.add(initialWalkerPosition);  //if no acceptable position was found within the allowed number of iterations, simply take position corresponding to initial state
        }
        return initialWalkerPositions;
    }
}
