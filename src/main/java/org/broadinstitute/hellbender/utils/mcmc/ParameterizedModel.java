package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a parameterized model.  The parameterized state of the model is represented by an
 * {@link ParameterizedState}, while the data is represented by an {@link DataCollection}.
 * Allows for sampling via Gibbs sampling (if {@link ParameterizedModel.GibbsBuilder} is used in construction)
 * or affine-invariant ensemble sampling (if {@link ParameterizedModel.EnsembleBuilder} is used in construction).
 * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of Gibbs sampling.
 * @param <S1>  type of the ParameterizedState
 * @param <T1>  type of the DataCollection
 */
public final class ParameterizedModel<V1 extends Enum<V1> & ParameterEnum, S1 extends ParameterizedState<V1>, T1 extends DataCollection> {
    private static final Logger logger = LogManager.getLogger(ParameterizedModel.class);

    private final S1 state;
    private final T1 dataCollection;
    private final Consumer<RandomGenerator> updateState;

    /**
     * Builder for constructing a {@link ParameterizedModel} to be Gibbs sampled using {@link ModelSampler}.
     * Given an initial instance "initialState" of a ConcreteParameterizedState (which extends
     * {@link ParameterizedState}) and an instance "dataset" of a ConcreteDataCollection (which extends
     * {@link DataCollection}), as well as i = 1,...,N {@link ParameterSampler} objects SAMPLER_i that return samples
     * of type TYPE_i, a {@link ParameterizedModel} model can be constructed using the Builder pattern as:
     *
     *  ParameterizedModel<ConcreteParameterizedState, ConcreteDataCollection> model =
     *      new ParameterizedModel.GibbsBuilder<>(initialState, dataset, ConcreteParameterizedState.class)
     *                            .addParameterSampler(SAMPLER_1, TYPE_1.class)
     *                            .
     *                            .
     *                            .
     *                            .addParameterSampler(SAMPLER_N, TYPE_N.class)
     *                            .build()
     *
     * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of use.
     * @param <V2>  type of the ParameterEnum
     * @param <S2>  type of the ParameterizedState
     * @param <T2>  type of the DataCollection
     */
    public static final class GibbsBuilder<V2 extends Enum<V2> & ParameterEnum, S2 extends ParameterizedState<V2>, T2 extends DataCollection> {
        private final S2 state;
        private final T2 dataCollection;
        private final Map<V2, ParameterSampler<?, V2, S2, T2>> samplerMap = new HashMap<>();

        /**
         * Constructor for {@link ParameterizedModel.GibbsBuilder}.
         * @param state             ParameterizedState held by the model
         * @param dataCollection    DataCollection used by the model
         */
        public GibbsBuilder(final S2 state, final T2 dataCollection) {
            Utils.nonNull(state);
            Utils.nonNull(dataCollection);
            this.state = state;
            this.dataCollection = dataCollection;
        }

        /**
         * Adds a {@link ParameterSampler} to the collection of parameter samplers using {@link ParameterizedModel.GibbsBuilder}.
         * @param parameterName         name of parameter to sample
         * @param parameterSampler      ParameterSampler that returns random samples of the parameter
         * @param parameterValueClass   class of the parameter value to sample
         * @param <U>                   type of the parameter value to sample
         */
        public <U> GibbsBuilder<V2, S2, T2> addParameterSampler(final V2 parameterName,
                                                                final ParameterSampler<U, V2, S2, T2> parameterSampler,
                                                                final Class<U> parameterValueClass) {
            Utils.nonNull(parameterName);
            Utils.nonNull(parameterSampler);
            Utils.nonNull(parameterValueClass);
            if (samplerMap.containsKey(parameterName)) {
                throw new UnsupportedOperationException("Cannot add more than one sampler per parameter.");
            }
            try {
                state.get(parameterName, parameterValueClass);
            } catch (final IllegalArgumentException e) {
                throw new IllegalArgumentException("Cannot add sampler for parameter that returns type different " +
                        "than that specified for parameter in initial state.");
            }
            samplerMap.put(parameterName, parameterSampler);
            return this;
        }

        /**
         * Builds the {@link ParameterizedModel} as specified via {@link ParameterizedModel.GibbsBuilder}.
         * @return {@link ParameterizedModel} as specified via GibbsBuilder
         * @throws UnsupportedOperationException if there is not a one-to-one mapping between Parameters in the
         *                                       {@link ParameterizedState} and the {@link ParameterSampler}s
         *                                       specified via GibbsBuilder
         */
        public ParameterizedModel<V2, S2, T2> build() {
            if (!samplerMap.keySet().equals(state.keySet())) {
                throw new UnsupportedOperationException("Each parameter must have a corresponding sampler specified.");
            }
            return new ParameterizedModel<>(this);
        }
    }

    /**
     * Constructor using {@link ParameterizedModel.GibbsBuilder}.
     */
    private ParameterizedModel(final GibbsBuilder<V1, S1, T1> builder) {
        state = builder.state;
        dataCollection = builder.dataCollection;
        updateState = rng -> doGibbsUpdate(builder, rng);
    }

    /**
     * Builder for constructing a {@link ParameterizedModel} to be ensemble sampled using {@link ModelSampler}.
     * Given a list of initial walker positions for all walkers in the ensemble, 
     * an instance "dataset" of a ConcreteDataCollection (which extends {@link DataCollection}), 
     * a function for transforming walker positions to ConcreteParameterizedStates, 
     * and a log-target function to be sampled, this builder holds the necessary state for
     * performing affine-invariant ensemble sampling according to Goodman & Weare 2010.
     * Walkers in the ensemble are iterated over in turn, with the ConcreteParameterizedState held internally 
     * by the {@link ParameterizedModel} corresponding to the currently selected walker; thus, a single iteration
     * over the entire ensemble corresponds to a number of sampled ConcreteParameterizedStates
     * held by the {@link ModelSampler} equal to the number of walkers.  The ConcreteParameterizedState
     * at which the maximum value of the log-target function is observed is also stored.
     * @param <V2>  type of the ParameterEnum
     * @param <S2>  type of the ParameterizedState
     * @param <T2>  type of the DataCollection
     */
    public static final class EnsembleBuilder<V2 extends Enum<V2> & ParameterEnum, S2 extends ParameterizedState<V2>, T2 extends DataCollection> {
        private final double scaleParameter;
        private final S2 state;
        private final T2 dataCollection;

        private int selectedWalkerIndex;
        private int numAccepted;
        private int numSamples;
        private final int numWalkers;
        private final Function<WalkerPosition, S2> transformWalkerPositionToState;
        private final Function<S2, Double> logTarget;
        private final List<WalkerTuple<V2, S2>> currentWalkerTuples;
        private final List<WalkerTuple<V2, S2>> nextWalkerTuples;
        private WalkerTuple<V2, S2> maxWalkerTuple;

        /**
         * Constructor for {@link ParameterizedModel.EnsembleBuilder}.
         * @param scaleParameter                    scale parameter for stretch-move proposal;
         *                                          see Goodman & Weare 2010, which recommends setting this to 2
         * @param initialWalkerPositions            initial positions of the walkers in N-dimensional space
         * @param dataCollection                    DataCollection used by the model
         * @param transformWalkerPositionToState    function that transforms a walker position to a state
         * @param logTarget                         log of the target function to sample
         */
        public EnsembleBuilder(final double scaleParameter,
                               final List<WalkerPosition> initialWalkerPositions,
                               final T2 dataCollection,
                               final Function<WalkerPosition, S2> transformWalkerPositionToState,
                               final Function<S2, Double> logTarget) {
            Utils.validateArg(scaleParameter > 1.,
                    "Scale parameter should be strictly greater than 1.");
            Utils.nonNull(initialWalkerPositions);
            Utils.nonNull(dataCollection);
            Utils.nonNull(transformWalkerPositionToState);
            Utils.nonNull(logTarget);
            Utils.validateArg(initialWalkerPositions.size() >= 2,
                    "Number of walkers must be greater than or equal to two.");
            Utils.validateArg(initialWalkerPositions.stream()
                    .map(WalkerPosition::numDimensions).allMatch(n -> n.equals(initialWalkerPositions.get(0).numDimensions())),
                    "Dimension of walker space must be identical for all walkers.");

            this.scaleParameter = scaleParameter;
            state = transformWalkerPositionToState.apply(initialWalkerPositions.get(0));
            this.dataCollection = dataCollection;
            selectedWalkerIndex = 0;
            numAccepted = 0;
            numSamples = 0;
            numWalkers = initialWalkerPositions.size();
            this.transformWalkerPositionToState = transformWalkerPositionToState;
            this.logTarget = logTarget;

            //initialize both current and next WalkerTuples using initial WalkerPositions;
            //next WalkerTuples will be updated on first iteration
            currentWalkerTuples = initialWalkerPositions.stream().map(wp -> new WalkerTuple<>(wp, transformWalkerPositionToState, logTarget)).collect(Collectors.toList());
            nextWalkerTuples = initialWalkerPositions.stream().map(wp -> new WalkerTuple<>(wp, transformWalkerPositionToState, logTarget)).collect(Collectors.toList());

            //set maxWalkerTuple (walker with maximum log target value) to first walker
            maxWalkerTuple = currentWalkerTuples.get(0);
        }

        private static final class WalkerTuple<V extends Enum<V> & ParameterEnum, S extends ParameterizedState<V>> {
            private final WalkerPosition walkerPosition;
            private final S state;
            private final double logTargetValue;

            WalkerTuple(final WalkerPosition walkerPosition,
                        final Function<WalkerPosition, S> transformWalkerPositionToState,
                        final Function<S, Double> logTarget) {
                this.walkerPosition = walkerPosition;
                state = transformWalkerPositionToState.apply(walkerPosition);
                logTargetValue = logTarget.apply(state);
            }
        }

        /**
         * Returns the {@link ParameterizedState} at which the maximum value of the log-target function was observed
         * over the entire sampling run.
         */
        public S2 getMaxLogTargetState() {
            return maxWalkerTuple.state;
        }

        /**
         * Returns the maximum value of the log-target function that was observed over the entire sampling run.
         */
        public double getMaxLogTarget() {
            return maxWalkerTuple.logTargetValue;
        }

        /**
         * Calculates and returns the current acceptance rate.
         */
        public double calculateAcceptanceRate() {
            return (double) numAccepted / numSamples;
        }

        /**
         * Builds the {@link ParameterizedModel} as specified via {@link ParameterizedModel.EnsembleBuilder}.
         * @return {@link ParameterizedModel} as specified via EnsembleBuilder
         */
        public ParameterizedModel<V2, S2, T2> build() {
            return new ParameterizedModel<>(this);
        }
    }

    /**
     * Constructor using {@link ParameterizedModel.EnsembleBuilder}.
     */
    private ParameterizedModel(final EnsembleBuilder<V1, S1, T1> builder) {
        state = builder.state;
        dataCollection = builder.dataCollection;
        updateState = rng -> doEnsembleUpdate(builder, rng);
    }

    /**
     * Returns a copy of the {@link ParameterizedState} held internally.
     * @return  copy of the {@link ParameterizedState} held internally
     */
    protected S1 state() {
        return state.copy();
    }

    /**
     * Updates the {@link ParameterizedState} held internally using the {@link ParameterSampler}s
     * and update method specified via the Builder pattern.
     * @param rng   {@link RandomGenerator} to pass to {@link ParameterSampler}s to generate samples
     */
    protected void update(final RandomGenerator rng) {
        updateState.accept(rng);
    }

    private void doGibbsUpdate(final GibbsBuilder<V1, S1, T1> builder, final RandomGenerator rng) {
        for (final V1 parameterName : state.keySet()) {
            state.update(parameterName, builder.samplerMap.get(parameterName).sample(rng, state, dataCollection));
        }
    }

    private void doEnsembleUpdate(final EnsembleBuilder<V1, S1, T1> builder, final RandomGenerator rng) {
        //pick a walker other than the one selected for an update
        int otherWalkerIndex;
        do {
            otherWalkerIndex = rng.nextInt(builder.numWalkers);
        } while (otherWalkerIndex == builder.selectedWalkerIndex);

        //get relevant current WalkerTuples
        final EnsembleBuilder.WalkerTuple<V1, S1> currentWalkerTupleSelected = builder.currentWalkerTuples.get(builder.selectedWalkerIndex);
        final EnsembleBuilder.WalkerTuple<V1, S1> currentWalkerTupleOther = builder.currentWalkerTuples.get(otherWalkerIndex);

        //propose a stretch move
        final int numDimensions = currentWalkerTupleSelected.walkerPosition.numDimensions();
        final double z = FastMath.pow((builder.scaleParameter - 1.) * rng.nextDouble() + 1, 2.) / builder.scaleParameter; //see Goodman & Weare 2010
        final WalkerPosition proposedWalkerPosition = new WalkerPosition(IntStream.range(0, numDimensions).boxed()
                .map(i -> currentWalkerTupleOther.walkerPosition.get(i) + z * (currentWalkerTupleSelected.walkerPosition.get(i) - currentWalkerTupleOther.walkerPosition.get(i)))
                .collect(Collectors.toList()));

        //create proposed WalkerTuple (calculates state and log target
        final EnsembleBuilder.WalkerTuple<V1, S1> proposedWalkerTuple = new EnsembleBuilder.WalkerTuple<>(
                proposedWalkerPosition, builder.transformWalkerPositionToState, builder.logTarget);

        //get states and output parameter values
        final S1 currentState = currentWalkerTupleSelected.state;
        currentState.values().forEach(p -> logger.debug("Current " + p.getName().name() + ": " + p.getValue()));
        final S1 proposedState = proposedWalkerTuple.state;
        proposedState.values().forEach(p -> logger.debug("Proposed " + p.getName().name() + ": " + p.getValue()));

        //get log targets and output
        final double currentLogTarget = currentWalkerTupleSelected.logTargetValue;
        logger.debug("Log target of current state: " + currentLogTarget);
        final double proposedLogTarget = proposedWalkerTuple.logTargetValue;
        logger.debug("Log target of proposed state: " + proposedLogTarget);

        //accept or reject
        final double acceptanceLogProbability = Math.min(0., (numDimensions - 1.) * FastMath.log(z) + proposedLogTarget - currentLogTarget); //see Goodman & Weare 2010
        builder.numSamples++;
        if (FastMath.log(rng.nextDouble()) < acceptanceLogProbability) {
            builder.numAccepted++;
            logger.debug("Proposed state accepted.");
            //update the next WalkerTuple for the selected walker
            builder.nextWalkerTuples.set(builder.selectedWalkerIndex, proposedWalkerTuple);
            //update the state held by the model using the accepted proposed state for the selected walker
            proposedState.values().forEach(p -> state.update(p.getName(), p.getValue()));
            //update the maximum target state if appropriate
            if (proposedLogTarget > builder.maxWalkerTuple.logTargetValue) {
                logger.debug("New maximum found.");
                builder.maxWalkerTuple = proposedWalkerTuple;
            }
        } else {
            //update the state held by the model using the previous state for the selected walker
            currentState.values().forEach(p -> state.update(p.getName(), p.getValue()));
        }
        logger.debug("Acceptance rate: " + (double) builder.numAccepted / builder.numSamples);
        state.values().forEach(p -> logger.debug("Sampled " + p.getName().name() + ": " + p.getValue()));

        //move to the next walker
        builder.selectedWalkerIndex = (builder.selectedWalkerIndex + 1) % builder.numWalkers;

        //if we have iterated over the entire ensemble, update all of the current WalkerTuples
        if (builder.selectedWalkerIndex == 0) {
            IntStream.range(0, builder.numWalkers).forEach(i -> builder.currentWalkerTuples.set(i, builder.nextWalkerTuples.get(i)));
        }
    }
}
