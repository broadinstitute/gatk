package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;

import java.util.HashMap;
import java.util.Map;

/**
 * Represents a parameterized model.  The parameterized state of the model is represented by an
 * {@link AbstractParameterizedState}, while the data is represented by an {@link DataCollection}.
 * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of use.
 * @param <S1>  type of the AbstractParameterizedState
 * @param <T1>  type of the DataCollection
 */
public final class ParameterizedModel<S1 extends AbstractParameterizedState, T1 extends DataCollection> {
    //enums for specifying method of updating parameters (currently, only Gibbs sampling is implemented)
    //updateMethod should be set accordingly within Builders and constructors corresponding to each update method
    private enum UpdateMethod {
        GIBBS
    }

    private final S1 state;
    private final T1 dataCollection;
    private final Class<S1> stateClass;
    private final Map<String, Sampler<?, S1, T1>> samplerMap;
    private final UpdateMethod updateMethod;

    /**
     * Builder for constructing a ParameterizedModel to be Gibbs sampled using {@link GibbsSampler}.
     * Given an initial instance "initialState" of a ConcreteParameterizedState (which extends
     * {@link AbstractParameterizedState}) and an instance "dataset" of a ConcreteDataCollection (which extends
     * {@link DataCollection}), as well as i = 1,...,N {@link Sampler} objects SAMPLER_i that return samples
     * of type TYPE_i, a ParameterizedModel model can be constructed using the Builder pattern as:
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
     * @param <S2>  type of the AbstractParameterizedState
     * @param <T2>  type of the DataCollection
     */
    public static final class GibbsBuilder<S2 extends AbstractParameterizedState, T2 extends DataCollection> {
        private final S2 state;
        private final T2 dataCollection;
        private final Class<S2> stateClass;
        private final Map<String, Sampler<?, S2, T2>> samplerMap = new HashMap<>();

        /**
         * Constructor for {@link ParameterizedModel.GibbsBuilder}.
         * @param state             AbstractParameterizedState held by the model
         * @param dataCollection    DataCollection used by the model
         * @param stateClass        class of AbstractParameterizedState held by the model
         */
        public GibbsBuilder(final S2 state, final T2 dataCollection, final Class<S2> stateClass) {
            this.state = state;
            this.dataCollection = dataCollection;
            this.stateClass = stateClass;
        }

        /**
         * Adds a {@link Sampler} to the collection of parameter samplers using {@link ParameterizedModel.GibbsBuilder}.
         * @param parameterName         name of parameter to sample
         * @param sampler               Sampler that returns random samples of the parameter
         * @param parameterValueClass   class of the parameter value to sample
         * @param <U>                   type of the parameter value to sample
         */
        public <U> GibbsBuilder<S2, T2> addParameterSampler(final String parameterName,
                                                            final Sampler<U, S2, T2> sampler,
                                                            final Class<U> parameterValueClass) {
            if (samplerMap.containsKey(parameterName)) {
                throw new UnsupportedOperationException("Cannot add more than one sampler per parameter.");
            }
            if (!state.parameterNames().contains(parameterName)) {
                throw new IllegalArgumentException("Cannot add sampler for parameter not specified in initial state.");
            }
            try {
                state.get(parameterName, parameterValueClass);
            } catch (final IllegalArgumentException e) {
                throw new IllegalArgumentException("Cannot add sampler for parameter that returns type different " +
                        "than that specified for parameter in initial state.");
            }
            samplerMap.put(parameterName, sampler);
            return this;
        }

        /**
         * Builds the ParameterizedModel as specified via {@link ParameterizedModel.GibbsBuilder}.
         * @return  ParameterizedModel as specified via GibbsBuilder
         * @throws UnsupportedOperationException if there is not a one-to-one mapping between Parameters in the
         *                                       AbstractParameterizedState and the Samplers specified via GibbsBuilder
         */
        public ParameterizedModel<S2, T2> build() {
            if (!(samplerMap.keySet().containsAll(state.parameterNames()) &&
                    state.parameterNames().containsAll(samplerMap.keySet()))) {
                throw new UnsupportedOperationException("Each parameter must have a corresponding sampler specified.");
            }
            return new ParameterizedModel<>(this);
        }
    }

    //Constructor for GibbsBuilder
    private ParameterizedModel(final GibbsBuilder<S1, T1> builder) {
        this.state = builder.state;
        this.dataCollection = builder.dataCollection;
        this.stateClass = builder.stateClass;
        this.samplerMap = builder.samplerMap;
        this.updateMethod = UpdateMethod.GIBBS;
    }

    /**
     * Returns a copy of the {@link AbstractParameterizedState} held internally, as given by the overridden
     * {@link AbstractParameterizedState#copy(Class)} method.
     * @return  copy of the AbstractParameterizedState held internally
     */
    protected S1 state() {
        return state.copy(stateClass);
    }

    /**
     * Updates the {@link AbstractParameterizedState} held internally using the Samplers and update method specified via
     * the Builder pattern.  Currently, only Gibbs sampling is implemented.
     * @param rng   RandomGenerator to pass to Samplers to generate samples
     */
    protected void update(final RandomGenerator rng) {
        if (updateMethod == UpdateMethod.GIBBS) {
            doGibbsUpdate(rng);
        }
    }

    private void doGibbsUpdate(final RandomGenerator rng) {
        for (final String parameterName : state.parameterNames()) {
            state.updateParameter(parameterName, samplerMap.get(parameterName).sample(rng, state, dataCollection));
        }
    }
}
