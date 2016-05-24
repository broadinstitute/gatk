package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.random.RandomGenerator;

import java.util.HashMap;
import java.util.Map;

/**
 * Represents a parameterized model.  The parameterized state of the model is represented by an
 * {@link ParameterizedState}, while the data is represented by an {@link DataCollection}.
 * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of use.
 * @param <S1>  type of the ParameterizedState
 * @param <T1>  type of the DataCollection
 */
public final class ParameterizedModel<V1 extends Enum<V1> & ParameterEnum, S1 extends ParameterizedState<V1>, T1 extends DataCollection> {
    //enums for specifying method of updating parameters (currently, only Gibbs sampling is implemented)
    //updateMethod should be set accordingly within Builders and constructors corresponding to each update method
    protected enum UpdateMethod {
        GIBBS
    }

    private final S1 state;
    private final T1 dataCollection;
    private final Map<V1, ParameterSampler<?, V1, S1, T1>> samplerMap;
    private final UpdateMethod updateMethod;

    /**
     * Builder for constructing a ParameterizedModel to be Gibbs sampled using {@link GibbsSampler}.
     * Given an initial instance "initialState" of a ConcreteParameterizedState (which extends
     * {@link ParameterizedState}) and an instance "dataset" of a ConcreteDataCollection (which extends
     * {@link DataCollection}), as well as i = 1,...,N {@link ParameterSampler} objects SAMPLER_i that return samples
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
            this.state = state;
            this.dataCollection = dataCollection;
        }

        /**
         * Adds a {@link ParameterSampler} to the collection of parameter samplers using {@link ParameterizedModel.GibbsBuilder}.
         * @param parameterName         name of parameter to sample
         * @param parameterSampler               ParameterSampler that returns random samples of the parameter
         * @param parameterValueClass   class of the parameter value to sample
         * @param <U>                   type of the parameter value to sample
         */
        public <U> GibbsBuilder<V2, S2, T2> addParameterSampler(final V2 parameterName,
                                                                final ParameterSampler<U, V2, S2, T2> parameterSampler,
                                                                final Class<U> parameterValueClass) {
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
         * Builds the ParameterizedModel as specified via {@link ParameterizedModel.GibbsBuilder}.
         * @return  ParameterizedModel as specified via GibbsBuilder
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

    //Constructor for GibbsBuilder
    private ParameterizedModel(final GibbsBuilder<V1, S1, T1> builder) {
        this.state = builder.state;
        this.dataCollection = builder.dataCollection;
        this.samplerMap = builder.samplerMap;
        this.updateMethod = UpdateMethod.GIBBS;
    }

    /**
     * Returns a copy of the {@link ParameterizedState} held internally.
     * @return  copy of the {@link ParameterizedState} held internally
     */
    protected S1 state() {
        return state.copy();
    }

    /**
     * Updates the {@link ParameterizedState} held internally using the {@link ParameterSampler}s and update method specified via
     * the Builder pattern.  Currently, only Gibbs sampling is implemented.
     * @param rng   {@link RandomGenerator} to pass to {@link ParameterSampler}s to generate samples
     */
    protected void update(final RandomGenerator rng) {
        if (updateMethod == UpdateMethod.GIBBS) {
            doGibbsUpdate(rng);
        }
    }

    protected UpdateMethod getUpdateMethod() {
        return updateMethod;
    }

    private void doGibbsUpdate(final RandomGenerator rng) {
        for (final V1 parameterName : state.keySet()) {
            state.update(parameterName, samplerMap.get(parameterName).sample(rng, state, dataCollection));
        }
    }
}
