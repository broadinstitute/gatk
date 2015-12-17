package org.broadinstitute.hellbender.utils.mcmc;

import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link ParameterizedModel}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ParameterizedModelTest {
    private static final ParameterizedState SIMPLE_STATE =
            new ParameterizedState(Collections.singletonList(new Parameter<>("parameter", 1.)));
    private static final TestDataCollection EMPTY_DATA =
            new TestDataCollection(Collections.emptyList());
    private static final TestDataCollection SIMPLE_DATA =
            new TestDataCollection(Collections.singletonList(Collections.singletonList(1.)));
    private static final Sampler<Double, ParameterizedState, DataCollection> SIMPLE_SAMPLER =
            (rng, state, dataCollection) -> 1.;
    private static final Sampler<Integer, ParameterizedState, DataCollection> BAD_SAMPLER =
            (rng, state, dataCollection) -> 1;

    private static final class TestDataCollection implements DataCollection {
        private final List<List<Double>> datasets;

        public TestDataCollection(final List<List<Double>> datasets) {
            this.datasets = new ArrayList<>(datasets);
        }
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testMissingSamplersException() {
        final ParameterizedModel.GibbsBuilder<ParameterizedState, DataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA, ParameterizedState.class);
        builder.build();
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testDuplicateSamplersException() {
        final ParameterizedModel.GibbsBuilder<ParameterizedState, DataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA, ParameterizedState.class);
        builder.addParameterSampler("parameter", SIMPLE_SAMPLER, Double.class)
                .addParameterSampler("parameter", SIMPLE_SAMPLER, Double.class)
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadSamplerNameException() {
        final ParameterizedModel.GibbsBuilder<ParameterizedState, DataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA, ParameterizedState.class);
        builder.addParameterSampler("notParameter", SIMPLE_SAMPLER, Double.class).build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadSamplerTypeException() {
        final ParameterizedModel.GibbsBuilder<ParameterizedState, DataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA, ParameterizedState.class);
        builder.addParameterSampler("parameter", BAD_SAMPLER, Integer.class).build();
    }
}