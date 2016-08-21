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
public final class ParameterizedModelUnitTest {
    private enum TestParameter implements ParameterEnum {
        PARAMETER_1
    }

    private static final ParameterizedState<TestParameter> SIMPLE_STATE =
            new ParameterizedState<>(Collections.singletonList(new Parameter<>(TestParameter.PARAMETER_1, 1.)));
    private static final TestDataCollection SIMPLE_DATA =
            new TestDataCollection(Collections.singletonList(Collections.singletonList(1.)));
    private static final ParameterSampler<Double, TestParameter, ParameterizedState<TestParameter>, TestDataCollection> SIMPLE_SAMPLER =
            (rng, state, dataCollection) -> 1.;
    private static final ParameterSampler<Integer, TestParameter, ParameterizedState<TestParameter>, TestDataCollection> BAD_SAMPLER_TYPE =
            (rng, state, dataCollection) -> 1;

    private static final class TestDataCollection implements DataCollection {
        private final List<List<Double>> datasets;

        TestDataCollection(final List<List<Double>> datasets) {
            this.datasets = new ArrayList<>(datasets);
        }
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testMissingSamplersException() {
        final ParameterizedModel.GibbsBuilder<TestParameter, ParameterizedState<TestParameter>, TestDataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA);
        builder.build();
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testDuplicateSamplersException() {
        final ParameterizedModel.GibbsBuilder<TestParameter, ParameterizedState<TestParameter>, TestDataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA);
        builder.addParameterSampler(TestParameter.PARAMETER_1, SIMPLE_SAMPLER, Double.class)
                .addParameterSampler(TestParameter.PARAMETER_1, SIMPLE_SAMPLER, Double.class)
                .build();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadSamplerTypeException() {
        final ParameterizedModel.GibbsBuilder<TestParameter, ParameterizedState<TestParameter>, TestDataCollection> builder =
                new ParameterizedModel.GibbsBuilder<>(SIMPLE_STATE, SIMPLE_DATA);
        builder.addParameterSampler(TestParameter.PARAMETER_1, BAD_SAMPLER_TYPE, Integer.class).build();
    }
}