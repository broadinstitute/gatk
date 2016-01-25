package org.broadinstitute.hellbender.utils.mcmc;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link AbstractParameterizedState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AbstractParameterizedStateUnitTest {
    private final class TestSubstate extends AbstractParameterizedState {
        private static final String SUBSTATE_PARAMETER_PREFIX = "theta";

        public TestSubstate(final TestSubstate state) {
            super(state);
        }

        @Override
        protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
            return stateClass.cast(new TestSubstate(this));
        }

        public <T> TestSubstate(final List<T> parameterValues) {
            super(SUBSTATE_PARAMETER_PREFIX, parameterValues);
        }

        public double getTheta(final int index) {
            return get(SUBSTATE_PARAMETER_PREFIX + index, Double.class);
        }
    }

    private final class TestState extends AbstractParameterizedState {
        private static final String PARAMETER_NAME = "parameter";
        private static final String SUBSTATE_NAME = "substate";

        public TestState(final TestState state) {
            super(state);
        }

        @Override
        protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
            return stateClass.cast(new TestState(this));
        }

        public TestState(final List<Parameter<?>> parameters) {
            super(parameters);
        }

        public TestState(final double parameterValue, final List<Double> thetaValues) {
            super(Arrays.asList(
                    new Parameter<>(PARAMETER_NAME, parameterValue),
                    new Parameter<>(SUBSTATE_NAME, new TestSubstate(thetaValues))));
        }

        public double getParameter() {
            return get(PARAMETER_NAME, Double.class);
        }

        public double getTheta(final int index) {
            return get(SUBSTATE_NAME, TestSubstate.class).getTheta(index);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicateParameterNamesException() {
        new TestState(Arrays.asList(new Parameter<>("parameter", 1.), new Parameter<>("parameter", 1.)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetBadParameterNameException() {
        final TestState state = new TestState(Collections.singletonList(new Parameter<>("parameter1", 1.)));
        state.get("parameter2", Double.class);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetBadParameterTypeException() {
        final TestState state = new TestState(Collections.singletonList(new Parameter<>("parameter1", 1.)));
        state.get("parameter1", Integer.class);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUpdateBadParameterNameException() {
        final TestState state = new TestState(Collections.singletonList(new Parameter<>("parameter1", 1.)));
        state.updateParameter("parameter2", 1.);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUpdateBadParameterTypeException() {
        final TestState state = new TestState(Collections.singletonList(new Parameter<>("parameter1", 1.)));
        state.updateParameter("parameter1", 1);
    }

    @Test
    public void testUpdateParameter() {
        final TestState state = new TestState(Collections.singletonList(new Parameter<>("parameter", 1.)));
        state.updateParameter("parameter", 10.);
        Assert.assertEquals(state.getParameter(), 10.);
    }

    @Test
    public void testInitializeParameterCopies() {
        final TestSubstate substate = new TestSubstate(Collections.nCopies(10, 1.));
        Assert.assertEquals(IntStream.range(0, 10).boxed().mapToDouble(substate::getTheta).sum(), 10.);
    }

    @Test
    public void testUpdateSubstateParameterAndCopyState() {
        final TestState state = new TestState(5., Collections.nCopies(10, 1.));
        final TestSubstate newSubstate = new TestSubstate(Collections.nCopies(10, 5.));
        state.updateParameter("substate", newSubstate);
        final TestState stateCopy = state.copy(TestState.class);
        Assert.assertEquals(IntStream.range(0, 10).boxed().mapToDouble(stateCopy::getTheta).sum(), 50.);
    }
}