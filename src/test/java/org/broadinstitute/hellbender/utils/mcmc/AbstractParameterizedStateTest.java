package org.broadinstitute.hellbender.utils.mcmc;

import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link AbstractParameterizedState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AbstractParameterizedStateTest {
    private final class TestState extends AbstractParameterizedState {
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

        public double get(final String parameterName) {
            return get(parameterName, Double.class);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicateParameterNamesException() {
        new TestState(Arrays.asList(new Parameter<>("parameter", 1.), new Parameter<>("parameter", 1.)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetBadParameterNameException() {
        final TestState state = new TestState(Collections.singletonList(new Parameter<>("parameter1", 1.)));
        state.get("parameter2");
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
}