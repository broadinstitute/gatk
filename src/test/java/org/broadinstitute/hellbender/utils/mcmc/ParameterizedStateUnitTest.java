package org.broadinstitute.hellbender.utils.mcmc;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link ParameterizedState}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ParameterizedStateUnitTest {
    private enum TestParameter implements ParameterEnum {
        PARAMETER_1, PARAMETER_2
    }

    private final class Substate extends ArrayList<Double> {
        private static final long serialVersionUID = 123654L;
        public Substate(final List<Double> other) {
            super(new ArrayList<>(other));
        }
    }

    private static final ParameterizedState<TestParameter> testState =
            new ParameterizedState<>(Arrays.asList(new Parameter<>(TestParameter.PARAMETER_1, 1.), new Parameter<>(TestParameter.PARAMETER_2, 1.)));

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyListException() {
        new ParameterizedState<>(Collections.emptyList());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDuplicateParameterNamesException() {
        new ParameterizedState<>(Arrays.asList(new Parameter<>(TestParameter.PARAMETER_1, 1.), new Parameter<>(TestParameter.PARAMETER_1, 1.)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIncompleteParameterListException() {
        new ParameterizedState<>(Collections.singletonList(new Parameter<>(TestParameter.PARAMETER_1, 1.)));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetBadParameterTypeException() {
        testState.get(TestParameter.PARAMETER_1, Integer.class);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUpdateBadParameterTypeException() {
        testState.update(TestParameter.PARAMETER_1, 1);
    }

    @Test
    public void testUpdateParameter() {
        final ParameterizedState<TestParameter> testState =
                new ParameterizedState<>(Arrays.asList(new Parameter<>(TestParameter.PARAMETER_1, 1.), new Parameter<>(TestParameter.PARAMETER_2, 1.)));
        testState.update(TestParameter.PARAMETER_1, 10.);
        Assert.assertEquals(testState.get(TestParameter.PARAMETER_1, Double.class), 10.);
    }

    @Test
    public void testUpdateSubstateParameterAndCopyState() {
        final Substate substate = new Substate(Collections.nCopies(10, 1.));
        final Substate newSubstate = new Substate(Collections.nCopies(10, 2.));
        final ParameterizedState<TestParameter> testStateWithSubstate =
                new ParameterizedState<>(Arrays.asList(new Parameter<>(TestParameter.PARAMETER_1, 1.), new Parameter<>(TestParameter.PARAMETER_2, substate)));
        final ParameterizedState<TestParameter> stateCopyBeforeUpdate = testStateWithSubstate.copy();
        testStateWithSubstate.update(TestParameter.PARAMETER_2, newSubstate);
        final ParameterizedState<TestParameter> stateCopyAfterUpdate = testStateWithSubstate.copy();
        double sumOfCopyBeforeUpdate = stateCopyBeforeUpdate.get(TestParameter.PARAMETER_2, Substate.class).stream().mapToDouble(Double::doubleValue).sum();
        double sumOfCopyAfterUpdate = stateCopyAfterUpdate.get(TestParameter.PARAMETER_2, Substate.class).stream().mapToDouble(Double::doubleValue).sum();
        Assert.assertEquals(sumOfCopyBeforeUpdate, 10.);
        Assert.assertEquals(sumOfCopyAfterUpdate, 20.);
    }
}