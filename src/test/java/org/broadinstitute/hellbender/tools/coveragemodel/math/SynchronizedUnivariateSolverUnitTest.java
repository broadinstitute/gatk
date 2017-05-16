package org.broadinstitute.hellbender.tools.coveragemodel.math;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.analysis.solvers.AbstractUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SynchronizedUnivariateSolver}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class SynchronizedUnivariateSolverUnitTest extends BaseTest {

    private static final Function<UnivariateSolverSpecifications,
            AbstractUnivariateSolver> SOLVER_FACTORY = desc -> new BrentSolver(desc.getAbsoluteAccuracy(),
            desc.getRelativeAccuracy(), desc.getFunctionValueAccuracy());

    @Test
    public void testFewEquations() throws InterruptedException {
        final Function<Map<Integer, Double>, Map<Integer, Double>> func = arg ->
            arg.entrySet().stream()
                    .map(entry -> {
                        final int index = entry.getKey();
                        final double x = entry.getValue();
                        switch (index) {
                            case 1:
                                return ImmutablePair.of(index, x*x - 3);
                            case 2:
                                return ImmutablePair.of(index, x*x*x - 4);
                            case 3:
                                return ImmutablePair.of(index, x - 5);
                            default:
                                return null;
                        }
                    }).collect(Collectors.toMap(p -> p.left, p -> p.right));


        final SynchronizedUnivariateSolver solver = new SynchronizedUnivariateSolver(func, SOLVER_FACTORY, 3);
        solver.add(1, 0, 4, 3.5, 1e-7, 1e-7, 20);
        solver.add(2, 0, 3, 0.5, 1e-7, 1e-7, 20);
        solver.add(3, 0, 10, 0.6, 1e-7, 1e-7, 20);

        final Map<Integer, SynchronizedUnivariateSolver.UnivariateSolverSummary> sol = solver.solve();
        Assert.assertEquals(sol.get(1).x, 1.732050, 1e-6);
        Assert.assertEquals(sol.get(2).x, 1.587401, 1e-6);
        Assert.assertEquals(sol.get(3).x, 5.000000, 1e-6);
    }

    @Test
    public void testManyEquations() throws InterruptedException {
        testManyEquationsInstance(10);
        testManyEquationsInstance(100);
        testManyEquationsInstance(1000);
    }

    private void testManyEquationsInstance(final int numEquations) throws InterruptedException {
        final Function<Map<Integer, Double>, Map<Integer, Double>> func = arg ->
                arg.entrySet().stream()
                        .map(entry -> {
                            final int index = entry.getKey();
                            final double x = entry.getValue();
                            return ImmutablePair.of(index, FastMath.pow(x, index) - index);
                        }).collect(Collectors.toMap(p -> p.left, p -> p.right));
        final SynchronizedUnivariateSolver solver = new SynchronizedUnivariateSolver(func, SOLVER_FACTORY, numEquations);
        for (int n = 1; n <= numEquations; n++) {
            solver.add(n, 0, 2, 0.5, 1e-7, 1e-7, 100);
        }
        final Map<Integer, SynchronizedUnivariateSolver.UnivariateSolverSummary> sol = solver.solve();
        for (int n = 1; n <= numEquations; n++) {
            Assert.assertEquals(sol.get(n).x, FastMath.pow(n, 1.0/n), 1e-6);
        }
    }
}
