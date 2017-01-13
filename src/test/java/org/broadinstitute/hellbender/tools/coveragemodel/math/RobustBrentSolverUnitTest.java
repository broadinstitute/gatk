package org.broadinstitute.hellbender.tools.coveragemodel.math;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

/**
 * Unit tests for {@link RobustBrentSolver}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class RobustBrentSolverUnitTest extends BaseTest {

    private static final double DEF_ABS_ACC = 1e-6;
    private static final double DEF_REL_ACC = 1e-6;
    private static final double DEF_F_ACC = 1e-15;

    @Test
    public void gridTest() {
        final RobustBrentSolver solver = new RobustBrentSolver(DEF_REL_ACC, DEF_REL_ACC, DEF_F_ACC);
        final double[] x = solver.createFractalSearchGrid(0, 1, 10, 2);
        Assert.assertEquals(x.length, 23);
        Assert.assertEquals(x[0], 0, 1e-12);
        Assert.assertEquals(x[22], 1, 1e-12);
    }

    @Test
    public void detectBracketsTest() {
        final List<RobustBrentSolver.Bracket> brackets = RobustBrentSolver.detectBrackets(
                new double[] {0, 1,  2,  3, 4, 5,  6, 7},
                new double[] {0, 1, -1, -5, 6, 7, -1, 0});
        Assert.assertEquals(brackets.size(), 5);
    }

    @Test
    public void simpleTest() {
        final UnivariateFunction objFunc = x -> 30 * x * (x - 1) * (x - 2) * (x - 3);
        final UnivariateFunction meritFunc = x -> 6 * FastMath.pow(x, 5) - 45 * FastMath.pow(x, 4) + 110 * FastMath.pow(x, 3) -
                90 * FastMath.pow(x, 2);
        final RobustBrentSolver solverRobust = new RobustBrentSolver(DEF_REL_ACC, DEF_REL_ACC, DEF_F_ACC);
        final BrentSolver solverSimple = new BrentSolver(DEF_REL_ACC, DEF_REL_ACC, DEF_F_ACC);
        final double xRobust = solverRobust.solve(100, objFunc, meritFunc, null, -1, 4, 4, 1);
        Assert.assertEquals(xRobust, 0, DEF_ABS_ACC);
        boolean simpleSolverFails = false;
        try {
            /* this will fail */
            solverSimple.solve(100, objFunc, -1, 4);
        } catch (final NoBracketingException ex) {
            simpleSolverFails = true;
        }
        Assert.assertTrue(simpleSolverFails);
    }

}
