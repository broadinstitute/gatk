package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.broadinstitute.hellbender.tools.coveragemodel.SubroutineSignal;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.testng.annotations.Test;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IterativeLinearSolverNDArrayUnitTest extends BaseTest {

    private static final int DEFAULT_MAX_ITERS = 10;
    private static final double DEFAULT_REL_TOL = 1e-6;
    private static final double DEFAULT_ABS_TOL = 1e-6;

    @Test
    public void testCGSmallSystem() {
        final INDArray A = Nd4j.create(new double[][]{{4, 1}, {1, 3}});
        final INDArray M = Nd4j.create(new double[][]{{0.25, 0}, {0, 1.0}});
        final INDArray b = Nd4j.create(new double[]{1, 2}).reshape(2, 1);
        final INDArray x0 = Nd4j.create(new double[]{2, 1}).reshape(2, 1);
        final INDArray xExpected = Nd4j.create(new double[]{1.0/11, 7.0/11}).reshape(2, 1);

        final GeneralLinearOperatorNDArray linop = new GeneralLinearOperatorNDArray(A);
        final GeneralLinearOperatorNDArray precond = new GeneralLinearOperatorNDArray(M);
        IterativeLinearSolverNDArray solver;
        SubroutineSignal sig;

        /* without preconditioning */
        solver = new IterativeLinearSolverNDArray(linop, b, null, DEFAULT_ABS_TOL, DEFAULT_REL_TOL, DEFAULT_MAX_ITERS,
                x -> x.normmaxNumber().doubleValue(), (x, y) -> x.mul(y).sumNumber().doubleValue(), true);
        sig = solver.cg(x0);
        Assert.assertTrue(sig.getObject("status") == IterativeLinearSolverNDArray.ExitStatus.SUCCESS_ABS_TOL ||
                sig.getObject("status") == IterativeLinearSolverNDArray.ExitStatus.SUCCESS_REL_TOL);
        Assert.assertArrayEquals(sig.getINDArray("x").data().asDouble(), xExpected.data().asDouble(), DEFAULT_ABS_TOL);

        /* with preconditioning */
        solver = new IterativeLinearSolverNDArray(linop, b, precond, DEFAULT_ABS_TOL, DEFAULT_REL_TOL, DEFAULT_MAX_ITERS,
                x -> x.normmaxNumber().doubleValue(), (x, y) -> x.mul(y).sumNumber().doubleValue(), true);
        sig = solver.cg(x0);
        Assert.assertTrue(sig.getObject("status") == IterativeLinearSolverNDArray.ExitStatus.SUCCESS_ABS_TOL ||
                sig.getObject("status") == IterativeLinearSolverNDArray.ExitStatus.SUCCESS_REL_TOL);
        Assert.assertArrayEquals(sig.getINDArray("x").data().asDouble(), xExpected.data().asDouble(), DEFAULT_ABS_TOL);
    }
}
