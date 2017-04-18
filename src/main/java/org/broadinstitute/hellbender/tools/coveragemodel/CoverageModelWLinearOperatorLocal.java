package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperatorNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.GeneralLinearOperator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

import javax.annotation.Nonnull;
import java.util.stream.IntStream;

/**
 * Implements the (real and symmetric) linear operator,
 *
 *      Q_{t}^{\mu\nu}\delta_{t,t'} + F_{t,t'}Z^{\mu\nu}.
 *
 * It acts on W_{t'}^{\nu} and returns an INDArray of similar shape.
 *
 * Note: this class is a local (driver node) implementation.
 *
 * TODO github/gatk-protected issue #701 -- this class is part of the CNV-avoiding regularizer and will undergo
 * significant changes soon. To the reviewer: let this go through as is for the time being.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelWLinearOperatorLocal implements GeneralLinearOperator<INDArray> {

    private static final Logger logger = LogManager.getLogger(CoverageModelWLinearOperatorLocal.class);

    private final int numLatents, numTargets;
    private final INDArray Q_tll, Z_ll;
    private final FourierLinearOperatorNDArray F_tt;

    public CoverageModelWLinearOperatorLocal(@Nonnull final INDArray Q_tll,
                                             @Nonnull final INDArray Z_ll,
                                             @Nonnull final FourierLinearOperatorNDArray F_tt) {
        numTargets = Q_tll.shape()[0];
        numLatents = Q_tll.shape()[1];
        if (Q_tll.rank() != 3 || Q_tll.shape()[2] != numLatents)
            throw new IllegalArgumentException("Malformed Q_tll.");
        if (Z_ll.shape()[0] != numLatents || Z_ll.shape()[1] != numLatents)
            throw new IllegalArgumentException("Malformed Z_ll.");
        if (F_tt.getRowDimension() != numTargets || F_tt.getColumnDimension() != numTargets)
            throw new IllegalArgumentException("Malformed F_tt.");
        this.Q_tll = Q_tll;
        this.Z_ll = Z_ll;
        this.F_tt = F_tt;
    }

    @Override
    public int getRowDimension() {
        return numLatents * numTargets;
    }

    @Override
    public int getColumnDimension() {
        return numLatents * numTargets;
    }

    @Override
    public INDArray operate(@Nonnull final INDArray W_tl) throws DimensionMismatchException {
        if (W_tl.rank() != 2 || W_tl.shape()[0] != numTargets || W_tl.shape()[1] != numLatents)
            throw new DimensionMismatchException(W_tl.length(), numTargets * numLatents);
        final INDArray fx = Nd4j.create(numTargets, numLatents);
        final INDArray res = Nd4j.create(numTargets, numLatents);

        /* F W */
        final long startTimeFW = System.nanoTime();
        IntStream.range(0, numLatents).parallel().forEach(li ->
            fx.get(NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                    F_tt.operate(W_tl.get(NDArrayIndex.all(), NDArrayIndex.point(li)))));
        final long endTimeFW = System.nanoTime();

        /* Q W + Z F W */
        final long startTimeCalc = System.nanoTime();
        IntStream.range(0, numTargets).parallel().forEach(ti ->
            res.get(NDArrayIndex.point(ti)).assign(
                    Q_tll.get(NDArrayIndex.point(ti))
                            .mmul(W_tl.get(NDArrayIndex.point(ti)).transpose())
                            .addi(Z_ll.mmul(fx.get(NDArrayIndex.point(ti)).transpose())).transpose()));
        final long endTimeCalc = System.nanoTime();

        logger.debug("Local [F] [W] timing: " + (endTimeFW - startTimeFW)/1000000 + " ms");
        logger.debug("Local [Q] [W] + [Z] [F] [W] timing: " + (endTimeCalc - startTimeCalc)/1000000 + " ms");

        return res;
    }

    @Override
    public INDArray operateTranspose(@Nonnull final INDArray x)
            throws DimensionMismatchException, UnsupportedOperationException {
        return operate(x);
    }

    @Override
    public boolean isTransposable() {
        return true;
    }
}
