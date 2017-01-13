package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.FourierLinearOperatorNDArray;
import org.broadinstitute.hellbender.tools.coveragemodel.linalg.GeneralLinearOperator;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

import javax.annotation.Nonnull;
import java.util.stream.IntStream;

/**
 * TODO github/gatk-protected issue #701 -- this class is part of the CNV-avoiding regularizer and will undergo
 * significant changes soon. To the reviewer: let this go through as is for the time being.
 *
 * @implNote (for 10 latents and ~ 200k targets) the forward + inverse Fourier transform part takes about
 * ~ 300ms together, and the matrix inversion takes about ~ 2.5s. Does it make sense
 * to sparkify the latter?
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelWPreconditionerLocal extends GeneralLinearOperator<INDArray> {

    private final Logger logger = LogManager.getLogger(CoverageModelWLinearOperatorLocal.class);

    private final int numLatents, numTargets, fftSize;
    private final INDArray Q_ll, Z_ll;
    private final FourierLinearOperatorNDArray F_tt;
    private double[] orderedFourierFactors;

    public CoverageModelWPreconditionerLocal(@Nonnull final INDArray Q_ll,
                                             @Nonnull final INDArray Z_ll,
                                             @Nonnull final FourierLinearOperatorNDArray F_tt,
                                             final int numTargets) {
        this.numTargets = ParamUtils.isPositive(numTargets, "Number of target must be positive");
        this.numLatents = Q_ll.shape()[0];
        this.fftSize = F_tt.getFFTSize();
        if (Q_ll.shape()[1] != numLatents)
            throw new IllegalArgumentException("Malformed Q_ll.");
        if (Z_ll.shape()[0] != numLatents || Z_ll.shape()[1] != numLatents)
            throw new IllegalArgumentException("Malformed Z_ll.");
        if (F_tt.getRowDimension() != numTargets || F_tt.getColumnDimension() != numTargets)
            throw new IllegalArgumentException("Malformed F_tt.");
        this.Q_ll = Q_ll;
        this.Z_ll = Z_ll;
        this.F_tt = F_tt;
        this.orderedFourierFactors = F_tt.getOrderedFourierFactors();
    }

    @Override
    public INDArray operate(@Nonnull final INDArray W_tl)
            throws DimensionMismatchException {
        if (W_tl.rank() != 2 || W_tl.shape()[0] != numTargets || W_tl.shape()[1] != numLatents) {
            throw new DimensionMismatchException(W_tl.length(), numTargets * numLatents);
        }
        /* take a Fourier transform in target space */
        long startTimeRFFT = System.nanoTime();
        final INDArray W_kl = Nd4j.create(fftSize, numLatents);
        IntStream.range(0, numLatents).parallel().forEach(li ->
            W_kl.get(NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                    Nd4j.create(F_tt.getForwardFFT(W_tl.get(NDArrayIndex.all(), NDArrayIndex.point(li))),
                            new int[]{fftSize, 1})));
        long endTimeRFFT = System.nanoTime();

        /* apply the preconditioner in the Fourier space */
        long startTimePrecond = System.nanoTime();
        IntStream.range(0, fftSize).parallel().forEach(k -> {
            final INDArray res = CoverageModelEMWorkspaceMathUtils.linsolve(
                    Z_ll.mul(orderedFourierFactors[k]).addi(Q_ll), W_kl.get(NDArrayIndex.point(k)));
            W_kl.get(NDArrayIndex.point(k)).assign(res);
        });
        long endTimePrecond = System.nanoTime();

        /* irfft */
        long startTimeIRFFT = System.nanoTime();
        final INDArray res = Nd4j.create(numTargets, numLatents);
        IntStream.range(0, numLatents).parallel().forEach(li ->
            res.get(NDArrayIndex.all(), NDArrayIndex.point(li)).assign(
                    F_tt.getInverseFFT(W_kl.get(NDArrayIndex.all(), NDArrayIndex.point(li)))));
        long endTimeIRFFT = System.nanoTime();

        logger.debug("Local FFT timing: " + (endTimeRFFT - startTimeRFFT + endTimeIRFFT - startTimeIRFFT)/1000000 + " ms");
        logger.debug("Local preconditioner application timing: " + (endTimePrecond - startTimePrecond)/1000000 + " ms");
        return res;
    }

    @Override
    public int getRowDimension() {
        return numTargets * numLatents;
    }

    @Override
    public int getColumnDimension() {
        return numTargets * numLatents;
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
