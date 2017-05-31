package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;

/**
 * A a subclass of {@link FourierLinearOperator} that implements the linear algebra operations
 * using Nd4j.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public final class FourierLinearOperatorNDArray extends FourierLinearOperator<INDArray> {

    public FourierLinearOperatorNDArray(final int dimension, final double[] fourierFactors,
                                        final boolean promoteToPowerOfTwo) {
        super(dimension, fourierFactors, promoteToPowerOfTwo);
    }

    @Override
    public INDArray operate(@Nonnull final INDArray x) throws DimensionMismatchException {
        if (x.length() != getDimension()) {
            throw new DimensionMismatchException(x.length(), getDimension());
        }
        final double[] xDoubleArray = createDoubleArrayFromINDArray(x);
        performForwardFFTInPlace(xDoubleArray);
        applyFilterInPlace(xDoubleArray);
        performInverseFFTInPlace(xDoubleArray);
        return Nd4j.create(xDoubleArray, new int[]{getDimension(), 1});
    }

    /**
     * @param x the input vector
     * @return rfft of x
     */
    public double[] getForwardFFT(@Nonnull final INDArray x) {
        if (x.length() != getDimension()) {
            throw new DimensionMismatchException(x.length(), getDimension());
        }
        final double[] xDoubleArray = createDoubleArrayFromINDArray(x);
        performForwardFFTInPlace(xDoubleArray);
        return xDoubleArray;
    }

    /**
     * @param xFFT the input vector
     * @return rfft of x
     */
    public INDArray getInverseFFT(@Nonnull final double[] xFFT) {
        if (xFFT.length != getFFTSize()) {
            throw new DimensionMismatchException(xFFT.length, getFFTSize());
        }
        final double[] xDoubleArray = xFFT.clone();
        performInverseFFTInPlace(xDoubleArray);
        return Nd4j.create(xDoubleArray, new int[]{getDimension(), 1});
    }

    /**
     * @param xFFT the input vector
     * @return rfft of x
     */
    public INDArray getInverseFFT(@Nonnull final INDArray xFFT) {
        final int fftSize = getFFTSize();
        if (xFFT.length() != fftSize) {
            throw new DimensionMismatchException(xFFT.length(), fftSize);
        }
        final double[] xFFTDoubleArray = new double[fftSize];
        for (int i = 0; i < fftSize; i++) {
            xFFTDoubleArray[i] = xFFT.getDouble(i);
        }
        performInverseFFTInPlace(xFFTDoubleArray);
        return Nd4j.create(xFFTDoubleArray, new int[]{getDimension(), 1});
    }

    private double[] createDoubleArrayFromINDArray(@Nonnull final INDArray x) {
        return zeroPad(new IndexRange(0, getDimension()).mapToDouble(x::getDouble));
    }

    private double[] zeroPad(final double[] x) {
        Utils.validateArg(x.length <= getFFTSize(), () -> String.format("Can not zero-pad an array that is already" +
                " longer than the FFT size; array length = %d, FFT size = %d", x.length, getFFTSize()));
        final double[] zeroPaddedArray = new double[getFFTSize()];
        System.arraycopy(x, 0, zeroPaddedArray, 0, x.length);
        return zeroPaddedArray;
    }

    private void performForwardFFTInPlace(final double[] xDoubleArray) {
        /* perform real forward FFT; automatically zero pads if fftSize > dimension */
        getFFTFactory().realForward(xDoubleArray);
    }

    private void performInverseFFTInPlace(final double[] xDoubleArray) {
        getFFTFactory().realInverse(xDoubleArray, true);
    }

    private void applyFilterInPlace(final double[] xDoubleArray) {
        /* apply filter */
        final double[] orderedFourierFactors = getOrderedFourierFactors();
        final int fftSize = getFFTSize();
        for (int k = 0; k < fftSize; k++) {
            xDoubleArray[k] *= orderedFourierFactors[k];
        }
    }
}
