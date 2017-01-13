package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;

/**
 * A a subclass of {@link RealLinearOperator} that defines the action of a real circulant
 * matrix operator $F(x,x') = F(x-x')$ by providing the DFT components of $F(x)$
 *
 * We use JTransforms for performing FFT. It is the fastest FFT implementation in pure Java,
 * but should also consider trying nd4j's native FFT backend or an FFTW wrapper.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */

public final class FourierLinearOperatorNDArray extends FourierLinearOperator<INDArray> {

    private final double[] orderedFourierFactors;

    public FourierLinearOperatorNDArray(final int dimension, final double[] fourierFactors,
                                        final boolean promoteToPowerOfTwo) {
        super(dimension, fourierFactors, promoteToPowerOfTwo);
        orderedFourierFactors = new double[fftSize];
        initializeOrderedFourierFactors();
    }

    /**
     * Implementation of the action of the Fourier linear operator on a given real vector
     * @param x an {@link INDArray} on which the linear operator is acted on
     * @return the transformed vector
     */
    @Override
    public INDArray operate(@Nonnull final INDArray x) throws DimensionMismatchException {
        if (x.length() != dimension) {
            throw new DimensionMismatchException(x.length(), dimension);
        }
        final double[] xDoubleArray = createDoubleArrayFromINDArray(x);
        performForwardFFTInPlace(xDoubleArray);
        applyFilterInPlace(xDoubleArray);
        performInverseFFTInPlace(xDoubleArray);
        return Nd4j.create(xDoubleArray, new int[]{dimension, 1});
    }

    /**
     * @param x the input vector
     * @return rfft of x
     */
    public double[] getForwardFFT(@Nonnull final INDArray x) {
        if (x.length() != dimension) {
            throw new DimensionMismatchException(x.length(), dimension);
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
        if (xFFT.length != fftSize) {
            throw new DimensionMismatchException(xFFT.length, fftSize);
        }
        final double[] xDoubleArray = xFFT.clone();
        performInverseFFTInPlace(xDoubleArray);
        return Nd4j.create(xDoubleArray, new int[]{dimension, 1});
    }

    /**
     * @param xFFT the input vector
     * @return rfft of x
     */
    public INDArray getInverseFFT(@Nonnull final INDArray xFFT) {
        if (xFFT.length() != fftSize) {
            throw new DimensionMismatchException(xFFT.length(), fftSize);
        }
        final double[] xFFTDoubleArray = new double[fftSize];
        for (int i = 0; i < fftSize; i++) {
            xFFTDoubleArray[i] = xFFT.getDouble(i);
        }
        performInverseFFTInPlace(xFFTDoubleArray);
        return Nd4j.create(xFFTDoubleArray, new int[]{dimension, 1});
    }

    public double[] getOrderedFourierFactors() {
        return orderedFourierFactors.clone();
    }

    private void initializeOrderedFourierFactors() {
        orderedFourierFactors[0] = fourierFactors[0];
        orderedFourierFactors[1] = fourierFactors[fftSize/2];
        for (int k = 1; k < fftSize/2; k++) {
            orderedFourierFactors[2*k] = fourierFactors[k];
            orderedFourierFactors[2*k+1] = fourierFactors[k];
        }
        if (fftSize % 2 == 1) {
            orderedFourierFactors[fftSize-1] = fourierFactors[(fftSize-1)/2];
        }
    }

    private double[] createDoubleArrayFromINDArray(@Nonnull final INDArray x) {
        /* transform to real space; entries > dimension will be ignored */
        final double[] xDoubleArray = new double[fftSize];
        for (int i = 0; i < dimension; i++) {
            xDoubleArray[i] = x.getDouble(i);
        }
        return xDoubleArray;
    }

    private void performForwardFFTInPlace(final double[] xDoubleArray) {
        /* perform real forward FFT; automatically zero pads if fftSize > dimension */
        fftFactory.realForward(xDoubleArray);
    }

    private void performInverseFFTInPlace(final double[] xDoubleArray) {
        fftFactory.realInverse(xDoubleArray, true);
    }

    private void applyFilterInPlace(final double[] xDoubleArray) {
        /* apply filter */
        for (int k = 0; k < fftSize; k++) {
            xDoubleArray[k] *= orderedFourierFactors[k];
        }
    }
}
