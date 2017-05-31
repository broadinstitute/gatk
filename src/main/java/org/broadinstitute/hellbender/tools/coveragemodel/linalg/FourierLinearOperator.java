package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.jtransforms.fft.DoubleFFT_1D;

import javax.annotation.Nonnull;

import static com.google.common.math.IntMath.isPowerOfTwo;
import static org.broadinstitute.hellbender.utils.GATKProtectedMathUtils.nearestNeighborUniform1DInterpolate;
import static org.broadinstitute.hellbender.utils.GATKProtectedMathUtils.smallestPowerOfTwoGreaterThan;

/**
 * An abstract class that represents a circulant linear operator $F(x,x') = F(|x-x'|)$ specified
 * by the DFT components of its first row (or column) $F(x)$, i.e. {@link #fourierFactors}[k] = $DFT_k[F(x)]$.
 *
 * The action of the operator on a vector $v(x)$ is defined as:
 *
 *      F[v](x) = FFT^{-1}[\sum_k F(k) FFT[v](k)](x)
 *
 * The FFT operations are implemented using JTransforms.
 *
 * Inheriting classes must implement {@link GeneralLinearOperator#operate(Object)}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class FourierLinearOperator<VECTOR> implements GeneralLinearOperator<VECTOR> {

    private final int dimension;
    private final int fftSize;
    private final DoubleFFT_1D fftFactory;
    private final double[] fourierFactors;

    /**
     * The properly ordered physical layout of the Fourier factors as required by
     * {@link DoubleFFT_1D#realForward(double[])}
     */
    private final double[] orderedFourierFactors;

    /**
     * Constructs the linear operator of a given {@code dimension} for given Fourier factors
     * Note: the length of {@code fourierFactors} must be equal to {@code floor(dimension/2) + 1}
     *
     * @param dimension dimension of the linear space on which the linear operator acts
     * @param fourierFactors DFT components of the first column of the operator in real space
     * @param promoteToPowerOfTwo promote the dimension to the power of two
     * @throws IllegalArgumentException if {@code fourierFactors.length != dimension/2 + 1}
     */
    public FourierLinearOperator(final int dimension,
                                 @Nonnull final double[] fourierFactors,
                                 final boolean promoteToPowerOfTwo) {
        ParamUtils.isPositive(dimension - 1, "The dimension of the linear operator " +
                "must be >= 2.");
        Utils.validateArg(fourierFactors.length == dimension/2 + 1, "The length of Fourier factors must be equal" +
                " to floor(dimension/2) + 1");
        if (!promoteToPowerOfTwo || isPowerOfTwo(dimension)) {
            /* dimension is already power of 2, turn off the flag */
            this.fftSize = dimension;
            this.dimension = dimension;
            this.fourierFactors = fourierFactors.clone();
        } else {
            this.fftSize = smallestPowerOfTwoGreaterThan(dimension);
            this.dimension = dimension;
            this.fourierFactors = nearestNeighborUniform1DInterpolate(fourierFactors, fftSize/2 + 1);
        }
        fftFactory = new DoubleFFT_1D(fftSize);
        orderedFourierFactors = new double[fftSize];
        initializeOrderedFourierFactors();
    }

    public double[] getOrderedFourierFactors() {
        return orderedFourierFactors.clone();
    }

    private void initializeOrderedFourierFactors() {
        final double[] fourierFactors = getFourierFactors();
        final int fftSize = getFFTSize();
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

    @Override
    public int getRowDimension() { return dimension; }

    @Override
    public int getColumnDimension() { return dimension; }

    @Override
    public boolean isTransposable() { return true; }

    public int getFFTSize() { return fftSize; }

    @Override
    public VECTOR operateTranspose(@Nonnull final VECTOR x)
            throws DimensionMismatchException, UnsupportedOperationException {
        return operate(x);
    }

    public double[] getFourierFactors() {
        return fourierFactors.clone();
    }

    public int getDimension() {
        return dimension;
    }

    public DoubleFFT_1D getFFTFactory() {
        return fftFactory;
    }
}
