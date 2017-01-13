package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.jtransforms.fft.DoubleFFT_1D;

import javax.annotation.Nonnull;
import java.util.stream.IntStream;

import static com.google.common.math.IntMath.isPowerOfTwo;
import static org.broadinstitute.hellbender.utils.GATKProtectedMathUtils.nearestNeighborUniform1DInterpolate;
import static org.broadinstitute.hellbender.utils.GATKProtectedMathUtils.smallestPowerOfTwoGreaterThan;

/**
 * A circulant linear operator specified by its eigenvalues in the Fourier space.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class FourierLinearOperator<V> extends GeneralLinearOperator<V> {

    protected final int dimension;
    protected final int fftSize;
    protected final double[] fourierFactors;
    protected final DoubleFFT_1D fftFactory;

    protected final boolean promoteToPowerOfTwo;

    /**
     * Constructs the linear operator of a given {@code dimension} for given Fourier factors
     * Note: the length of {@code fourierFactors} must be equal to {@code floor(dimension/2) + 1}
     *
     * @param dimension dimension of the linear space on which the linear operator acts
     * @param fourierFactors Fourier factors of the linear operator
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
            this.promoteToPowerOfTwo = false;
            this.fftSize = dimension;
            this.dimension = dimension;
            this.fourierFactors = fourierFactors.clone();
        } else {
            this.promoteToPowerOfTwo = true;
            this.fftSize = smallestPowerOfTwoGreaterThan(dimension);
            this.dimension = dimension;
            this.fourierFactors = nearestNeighborUniform1DInterpolate(fourierFactors, fftSize/2 + 1);
        }
        fftFactory = new DoubleFFT_1D(fftSize);
    }

    @Override
    public int getRowDimension() { return dimension; }

    @Override
    public int getColumnDimension() { return dimension; }

    @Override
    public boolean isTransposable() {
        return true;
    }

    public int getFFTSize() { return fftSize; }

    @Override
    public V operateTranspose(@Nonnull final V x)
            throws DimensionMismatchException, UnsupportedOperationException {
        return operate(x);
    }

    public double[] getFourierFactors() {
        return fourierFactors.clone();
    }

    /**
     * Generate Fourier factors for a midpass filter
     *
     * @param dimension dimension of the signal
     * @param freqLowerCutoff lower frequency cutoff (excluded)
     * @param freqUpperCutoff upper frequency cutoff (included)
     * @return Fourier factors
     */
    public static double[] getMidpassFilterFourierFactors(final int dimension, final int freqLowerCutoff,
                                                          final int freqUpperCutoff) {
        ParamUtils.isPositiveOrZero(freqUpperCutoff - freqLowerCutoff, "Upper cutoff must be >= lower cutoff.");
        ParamUtils.isPositive(dimension/2 + 1 - freqUpperCutoff, "For dimension " + dimension +
                ", the upper frequency cutoff must be <= " + (dimension/2));
        return IntStream.range(0, dimension/2 + 1)
                .mapToDouble(i -> i > freqLowerCutoff && i <= freqUpperCutoff ? 1.0 : 0).toArray();
    }
}
