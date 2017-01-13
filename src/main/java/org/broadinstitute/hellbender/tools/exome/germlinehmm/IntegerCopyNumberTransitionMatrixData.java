package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jApacheAdapterUtils;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.nd4j.linalg.api.ndarray.INDArray;

import javax.annotation.Nonnull;
import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * This class represents the per-base transition prior probability from one integer copy number
 * state to another.
 *
 * Note: we use the convention that the columns and rows represent "FROM" and "TO" states, respectively.
 *
 * The implementation must observe that this class is intended to be immutable.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionMatrixData {

    private final Logger logger = LogManager.getLogger(IntegerCopyNumberTransitionMatrixData.class);

    private static final double PROBABILITY_NORMALIZATION_TOL = 1e-6;

    /**
     * Maximum copy number allowed by this prior
     */
    final private int maxCopyNumber;

    /**
     * Padding with extra hidden states
     */
    final private int padding;

    /**
     * The transition matrix
     */
    final private RealMatrix transitionMatrix;

    /**
     * Public constructor.
     *
     * The rows and columns of the transition matrix denote departure and arrival states,
     * respectively.
     *
     * @param transitionMatrix the transition matrix; must be non-null and square
     * @param padding padding with extra states
     * @throws IllegalArgumentException if the transition matrix is not square
     * @throws IllegalArgumentException if padding is negative
     */
    public IntegerCopyNumberTransitionMatrixData(@Nonnull final RealMatrix transitionMatrix,
                                                 final int padding) {
        this.transitionMatrix = getPaddedTransitionMatrix(transitionMatrix, padding); /* checks are performed here */
        this.maxCopyNumber = this.transitionMatrix.getColumnDimension() - 1;
        this.padding = padding;
        enforceTransitionMatrixValidity();
    }

    /**
     * Reads a single integer copy number transition matrix from a tab-separated table written in the style of
     * {@link org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils#writeNDArrayToTextFile(INDArray, File, List, List)}
     *
     * @param inputFile the input tab-separated file
     * @throws UserException.CouldNotReadInputFile if the input file could not be read
     * @throws UserException.BadInput if the table is malformed or the matrix is not square
     * @return an instance of {@link IntegerCopyNumberTransitionMatrixData}
     */
    public static IntegerCopyNumberTransitionMatrixData read(@Nonnull final File inputFile, final int padding) {
        final RealMatrix transitionMatrix = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(
                Nd4jIOUtils.readNDArrayFromTextFile(inputFile));
        return new IntegerCopyNumberTransitionMatrixData(transitionMatrix, padding);
    }

    /**
     * Return a clone of the transition matrix
     * @return the transition matrix
     */
    public RealMatrix getTransitionMatrix() {
        return transitionMatrix.copy();
    }

    public int getMaxCopyNumber() {
        return maxCopyNumber;
    }

    public int getPadding() {
        return padding;
    }

    private void enforceTransitionMatrixValidity() {
        for (int colIndex = 0; colIndex <= maxCopyNumber; colIndex++) {
            final double[] col = transitionMatrix.getColumn(colIndex);
            if (Arrays.stream(col).anyMatch(d -> d < 0)) {
                throw new UserException.BadInput(String.format("Column %d of the transition matrix has a negative" +
                        " entry. Only zero or positive entries are allowed.", colIndex));
            }
            final double colSum = Arrays.stream(col).reduce(Double::sum).getAsDouble();
            if (FastMath.abs(colSum - 1.0) > PROBABILITY_NORMALIZATION_TOL) {
                logger.warn(String.format("Column %d of the transition matrix is not properly normalized (sum = %f)." +
                        " Enforcing normalization. Please check the input data for consistency.", colIndex, colSum));
                transitionMatrix.setColumn(colIndex, Arrays.stream(col).map(d -> d / colSum).toArray());
            }
        }
    }

    /**
     * Pad the transition matrix with extra hidden states. The old transition matrix will be embedded on the upper left
     * corner of an identity matrix. As a result, there will be no mixing between the existing states and the hidden
     * states.
     *
     * Note: padding is performed for convenience and the extra hidden states are meant to be inaccessible. It is the
     * user's responsibility to make sure that the prior probabilities prohibits the occupancy of these states.
     *
     * @param padding the non-negative number of extra hidden states
     * @throws IllegalArgumentException if the transition matrix is not square
     * @throws IllegalArgumentException if padding is negative
     * @return an instance of {@link IntegerCopyNumberTransitionMatrixData}
     */
    private static RealMatrix getPaddedTransitionMatrix(final RealMatrix originalTransitionMatrix, final int padding) {
        if (!Utils.nonNull(originalTransitionMatrix, "The transition matrix must be non-null").isSquare()) {
            throw new IllegalArgumentException("The transition matrix must be square");
        }
        final int maxCopyNumber = originalTransitionMatrix.getColumnDimension() - 1;
        final int newMaxCopyNumber = maxCopyNumber + ParamUtils.isPositiveOrZero(padding,
                "Padding must be non-negative");
        if (padding > 0) {
            final RealMatrix paddedTransitionMatrix = MatrixUtils.createRealIdentityMatrix(newMaxCopyNumber + 1);
            for (int i = 0; i <= maxCopyNumber; i++) {
                for (int j = 0; j <= maxCopyNumber; j++) {
                    paddedTransitionMatrix.setEntry(i, j, originalTransitionMatrix.getEntry(i, j));
                }
            }
            return paddedTransitionMatrix;
        } else {
            return originalTransitionMatrix;
        }
    }

}
