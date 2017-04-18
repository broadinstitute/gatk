package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jApacheAdapterUtils;
import org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils.Nd4jIOUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.HMM;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.File;
import java.util.Arrays;

/**
 * This class represents the per-base prior transition probability from one integer copy number
 * state to another.
 *
 * Notes:
 *
 *   - The columns and rows of the transition matrix represent "FROM" and "TO" states, respectively.
 *
 *   - The transition matrix does not need to be normalized (i.e. each column of the matrix sum to 1)
 *     and it will be automatically enforced by the constructor (along with a warning of the deviation
 *     of the probability lies outside of a {@link #PROBABILITY_NORMALIZATION_TOL} neighborhood of 1.
 *     Nevertheless, the user must keep in mind the row/column semantics of the transition matrix
 *     when generating input files.
 *
 *   - The user has the option to "pad" the transition matrix with extra states with ascending copy
 *     numbers provided that {@link #padding} > 0. The necessity for state padding is a technical one:
 *     since the interface {@link HMM} requires the implementation of a single method
 *     {@link HMM#hiddenStates()}, all targets share the same list of hidden states.
 *     Therefore, even though certain contigs require fewer states, the set of states must be appropriately
 *     padded to bring up the maximum copy number state on all contigs to the same level.
 *
 *   - In practice, the user does not need to pad the states manually: a collection of prior
 *     transition matrices can be loaded using {@link IntegerCopyNumberTransitionProbabilityCacheCollection},
 *     and the latter provides the automatic padding functionality with the minimum number of
 *
 * @implNote This class is intended to be immutable.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionMatrix {

    private static final Logger logger = LogManager.getLogger(IntegerCopyNumberTransitionMatrix.class);

    private static final double PROBABILITY_NORMALIZATION_TOL = 1e-6;

    /**
     * Maximum copy number allowed by this prior
     */
    private final int maxCopyNumber;

    /**
     * Padding with extra hidden states
     */
    private final int padding;

    /**
     * The transition matrix
     */
    private final RealMatrix transitionMatrix;

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
    public IntegerCopyNumberTransitionMatrix(@Nonnull final RealMatrix transitionMatrix,
                                             final int padding) {
        Utils.nonNull(transitionMatrix, "The transition matrix must be non-null");
        Utils.validateArg(transitionMatrix.isSquare(), "The transition matrix must be square");
        ParamUtils.isPositiveOrZero(padding, "The padding value must be a non-negative integer");
        this.transitionMatrix = getPaddedTransitionMatrix(transitionMatrix, padding);
        this.maxCopyNumber = this.transitionMatrix.getColumnDimension() - 1;
        this.padding = padding;
        enforceTransitionMatrixValidity();
    }

    /**
     * Reads a single integer copy number transition matrix from a tab-separated table written in the style of
     * {@link Nd4jIOUtils#writeNDArrayMatrixToTextFile}
     *
     * @param inputFile the input tab-separated file
     * @throws UserException.CouldNotReadInputFile if the input file could not be read
     * @throws UserException.BadInput if the table is malformed or the matrix is not square
     * @return an instance of {@link IntegerCopyNumberTransitionMatrix}
     */
    public static IntegerCopyNumberTransitionMatrix read(@Nonnull final File inputFile, final int padding) {
        final RealMatrix transitionMatrix = Nd4jApacheAdapterUtils.convertINDArrayToApacheMatrix(
                Nd4jIOUtils.readNDArrayMatrixFromTextFile(inputFile));
        return new IntegerCopyNumberTransitionMatrix(transitionMatrix, padding);
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
     * @param originalTransitionMatrix the original non-padded transition matrix
     * @param padding non-negative number of extra hidden states to pad
     * @return an instance of {@link IntegerCopyNumberTransitionMatrix}
     */
    private static RealMatrix getPaddedTransitionMatrix(final RealMatrix originalTransitionMatrix, final int padding) {
        if (padding > 0) {
            final int maxCopyNumber = originalTransitionMatrix.getColumnDimension() - 1;
            final int newMaxCopyNumber = maxCopyNumber + padding;
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
