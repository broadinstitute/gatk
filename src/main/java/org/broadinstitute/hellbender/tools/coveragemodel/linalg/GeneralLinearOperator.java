package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;

import javax.annotation.Nonnull;

/**
 * An abstract linear operator
 *
 * @param <VECTOR> a vector type

 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public interface GeneralLinearOperator<VECTOR> {

    /**
     * Returns the dimension of the codomain of this operator.
     *
     * @return the number of rows of the underlying matrix
     */
    int getRowDimension();

    /**
     * Returns the dimension of the domain of this operator.
     *
     * @return the number of columns of the underlying matrix
     */
    int getColumnDimension();

    /**
     * Returns the result of multiplying {@code this} by the vector {@code x}.
     *
     * @param x the vector to operate on
     * @return the product of {@code this} instance with {@code x}
     * @throws DimensionMismatchException if the column dimension does not match
     * the size of {@code x}
     */
    VECTOR operate(@Nonnull final VECTOR x)
            throws DimensionMismatchException;

    /**
     * Returns the result of multiplying the transpose of {@code this} operator
     * by the vector {@code x} (optional operation). The default implementation
     * throws an {@link UnsupportedOperationException}. Users overriding this
     * method must also override {@link #isTransposable()}.
     *
     * @param x the vector to operate on
     * @return the product of the transpose of {@code this} instance with
     * {@code x}
     * @throws DimensionMismatchException
     * if the row dimension does not match the size of {@code x}
     * @throws UnsupportedOperationException if this operation is not supported
     * by {@code this} operator
     */
    default VECTOR operateTranspose(@Nonnull final VECTOR x)
            throws DimensionMismatchException, UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns {@code true} if this operator supports
     * {@link #operateTranspose(VECTOR)}. If {@code true} is returned,
     * {@link #operateTranspose(VECTOR)} should not throw
     * {@code UnsupportedOperationException}. The default implementation returns
     * {@code false}.
     *
     * @return {@code false}
     */
    default boolean isTransposable() {
        return false;
    }

    /**
     * Updates a row of the matrix representation of the linear operator
     * @param rowIndex row index
     * @param row new row
     * @throws DimensionMismatchException
     * @throws UnsupportedOperationException
     */
    default void updateRow(final int rowIndex, @Nonnull final VECTOR row)
            throws DimensionMismatchException, UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    /**
     * Updates a column of the matrix representation of the linear operator
     * @param colIndex column index
     * @param col new column
     * @throws DimensionMismatchException
     * @throws UnsupportedOperationException
     */
    default void updateColumn(final int colIndex, @Nonnull final VECTOR col)
            throws DimensionMismatchException, UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    /**
     * Performs clean up after a call to {@link #operate(Object)} such as null-referencing
     * unnecessary members, unpersisting RDD caches, etc)
     */
    default void cleanupAfter() {}
}
