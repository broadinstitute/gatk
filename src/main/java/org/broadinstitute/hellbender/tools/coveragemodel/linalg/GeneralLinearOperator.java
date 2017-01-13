package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;

import javax.annotation.Nonnull;

/**
 * An abstract linear operator
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class GeneralLinearOperator<V> {

    /**
     * Returns the dimension of the codomain of this operator.
     *
     * @return the number of rows of the underlying matrix
     */
    public abstract int getRowDimension();

    /**
     * Returns the dimension of the domain of this operator.
     *
     * @return the number of columns of the underlying matrix
     */
    public abstract int getColumnDimension();

    /**
     * Returns the result of multiplying {@code this} by the vector {@code x}.
     *
     * @param x the vector to operate on
     * @return the product of {@code this} instance with {@code x}
     * @throws DimensionMismatchException if the column dimension does not match
     * the size of {@code x}
     */
    public abstract V operate(@Nonnull final V x)
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
     * @throws org.apache.commons.math3.exception.DimensionMismatchException
     * if the row dimension does not match the size of {@code x}
     * @throws UnsupportedOperationException if this operation is not supported
     * by {@code this} operator
     */
    public V operateTranspose(@Nonnull final V x)
            throws DimensionMismatchException, UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns {@code true} if this operator supports
     * {@link #operateTranspose(V)}. If {@code true} is returned,
     * {@link #operateTranspose(V)} should not throw
     * {@code UnsupportedOperationException}. The default implementation returns
     * {@code false}.
     *
     * @return {@code false}
     */
    public boolean isTransposable() {
        return false;
    }

    /**
     * Updates a row of the matrix representation of the linear operator
     * @param rowIndex row index
     * @param row new row
     * @throws DimensionMismatchException
     * @throws UnsupportedOperationException
     */
    public void updateRow(final int rowIndex, @Nonnull final V row)
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
    public void updateColumn(final int colIndex, @Nonnull final V col)
            throws DimensionMismatchException, UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    public void cleanupAfter() {}

}
