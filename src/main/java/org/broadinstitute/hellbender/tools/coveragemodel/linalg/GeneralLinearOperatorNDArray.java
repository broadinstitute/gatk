package org.broadinstitute.hellbender.tools.coveragemodel.linalg;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.broadinstitute.hellbender.utils.Utils;
import org.nd4j.linalg.api.ndarray.INDArray;

import javax.annotation.Nonnull;

/**
 * A general linear operator specified by an {@link INDArray}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class GeneralLinearOperatorNDArray implements GeneralLinearOperator<INDArray> {

    private final INDArray mat;

    public GeneralLinearOperatorNDArray(@Nonnull final INDArray mat) {
        Utils.validateArg(mat.rank() == 2, "The provided INDArray must be rank-2.");
        this.mat = mat;
    }

    @Override
    public int getRowDimension() {
        return mat.shape()[0];
    }

    @Override
    public int getColumnDimension() {
        return mat.shape()[1];
    }

    @Override
    public INDArray operate(@Nonnull final INDArray x) throws DimensionMismatchException {
        return mat.mmul(x);
    }

    @Override
    public INDArray operateTranspose(@Nonnull final  INDArray x)
            throws DimensionMismatchException, UnsupportedOperationException {
        return mat.transpose().mmul(x);
    }

    @Override
    public boolean isTransposable() {
        return true;
    }
}
