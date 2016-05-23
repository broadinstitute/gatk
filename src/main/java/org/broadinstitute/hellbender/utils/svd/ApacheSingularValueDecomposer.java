package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Perform singular value decomposition (and pseudoinverse calculation) in pure Java, Commons Math.
 */
public final class ApacheSingularValueDecomposer implements SingularValueDecomposer{

    /** Create a SVD instance using Apache Commons Math.
     *
     * @param m matrix that is not {@code null}
     * @return SVD instance that is never {@code null}
     */
    @Override
    public SVD createSVD(final RealMatrix m) {

        Utils.nonNull(m, "Cannot create SVD on a null matrix.");

        final SingularValueDecomposition svd = new SingularValueDecomposition(m);
        final RealMatrix pinv = svd.getSolver().getInverse();
        return new SimpleSVD(svd.getU(), svd.getSingularValues(), svd.getV(), pinv);
    }
}
