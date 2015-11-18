package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Perform singular value decompoisition (and pseudoinverse calculation) in pure Java, Commons Math.
 */
class ApacheSingularValueDecomposer {

    /** Create a SVD instance using Apache Commons Math.
     *
     * @param m matrix that is not {@code null}
     * @return SVD instance that is never {@code null}
     */
    public static SVD createSVD(final RealMatrix m) {

        Utils.nonNull(m, "Cannot create SVD on a null matrix.");

        final SingularValueDecomposition svd = new SingularValueDecomposition(m);
        final RealMatrix pinv = svd.getSolver().getInverse();
        return new SimpleSVD(svd.getU(), svd.getSingularValues(), svd.getV(), pinv);
    }
}
