package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;

/**
 *  Perform singular value decomposition (and pseudoinverse calculation).
 */
public interface SingularValueDecomposer {
    public SVD createSVD(final RealMatrix m);
}
