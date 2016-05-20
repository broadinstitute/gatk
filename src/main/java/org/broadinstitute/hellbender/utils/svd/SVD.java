package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Interface for SVD implementation.  Leverages Apache Commons for data fields.
 */
public interface SVD {

    /**
     * Get the V matrix of a Singular Value Decomposition
     */
    RealMatrix getV();

    /**
     * Get the U matrix of a Singular Value Decomposition
     */
    RealMatrix getU();

    /**
     * Get the pseudoinverse as calculated using the SVD
     */
    RealMatrix getPinv();

    /**
     * Get the singular values as an array.
     */
    double [] getSingularValues();
}
