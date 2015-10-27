package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Simple implementation of the SVD interface for storing the matrices (and vector) of a SVD result.
 */
public class SimpleSVD implements SVD {
    private RealMatrix v;
    private RealMatrix u;
    private double [] singularValues;
    private RealMatrix pinv;


    public SimpleSVD(final RealMatrix U, double [] singularValues, final RealMatrix V, final RealMatrix pinv) {
        this.v = V;
        this.u = U;
        this.singularValues = singularValues;
        this.pinv = pinv;
    }

    @Override
    public RealMatrix getV() {
        return v;
    }

    @Override
    public RealMatrix getU() {
        return u;
    }

    @Override
    public double[] getSingularValues() {
        return singularValues;
    }

    @Override
    public RealMatrix getPinv() {
        return pinv;
    }
}
