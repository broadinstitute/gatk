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

    public SimpleSVD(final RealMatrix u, double [] singularValues, final RealMatrix v, final RealMatrix pinv) {
        this.v = v;
        this.u = u;
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
