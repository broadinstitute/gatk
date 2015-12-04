package org.broadinstitute.hellbender.tools.exome;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Structure to hold the result of the SVD/reduction steps.
 */
public final class ReductionResult {
    final private RealMatrix pseudoInverse;
    final private RealMatrix reducedCounts;
    final private RealMatrix reducedInverse;
    final private double[] allSingularValues;

    public ReductionResult(final RealMatrix pseudoInverse, final RealMatrix reducedCounts, final RealMatrix reduceInverse, final double[] allSingularValues) {
        this.pseudoInverse = pseudoInverse;
        this.reducedCounts = reducedCounts;
        this.reducedInverse = reduceInverse;
        this.allSingularValues = allSingularValues;
    }

    public RealMatrix getPseudoInverse() {
        return pseudoInverse;
    }

    public RealMatrix getReducedCounts() {
        return reducedCounts;
    }

    public RealMatrix getReducedInverse() {
        return reducedInverse;
    }

    public double[] getAllSingularValues() {
        return allSingularValues;
    }
}
