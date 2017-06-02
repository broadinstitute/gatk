package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Structure to hold the result of the PoN reduction steps.
 */
final class ReductionResult {
    private final RealMatrix pseudoInverse;
    private final RealMatrix reducedCounts;
    private final RealMatrix reducedPseudoInverse;
    private final double[] allSingularValues;

    ReductionResult(final RealMatrix pseudoInverse, final RealMatrix reducedCounts, final RealMatrix reducedPseudoInverse, final double[] allSingularValues) {
        this.pseudoInverse = pseudoInverse;
        this.reducedCounts = reducedCounts;
        this.reducedPseudoInverse = reducedPseudoInverse;
        this.allSingularValues = allSingularValues;
    }

    RealMatrix getPseudoInverse() {
        return pseudoInverse;
    }

    RealMatrix getReducedCounts() {
        return reducedCounts;
    }

    RealMatrix getReducedPseudoInverse() {
        return reducedPseudoInverse;
    }

    double[] getAllSingularValues() {
        return allSingularValues;
    }
}
