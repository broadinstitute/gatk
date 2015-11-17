package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Coverage units supported.
 */
public enum CoverageUnit {

    /**
     * Each qualifying read contributes 1 to the coverage.
     */
    OVERLAPPING_READ(CoverageDatum.READ, CoverageTransformation.NONE),

    /**
     * Average number of reads stacked up across the target.
     * <p>
     *  Notice that deletions in the read does not contribute to the depth on the overlapping site.
     * </p>
     */
    AVERAGE_DEPTH(CoverageDatum.BASE_CALL, CoverageTransformation.AVERAGE_PER_BP);

    /**
     * Reference to the corresponding datum type.
     */
    final CoverageDatum datum;

    /**
     * Reference to the transformation to be carried out on the absolute data count to
     * obtain the actual coverage value to report.
     */
    final CoverageTransformation transformation;

    /**
     * Creates a new unit given the relevant datum and transformation.
     * @param datum the datum type the unit is based on.
     * @param transformation the transformation that need to be performed to the absolute count.
     */
    CoverageUnit(final CoverageDatum datum, final CoverageTransformation transformation) {
        this.datum = Utils.nonNull(datum);
        this.transformation = Utils.nonNull(transformation);
    }
}

