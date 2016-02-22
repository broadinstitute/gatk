package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Coverage units supported.
 */
public enum CoverageUnit {

    /**
     * Each qualifying read contributes 1 to the coverage.
     */
    OVERLAPPING_READ(CoverageDatum.READ, CoverageTransformation.NONE, false),

    /**
     * Each qualifying fragment contributes 1 to the coverage.
     * <p>
     *     Qualifying fragments are the ones that have at least one qualifying read.
     * </p>
     */
    OVERLAPPING_FRAGMENT(CoverageDatum.READ, CoverageTransformation.NONE, true),

    /**
     * Average number of reads stacked up across the target.
     * <p>
     *  Notice that deletions in the read does not contribute to the depth on the overlapping site.
     * </p>
     */
    AVERAGE_READ_DEPTH(CoverageDatum.BASE_CALL, CoverageTransformation.AVERAGE_PER_BP, false),

    /**
     * Average number of fragments stacked up across the target.
     * <p>
     *  Notice that deletions in reads belonging to the target or gaps between read pairs are not counted.
     * </p>
     */
    AVERAGE_FRAGMENT_DEPTH(CoverageDatum.BASE_CALL, CoverageTransformation.AVERAGE_PER_BP, true);

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
     * Whether it counts evidences based on fragments rather than reads.
     */
    final boolean basedOnFragments;

    /**
     * Creates a new unit given the relevant datum and transformation.
     * @param datum the datum type the unit is based on.
     * @param transformation the transformation that need to be performed to the absolute count.
     */
    CoverageUnit(final CoverageDatum datum, final CoverageTransformation transformation, final boolean basedOnFragments) {
        this.datum = Utils.nonNull(datum);
        this.transformation = Utils.nonNull(transformation);
        this.basedOnFragments = basedOnFragments;
    }
}

