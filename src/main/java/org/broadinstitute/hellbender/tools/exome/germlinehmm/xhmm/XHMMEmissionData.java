package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

/**
 * This class represents the data required to calculate the emission probability at a
 * given target according to the XHMM model.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class XHMMEmissionData {

    private final double coverageZScore;

    public XHMMEmissionData(final double coverageZScore) {
        this.coverageZScore = coverageZScore;
    }

    public double getCoverageZScore() {
        return coverageZScore;
    }
}
