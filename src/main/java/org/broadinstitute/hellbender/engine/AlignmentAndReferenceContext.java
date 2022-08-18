package org.broadinstitute.hellbender.engine;

/**
 * Bundles together and AlignmentContext and a ReferenceContext
 */
public class AlignmentAndReferenceContext {

    private final AlignmentContext alignmentContext;
    private final ReferenceContext referenceContext;
    private double activityScore = 0.0;

    public AlignmentAndReferenceContext(final AlignmentContext alignmentContext,
                                        final ReferenceContext referenceContext) {
        this.alignmentContext = alignmentContext;
        this.referenceContext = referenceContext;
    }

    /**
     * getter for the AlignmentContect
     * @return the alignmentContext
     */
    public AlignmentContext getAlignmentContext() {
        return alignmentContext;
    }

    /**
     * getter for the ReferenceContect
     * @return the referenceContext
     */
    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public double getActivityScore() {
        return activityScore;
    }

    public void setActivityScore(double activityScore) {
        this.activityScore = activityScore;
    }
}
