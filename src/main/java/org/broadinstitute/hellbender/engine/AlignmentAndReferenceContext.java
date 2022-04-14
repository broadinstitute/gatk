package org.broadinstitute.hellbender.engine;

/**
 * Bundles together and AlignmentContext and a ReferenceContext
 */
public class AlignmentAndReferenceContext {

    private final AlignmentContext alignmentContext;
    private final ReferenceContext referenceContext;

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
}
