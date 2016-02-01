package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

/**
 * Represents a segment with copy-ratio--segment-mean and minor-allele-fraction posterior summaries.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class ACNVModeledSegment extends ModeledSegment {

    private final PosteriorSummary segmentMeanPosteriorSummary;
    private final PosteriorSummary minorAlleleFractionPosteriorSummary;

    public ACNVModeledSegment(final SimpleInterval interval,
                              final PosteriorSummary segmentMeanPosteriorSummary,
                              final PosteriorSummary minorAlleleFractionPosteriorSummary) {
        super(interval, ModeledSegment.NO_CALL, 0, segmentMeanPosteriorSummary.center());
        this.segmentMeanPosteriorSummary = segmentMeanPosteriorSummary;
        this.minorAlleleFractionPosteriorSummary = minorAlleleFractionPosteriorSummary;
    }


    /**
     *  Get segment mean in log2 space
     * @return estimate of the mean for this modeled segment
     */
    @Override
    public double getSegmentMean() {
        return segmentMeanPosteriorSummary.center();
    }

    /**
     * Get the segment mean in non-logged space
     *
     * @return estimate of the mean for this modeled segment
     */
    @Override
    public double getSegmentMeanInCRSpace() {
        return Math.pow(2, segmentMeanPosteriorSummary.center());
    }

    /**
     * Not supported in this class
     *
     * @param segmentMean desired log2 mean
     */
    @Override
    public void setSegmentMean(final double segmentMean) {
        throw new UnsupportedOperationException("Cannot set the segment mean for ACNVModeledSegment.  Set the segmentMeanPosteriorSummary, instead.");
    }

    /**
     * Not supported in this class
     *
     * @param segmentMean desired mean
     */
    @Override
    public void setSegmentMeanInCRSpace(final double segmentMean) {
        throw new UnsupportedOperationException("Cannot set the segment mean for ACNVModeledSegment.  Set the segmentMeanPosteriorSummary, instead.");
    }


    public PosteriorSummary getSegmentMeanPosteriorSummary() {
        return segmentMeanPosteriorSummary;
    }

    public PosteriorSummary getMinorAlleleFractionPosteriorSummary() {
        return minorAlleleFractionPosteriorSummary;
    }
}
