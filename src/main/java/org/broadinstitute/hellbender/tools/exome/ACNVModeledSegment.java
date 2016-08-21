package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.posteriorsummary.PosteriorSummary;

import java.io.Serializable;

/**
 * Represents a segment with copy-ratio--segment-mean and minor-allele-fraction posterior summaries.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class ACNVModeledSegment implements Locatable, Serializable {

    static final long serialVersionUID = 3373373371L;

    private final SimpleInterval interval;
    private final PosteriorSummary segmentMeanPosteriorSummary;
    private final PosteriorSummary minorAlleleFractionPosteriorSummary;

    public ACNVModeledSegment(final SimpleInterval interval,
                              final PosteriorSummary segmentMeanPosteriorSummary,
                              final PosteriorSummary minorAlleleFractionPosteriorSummary) {
        this.interval = interval;
        this.segmentMeanPosteriorSummary = segmentMeanPosteriorSummary;
        this.minorAlleleFractionPosteriorSummary = minorAlleleFractionPosteriorSummary;
    }


    public SimpleInterval getInterval() {
        return interval;
    }

    public PosteriorSummary getSegmentMeanPosteriorSummary() {
        return segmentMeanPosteriorSummary;
    }

    public PosteriorSummary getMinorAlleleFractionPosteriorSummary() {
        return minorAlleleFractionPosteriorSummary;
    }

    public double getSegmentMeanInCRSpace() {
        return Math.pow(2, segmentMeanPosteriorSummary.getCenter());
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
