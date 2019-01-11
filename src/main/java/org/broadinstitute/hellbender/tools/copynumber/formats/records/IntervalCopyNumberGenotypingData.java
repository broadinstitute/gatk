package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * The bundle of integer copy-number posterior distribution and baseline integer copy-number state
 * for an interval.
 *
 * Note: the copy-number posterior distribution is stored in the natural log space.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntervalCopyNumberGenotypingData implements Locatable {
    private final SimpleInterval interval;
    private final CopyNumberPosteriorDistribution copyNumberPosteriorDistribution;
    private final IntegerCopyNumberState baselineCopyNumberState;

    public IntervalCopyNumberGenotypingData(final SimpleInterval interval,
                                            final CopyNumberPosteriorDistribution copyNumberStatePosteriors,
                                            final IntegerCopyNumberState baselineCopyNumberState) {
        this.interval = Utils.nonNull(interval);
        this.copyNumberPosteriorDistribution = Utils.nonNull(copyNumberStatePosteriors);
        this.baselineCopyNumberState = Utils.nonNull(baselineCopyNumberState);
    }

    public CopyNumberPosteriorDistribution getCopyNumberPosteriorDistribution() {
        return copyNumberPosteriorDistribution;
    }

    public IntegerCopyNumberState getBaselineIntegerCopyNumberState() {
        return baselineCopyNumberState;
    }

    public SimpleInterval getInterval() {
        return interval;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof IntervalCopyNumberGenotypingData)) {
            return false;
        }

        final IntervalCopyNumberGenotypingData that = (IntervalCopyNumberGenotypingData) o;
        return that.interval.equals(this.interval) &&
                that.copyNumberPosteriorDistribution.equals(this.copyNumberPosteriorDistribution) &&
                that.baselineCopyNumberState.equals(this.baselineCopyNumberState);
    }

    @Override
    public int hashCode() {
        int result = interval != null ? interval.hashCode() : 0;
        result = 31 * result + (copyNumberPosteriorDistribution != null ? copyNumberPosteriorDistribution.hashCode() : 0);
        result = 31 * result + (baselineCopyNumberState != null ? baselineCopyNumberState.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "IntervalCopyNumberGenotypingData{" +
                "interval=" + interval +
                ", copyNumberPosteriorDistribution=" + copyNumberPosteriorDistribution +
                ", baselineCopyNumberState=" + baselineCopyNumberState +
                '}';
    }
}
