package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * A single integer copy-number posterior record for a specific interval.
 *
 * Note: the posterior distribution is stored in the natural log space.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class LocatableCopyNumberPosteriorDistribution implements Locatable {
    private final SimpleInterval interval;
    private final CopyNumberPosteriorDistribution copyNumberPosteriorDistribution;

    public LocatableCopyNumberPosteriorDistribution(final SimpleInterval interval,
                                                    final CopyNumberPosteriorDistribution copyNumberStatePosteriors) {
        this.interval = Utils.nonNull(interval);
        this.copyNumberPosteriorDistribution = Utils.nonNull(copyNumberStatePosteriors);
    }

    /**
     * Get the copy number posteriors for this record
     */
    public CopyNumberPosteriorDistribution getCopyNumberPosteriors() {
        return copyNumberPosteriorDistribution;
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
        if (!(o instanceof LocatableCopyNumberPosteriorDistribution)) {
            return false;
        }

        final LocatableCopyNumberPosteriorDistribution locatableDistribution =
                (LocatableCopyNumberPosteriorDistribution) o;
        return locatableDistribution.interval.equals(this.interval) &&
                locatableDistribution.copyNumberPosteriorDistribution.equals(this.copyNumberPosteriorDistribution);
    }

    @Override
    public int hashCode() {
        return interval.hashCode() + 31 * copyNumberPosteriorDistribution.hashCode();
    }

    @Override
    public String toString() {
        return "LocatableCopyNumberPosteriorDistribution{interval=" + interval.toString() +
                "," + copyNumberPosteriorDistribution.toString() + "}";
    }
}
