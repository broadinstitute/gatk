package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.Objects;

public class LocatableCopyNumberPosteriorDistribution extends CopyNumberPosteriorDistribution implements Locatable {
    final SimpleInterval interval;
    final String sample;

    public LocatableCopyNumberPosteriorDistribution(final Map<IntegerCopyNumberState, Double> copyNumberPosteriorDistribution,
                                                    final String sample,
                                                    final SimpleInterval interval) {
        super(copyNumberPosteriorDistribution);
        Utils.nonNull(sample, "Sample cannot be null");
        Utils.nonNull(interval, "Interval cannot be null");
        this.sample = sample;
        this.interval = interval;
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

    public String getSample() {
        return sample;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof LocatableCopyNumberPosteriorDistribution)) return false;
        if (!super.equals(o)) return false;
        LocatableCopyNumberPosteriorDistribution that = (LocatableCopyNumberPosteriorDistribution) o;
        return interval.equals(that.interval) &&
                sample.equals(that.sample);
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), interval, sample);
    }
}
