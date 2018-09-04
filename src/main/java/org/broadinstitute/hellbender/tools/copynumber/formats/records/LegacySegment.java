package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

public class LegacySegment implements Locatable {
    private final String sampleName;
    private final SimpleInterval interval;
    private final int numProbes;
    private final double segmentMean;

    public LegacySegment(final String sampleName,
                         final SimpleInterval interval,
                         final int numProbes,
                         final double segmentMean) {
        Utils.nonEmpty(sampleName);
        Utils.nonNull(interval);
        ParamUtils.isPositiveOrZero(numProbes, "Number of probes must be non-negative.");
        this.sampleName = sampleName;
        this.interval = interval;
        this.numProbes = numProbes;
        this.segmentMean = segmentMean;
    }

    public String getSampleName() {
        return sampleName;
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

    public SimpleInterval getInterval() {
        return interval;
    }

    public int getNumProbes() {
        return numProbes;
    }

    public double getSegmentMean() {
        return segmentMean;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final LegacySegment that = (LegacySegment) o;
        return numProbes == that.numProbes &&
                Double.compare(that.segmentMean, segmentMean) == 0 &&
                sampleName.equals(that.sampleName) &&
                interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = sampleName.hashCode();
        result = 31 * result + interval.hashCode();
        result = 31 * result + numProbes;
        temp = Double.doubleToLongBits(segmentMean);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "LegacySegment{" +
                "sampleName='" + sampleName + '\'' +
                ", interval=" + interval +
                ", numProbes=" + numProbes +
                ", segmentMean=" + segmentMean +
                '}';
    }
}
