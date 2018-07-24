package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

public class LegacyCopyRatioSegment extends CopyRatioSegment {
    final protected String sampleName;

    public LegacyCopyRatioSegment(final String sampleName, final SimpleInterval interval, final int numPoints, final double meanLog2CopyRatio) {
        super(interval, numPoints, meanLog2CopyRatio);
        this.sampleName = sampleName;
    }

    public LegacyCopyRatioSegment(final String sampleName, final SimpleInterval interval, final List<CopyRatio> denoisedLog2CopyRatios) {
        super(interval, denoisedLog2CopyRatios);
        this.sampleName = sampleName;
    }

    public String getSampleName() {
        return sampleName;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        LegacyCopyRatioSegment that = (LegacyCopyRatioSegment) o;

        return sampleName != null ? sampleName.equals(that.sampleName) : that.sampleName == null;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (sampleName != null ? sampleName.hashCode() : 0);
        return result;
    }
}
