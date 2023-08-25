package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;

// Mini class to package SV type and interval into one object
public class SVSegment implements Locatable {
    protected final GATKSVVCFConstants.StructuralVariantAnnotationType intervalSVType;
    protected final SimpleInterval interval;

    public SVSegment(final GATKSVVCFConstants.StructuralVariantAnnotationType svType, final SimpleInterval interval) {
        Utils.nonNull(interval);
        this.intervalSVType = svType;
        this.interval = interval;
    }

    public GATKSVVCFConstants.StructuralVariantAnnotationType getIntervalSVType() {
        return intervalSVType;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final SVSegment svSegment = (SVSegment) o;
        return intervalSVType == svSegment.intervalSVType && interval.equals(svSegment.interval);
    }

    @Override
    public int hashCode() {
        return Objects.hash(intervalSVType, interval);
    }
}
