package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;

// Mini class to package SV type and interval into one object
public class SVSegment implements Locatable, Comparable<SVSegment> {
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

    /**
     * Compares to another SVSegment using contig (alphabetical), then start position, then end position, then SV type
     * @param that SVSegment
     * @return 0 if first and second are equal, a negative value if first < second, or a positive value if first > second
     */
    public int compareTo(final SVSegment that) {
        Utils.nonNull(that);

        int result = this.getContig().compareTo(that.getContig());
        if ( result == 0 ) {
            result = Integer.compare(this.getStart(), that.getStart());
            if ( result == 0 ) {
                result = Integer.compare(this.getEnd(), that.getEnd());
                if ( result == 0 ) result = this.getIntervalSVType().compareTo(that.getIntervalSVType());
            }
        }
        return result;
    }
}
