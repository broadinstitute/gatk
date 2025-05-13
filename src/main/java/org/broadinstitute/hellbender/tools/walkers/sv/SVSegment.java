package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Comparator;
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

    /**
     * Get comparator for {@link SVSegment} objects. SVsegment is not implemented as a {@link Comparable} due to the
     * need for a common sequence dictionary.
     * @param dictionary SAMSequenceDictionary pertaining to both objects
     * @param <T> object class
     * @return comparator
     */
    public static <T extends SVSegment> Comparator<T> getSVSegmentComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareSVSegments(o1, o2, dictionary);
    }

    /**
     * Compares two SVSegments using sequence dictionary based on contig, then start position, then end position, then SV type
     * @param first SVSegment
     * @param second SVSegment
     * @param dictionary SAMSequenceDictionary to dictate comparison order
     * @return 0 if first and second are equal, a negative value if first < second, or a positive value if first > second
     */
    public static int compareSVSegments(final SVSegment first, final SVSegment second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        final Comparator<Locatable> locatableComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final int comparePositions = locatableComparator.compare(first.getInterval(), second.getInterval());
        if (comparePositions != 0) return comparePositions;
        return first.getIntervalSVType().compareTo(second.getIntervalSVType());
    }
}
