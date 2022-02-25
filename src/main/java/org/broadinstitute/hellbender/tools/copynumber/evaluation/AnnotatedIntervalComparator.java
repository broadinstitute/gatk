package org.broadinstitute.hellbender.tools.copynumber.evaluation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalCoordinateComparator;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;

import java.util.Comparator;

public class AnnotatedIntervalComparator implements Comparator<AnnotatedInterval> {
    final SAMFileHeader header;
    final IntervalCoordinateComparator intervalCoordinateComparator;

    public AnnotatedIntervalComparator(final SAMFileHeader header) {
        intervalCoordinateComparator = new IntervalCoordinateComparator(header);
        this.header = header;
    }

    @Override
    public int compare(AnnotatedInterval left, AnnotatedInterval right) {
        // This is fine for now
        final Interval leftInterval = new Interval(left);
        final Interval rightInterval = new Interval(right);

        return intervalCoordinateComparator.compare(leftInterval, rightInterval);
    }
}
