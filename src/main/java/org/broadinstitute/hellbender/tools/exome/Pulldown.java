package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.stream.Collectors;

/**
 * Simple data structure to pass pulldown results.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Pulldown extends AllelicCountCollection {
    private final SAMFileHeader header;

    public Pulldown(final SAMFileHeader header) {
        this.header = Utils.nonNull(header, "SAMFileHeader must be supplied.");;
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file and
     * uses external SAMFile header to construct Pulldown.
     * @param inputFile     file to read from
     * @param header        SAMFile header for IntervalList
     */
    public Pulldown(final File inputFile, final SAMFileHeader header) {
        super(inputFile);
        this.header = Utils.nonNull(header, "SAMFileHeader must be supplied.");
    }

    /** Returns a new instance of an IntervalList, constructed from the intervals of the internally held
     * AllelicCounts.  This IntervalList is modifiable and does not change with the state of the Pulldown.   */
    public IntervalList getIntervals() {
        final IntervalList intervals = new IntervalList(header);
        intervals.addall(getCounts().stream().map(AllelicCount::getInterval)
                .map(si -> new Interval(si.getContig(), si.getStart(), si.getEnd())).collect(Collectors.toList()));
        return intervals;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof Pulldown)) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        final Pulldown pulldown = (Pulldown) o;
        return header.equals(pulldown.header);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + header.hashCode();
        return result;
    }
}
