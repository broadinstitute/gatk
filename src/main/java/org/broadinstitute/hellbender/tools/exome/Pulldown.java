package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.stream.Collectors;

/**
 * Simple data structure to pass and write pulldown results.  Should probably replace with a more generic class later.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Pulldown extends AllelicCountCollection {
    private final SAMFileHeader header;

    public Pulldown(final SAMFileHeader header) {
        Utils.nonNull(header, "SAMFileHeader must be supplied.");
        this.header = header;
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file and
     * uses external SAMFile header to construct Pulldown.
     * @param inputFile     file to read from
     * @param header        SAMFile header for IntervalList
     * TODO remove dependency from IntervalList/SamLocusIterator on external header once LocusWalker implemented?
     */
    public Pulldown(final File inputFile, final SAMFileHeader header) {
        super(inputFile);
        Utils.nonNull(header, "SAMFileHeader must be supplied.");
        this.header = header;
    }

    /** Returns the IntervalList of SNP sites.   */
    public IntervalList getIntervals() {
        final IntervalList intervals = new IntervalList(header);
        intervals.addall(getCounts().stream().map(count -> count.getInterval()).collect(Collectors.toList()));
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
