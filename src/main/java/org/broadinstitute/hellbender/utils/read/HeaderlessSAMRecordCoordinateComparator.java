package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.io.Serializable;
import java.util.Comparator;

/**
 * A comparator for headerless SAMRecords that exactly matches the ordering of the
 * {@link htsjdk.samtools.SAMRecordCoordinateComparator}
 */
public final class HeaderlessSAMRecordCoordinateComparator implements Comparator<SAMRecord>, Serializable {
    private static final long serialVersionUID = 1L;

    private final SAMFileHeader header;
    private transient SAMRecordCoordinateComparator samComparator;

    public HeaderlessSAMRecordCoordinateComparator( final SAMFileHeader header ) {
        this.header = header;
    }

    @Override
    public int compare( SAMRecord firstRead, SAMRecord secondRead ) {
        // Restore our transient samComparator if it was lost due to serialization
        if ( samComparator == null ) {
            samComparator = new SAMRecordCoordinateComparator();
        }

        // Temporarily set the headers on the reads to our header so that they can be compared
        // using the existing SAMRecordCoordinateComparator
        firstRead.setHeader(header);
        secondRead.setHeader(header);

        final int result = samComparator.compare(firstRead, secondRead);

        // Set the headers on the reads back to null
        firstRead.setHeader(null);
        secondRead.setHeader(null);

        return result;
    }
}
