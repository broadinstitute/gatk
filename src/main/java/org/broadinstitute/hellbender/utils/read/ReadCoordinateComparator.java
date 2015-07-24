package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Comparator for sorting Reads by coordinate. Note that a header is required in
 * order to meaningfully compare contigs.
 *
 * Uses the various other fields in a read to break ties for reads that share
 * the same location.
 *
 * Based loosely on the SAMRecordCoordinateComparator from htsjdk.
 */
public final class ReadCoordinateComparator implements Comparator<GATKRead>, Serializable {
    private static final long serialVersionUID = 1L;

    private final SAMFileHeader header;

    public ReadCoordinateComparator( final SAMFileHeader header ) {
        if ( header == null ) {
            throw new IllegalArgumentException("header must be non-null");
        }
        this.header = header;
    }

    @Override
    public int compare( GATKRead first, GATKRead second ) {
        int result = compareCoordinates(first, second);
        if ( result != 0 ) {
            return result;
        }

        result = first.getName().compareTo(second.getName());
        if ( result != 0 ) { return result; }
        result = Integer.compare(ReadUtils.getSAMFlagsForRead(first), ReadUtils.getSAMFlagsForRead(second));
        if ( result != 0 ) { return result; }
        result = Integer.compare(first.getMappingQuality(), second.getMappingQuality());
        if ( result != 0 ) { return result; }
        result = Integer.compare(ReadUtils.getMateReferenceIndex(first, header), ReadUtils.getMateReferenceIndex(second, header));
        if ( result != 0 ) { return result; }
        result = Integer.compare(first.getMateStart(), second.getMateStart());
        if ( result != 0 ) { return result; }
        result = Integer.compare(first.getFragmentLength(), second.getFragmentLength());

        return result;
    }

    private int compareCoordinates( final GATKRead first, final GATKRead second ) {
        final int firstRefIndex = ReadUtils.getReferenceIndex(first, header);
        final int secondRefIndex = ReadUtils.getReferenceIndex(second, header);

        if ( first.isUnmapped() ) {
            return second.isUnmapped() ? 0 : 1;
        }
        else if ( second.isUnmapped() ) {
            return -1;
        }

        final int refIndexComparison = firstRefIndex - secondRefIndex;
        if ( refIndexComparison != 0 ) {
            return refIndexComparison;
        }

        return Integer.compare(first.getStart(), second.getStart());
    }
}
