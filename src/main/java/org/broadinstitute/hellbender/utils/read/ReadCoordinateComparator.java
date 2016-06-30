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
 * Ordering is almost identical to the {@link htsjdk.samtools.SAMRecordCoordinateComparator},
 * modulo a few subtle differences in tie-breaking rules for reads that share the same
 * position. This comparator will produce an ordering consistent with coordinate ordering
 * in a bam file, including interleaving unmapped reads assigned the positions of their
 * mates with the mapped reads.
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
        int result = compareCoordinates(first, second, header);
        if ( result != 0 ) {
            return result;
        }

        //This is done to mimic SAMRecordCoordinateComparator's behavior
        if (first.isReverseStrand() != second.isReverseStrand()) {
            return first.isReverseStrand()? 1: -1;
        }

        if ( first.getName() != null && second.getName() != null ) {
            result = first.getName().compareTo(second.getName());
            if ( result != 0 ) { return result; }
        }
        result = Integer.compare(ReadUtils.getSAMFlagsForRead(first), ReadUtils.getSAMFlagsForRead(second));
        if ( result != 0 ) { return result; }
        result = Integer.compare(first.getMappingQuality(), second.getMappingQuality());
        if ( result != 0 ) { return result; }
        if (first.isPaired() && second.isPaired()) {
            result = Integer.compare(ReadUtils.getMateReferenceIndex(first, header), ReadUtils.getMateReferenceIndex(second, header));
            if ( result != 0 ) { return result; }
            result = Integer.compare(first.getMateStart(), second.getMateStart());
            if ( result != 0 ) { return result; }
        }
        result = Integer.compare(first.getFragmentLength(), second.getFragmentLength());

        return result;
    }

    public static int compareCoordinates( final GATKRead first, final GATKRead second, final SAMFileHeader header ) {
        final int firstRefIndex = ReadUtils.getAssignedReferenceIndex(first, header);
        final int secondRefIndex = ReadUtils.getAssignedReferenceIndex(second, header);

        if ( firstRefIndex == -1 ) {
            return (secondRefIndex == -1 ? 0 : 1);
        }
        else if ( secondRefIndex == -1 ) {
            return -1;
        }

        final int refIndexDifference = firstRefIndex - secondRefIndex;
        if ( refIndexDifference != 0 ) {
            return refIndexDifference;
        }

        return Integer.compare(first.getAssignedStart(), second.getAssignedStart());
    }
}
