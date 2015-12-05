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
 * Based loosely on the {@link htsjdk.samtools.SAMRecordCoordinateComparator}, but does not
 * match its order for all reads (see explanation below). Use {@link HeaderlessSAMRecordCoordinateComparator}
 * if an exact match against the ordering of {@link htsjdk.samtools.SAMRecordCoordinateComparator}
 * is required (eg., you are writing a bam).
 *
 * This comparator, when given two GATKReads that are backed by SAMRecords will produce the same exact order as {@link htsjdk.samtools.SAMRecordCoordinateComparator}.
 * EXCEPT for read(s) that are in a state that is disallowed by the GATKRead interface, eg: a read that is
 * unmapped but has a reference index that is not {@link htsjdk.samtools.SAMRecord#NO_ALIGNMENT_REFERENCE_INDEX}.
 * For such reads, the order will not be the same as the {@link htsjdk.samtools.SAMRecordCoordinateComparator}.
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

        //This is done to mimic SAMRecordCoordinateComparator's behavior
        if (first.isReverseStrand() != second.isReverseStrand()){
            return first.isReverseStrand()? 1: -1;
        }

        result = first.getName().compareTo(second.getName());
        if ( result != 0 ) { return result; }
        result = Integer.compare(ReadUtils.getSAMFlagsForRead(first), ReadUtils.getSAMFlagsForRead(second));
        if ( result != 0 ) { return result; }
        result = Integer.compare(first.getMappingQuality(), second.getMappingQuality());
        if ( result != 0 ) { return result; }
        if (first.isPaired()) {
            result = Integer.compare(ReadUtils.getMateReferenceIndex(first, header), ReadUtils.getMateReferenceIndex(second, header));
            if ( result != 0 ) { return result; }
            result = Integer.compare(first.getMateStart(), second.getMateStart());
            if ( result != 0 ) { return result; }
        }
        result = Integer.compare(first.getFragmentLength(), second.getFragmentLength());

        return result;
    }

    private int compareCoordinates( final GATKRead first, final GATKRead second ) {
        if ( first.isUnmapped() ) {
            return second.isUnmapped() ? 0 : 1;
        }
        if ( second.isUnmapped() ) {
            return -1;
        }

        final int firstRefIndex = ReadUtils.getReferenceIndex(first, header);
        final int secondRefIndex = ReadUtils.getReferenceIndex(second, header);

        final int refIndexComparison = firstRefIndex - secondRefIndex;
        if ( refIndexComparison != 0 ) {
            return refIndexComparison;
        }

        return Integer.compare(first.getStart(), second.getStart());
    }
}
