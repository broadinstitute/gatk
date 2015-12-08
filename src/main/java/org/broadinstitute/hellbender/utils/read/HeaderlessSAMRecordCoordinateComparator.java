package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import java.io.Serializable;
import java.util.Comparator;

/**
 * A comparator for headerless SAMRecords that exactly matches the ordering of the
 * {@link htsjdk.samtools.SAMRecordCoordinateComparator}
 */
public final class HeaderlessSAMRecordCoordinateComparator implements Comparator<SAMRecord>, Serializable {
    private static final long serialVersionUID = 1L;

    private final SAMFileHeader header;

    public HeaderlessSAMRecordCoordinateComparator( final SAMFileHeader header ) {
        this.header = header;
    }

    @Override
    public int compare( final SAMRecord samRecord1, final SAMRecord samRecord2 ) {
        int cmp = compareCoordinates(samRecord1, samRecord2);
        if ( cmp != 0 ) {
            return cmp;
        }

        // Test of negative strand flag is not really necessary, because it is tested
        // via getFlags(), but it is left here because that is the way it is done
        // in {@link htsjdk.samtools.SAMRecordCoordinateComparator}
        if ( samRecord1.getReadNegativeStrandFlag() == samRecord2.getReadNegativeStrandFlag() ) {
            cmp = samRecord1.getReadName().compareTo(samRecord2.getReadName());
            if ( cmp != 0 ) {
                return cmp;
            }
            cmp = Integer.compare(samRecord1.getFlags(), samRecord2.getFlags());
            if ( cmp != 0 ) {
                return cmp;
            }
            cmp = Integer.compare(samRecord1.getMappingQuality(), samRecord2.getMappingQuality());
            if ( cmp != 0 ) {
                return cmp;
            }
            cmp = Integer.compare(header.getSequenceIndex(samRecord1.getMateReferenceName()), header.getSequenceIndex(samRecord2.getMateReferenceName()));
            if ( cmp != 0 ) {
                return cmp;
            }
            cmp = Integer.compare(samRecord1.getMateAlignmentStart(), samRecord2.getMateAlignmentStart());
            if ( cmp != 0 ) {
                return cmp;
            }
            cmp = Integer.compare(samRecord1.getInferredInsertSize(), samRecord2.getInferredInsertSize());
            return cmp;
        }
        else {
            return samRecord1.getReadNegativeStrandFlag() ? 1: -1;
        }
    }

    /**
     * Compare the coordinates of two reads. If a read is paired and unmapped, use its mate mapping
     * as its position.
     *
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
    private int compareCoordinates( final SAMRecord samRecord1, final SAMRecord samRecord2 ) {
        final int refIndex1 = header.getSequenceIndex(samRecord1.getReferenceName());
        final int refIndex2 = header.getSequenceIndex(samRecord2.getReferenceName());

        if ( refIndex1 == -1 ) {
            return refIndex2 == -1 ? 0: 1;
        } else if ( refIndex2 == -1 ) {
            return -1;
        }
        final int cmp = refIndex1 - refIndex2;
        if ( cmp != 0 ) {
            return cmp;
        }
        return samRecord1.getAlignmentStart() - samRecord2.getAlignmentStart();
    }
}
