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

    public HeaderlessSAMRecordCoordinateComparator( final SAMFileHeader header ) {
        this.header = header;
    }

    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = fileOrderCompare(samRecord1, samRecord2);
        if (cmp != 0) {
            return cmp;
        }
        // Test of negative strand flag is not really necessary, because it is tested
        // with cmp if getFlags, but it is left here because that is the way it was done
        // in the past.
        if (samRecord1.getReadNegativeStrandFlag() == samRecord2.getReadNegativeStrandFlag()) {
            cmp = samRecord1.getReadName().compareTo(samRecord2.getReadName());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getFlags(), samRecord2.getFlags());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getMappingQuality(), samRecord2.getMappingQuality());
            if (cmp != 0) return cmp;
            cmp = compareInts(header.getSequenceIndex(samRecord1.getMateReferenceName()), header.getSequenceIndex(samRecord2.getMateReferenceName()));
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getMateAlignmentStart(), samRecord2.getMateAlignmentStart());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getInferredInsertSize(), samRecord2.getInferredInsertSize());
            return cmp;

        }
        else return (samRecord1.getReadNegativeStrandFlag()? 1: -1);
    }

    private int compareInts(int i1, int i2) {
        if (i1 < i2) return -1;
        else if (i1 > i2) return 1;
        else return 0;
    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.  If read is paired and unmapped, use the mate mapping to sort.
     * Records being compared must have non-null SAMFileHeaders.
     *
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        final int refIndex1 = header.getSequenceIndex(samRecord1.getReferenceName());
        final int refIndex2 = header.getSequenceIndex(samRecord2.getReferenceName());

        if (refIndex1 == -1) {
            return (refIndex2 == -1? 0: 1);
        } else if (refIndex2 == -1) {
            return -1;
        }
        final int cmp = refIndex1 - refIndex2;
        if (cmp != 0) {
            return cmp;
        }
        return samRecord1.getAlignmentStart() - samRecord2.getAlignmentStart();
    }
}
