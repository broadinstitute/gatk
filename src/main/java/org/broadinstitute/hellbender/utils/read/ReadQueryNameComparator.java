package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMTag;

import java.io.Serializable;
import java.util.Comparator;

/**
 * compare {@link GATKRead} by queryname
 * duplicates the exact ordering of {@link SAMRecordQueryNameComparator}
 */
public class ReadQueryNameComparator implements Comparator<GATKRead>, Serializable {
    private static final long serialVersionUID = 1L;

    @Override
    public int compare(final GATKRead read1, final GATKRead read2) {
        int cmp = compareReadNames(read1, read2);
        if (cmp != 0) {
            return cmp;
        }

        final boolean r1Paired = read1.isPaired();
        final boolean r2Paired = read2.isPaired();

        if (r1Paired || r2Paired) {
            if (!r1Paired) return 1;
            else if (!r2Paired) return -1;
            else if (read1.isFirstOfPair()  && read2.isSecondOfPair()) return -1;
            else if (read1.isSecondOfPair() && read2.isFirstOfPair()) return 1;
        }

        if (read1.isReverseStrand() != read2.isReverseStrand()) {
            return (read1.isReverseStrand()? 1: -1);
        }
        if (read1.isSecondaryAlignment() != read2.isSecondaryAlignment()) {
            return read2.isSecondaryAlignment()? -1: 1;
        }
        if (read1.isSupplementaryAlignment() != read2.isSupplementaryAlignment()) {
            return read2.isSupplementaryAlignment() ? -1 : 1;
        }
        final Integer hitIndex1 = read1.getAttributeAsInteger(SAMTag.HI.name());
        final Integer hitIndex2 = read2.getAttributeAsInteger(SAMTag.HI.name());
        if (hitIndex1 != null) {
            if (hitIndex2 == null) return 1;
            else {
                cmp = hitIndex1.compareTo(hitIndex2);
                if (cmp != 0) return cmp;
            }
        } else if (hitIndex2 != null) return -1;
        return 0;
    }

    /**
     * compare read names lexicographically without any additional tie breakers
     */
    public int compareReadNames(final GATKRead read1, final GATKRead read2) {
        return read1.getName().compareTo(read2.getName());
    }
}

