package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicatesGATK.
 */
public final class ReadsKey {

    /**
     * Makes a unique key for the fragment.
     */
    public static int hashKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, String library) {

        int key = library != null ? library.hashCode() : 1;
        key = key * 31 + referenceIndex;
        key = key * 31 + strandedUnclippedStart;
        return key * 31 + (reverseStrand ? 0 : 1);
    }

    /**
     * Makes a unique key for the paired reads.
     */
    public static int hashKeyForPair(final SAMFileHeader header, final GATKRead first, final GATKRead second) {
        int key = hashKeyForFragment(ReadUtils.getStrandedUnclippedStart(first), first.isReverseStrand(),
                                     ReadUtils.getReferenceIndex(first, header), ReadUtils.getLibrary(first, header));
        if (second == null) {
            return key;
        }

        key = 31 * key + ReadUtils.getReferenceIndex(second, header);
        key = 31 * key + ReadUtils.getStrandedUnclippedStart(second);
        return 31 * key + (second.isReverseStrand() ? 0 : 1);
    }

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final GATKRead read) {
        return read.getName();
    }
}
