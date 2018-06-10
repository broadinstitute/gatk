package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicatesGATK.
 *
 * This class was changed to primarily operate on key hashing instead of generating long string keys as it was discovered
 * that it had performance implications for serialization in MarkDuplicatesSpark. 
 */
public final class ReadsKey {

    /**
     * Makes a hash key for the fragment.
     */
    public static String hashKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, String library) {

//        int key = library != null ? library.hashCode() : 1;
//        key = key * 31 + referenceIndex;
//        key = key * 31 + strandedUnclippedStart;
//        return key * 31 + (reverseStrand ? 0 : 1);
//
        return new StringBuilder(library!=null?library:"").append('|').append(referenceIndex).append('|').append(strandedUnclippedStart).append(reverseStrand ? 'r' : 'f').toString();
    }

    /**
     * Makes a hash key for the paired reads.
     */
    public static String hashKeyForPair(final SAMFileHeader header, final GATKRead first, final GATKRead second) {
        String key = hashKeyForFragment(ReadUtils.getStrandedUnclippedStart(first), first.isReverseStrand(),
                                     ReadUtils.getReferenceIndex(first, header), ReadUtils.getLibrary(first, header));
        if (second == null) {
            return key;
        }
        return new StringBuffer(key).append(ReadUtils.getReferenceIndex(second, header)).append('|').append(ReadUtils.getStrandedUnclippedStart(second)).append(second.isReverseStrand() ? 'r' : 'f').toString();

//        key = 31 * key + ReadUtils.getReferenceIndex(second, header);
//        key = 31 * key + ReadUtils.getStrandedUnclippedStart(second);
//        return 31 * key + (second.isReverseStrand() ? 0 : 1);
    }

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final GATKRead read) {
        return read.getName();
    }

    /**
     * Makes a hash key for the read.
     */
    public static int hashKeyForRead(final GATKRead read) {
        return read.getName().hashCode();
    }
}
