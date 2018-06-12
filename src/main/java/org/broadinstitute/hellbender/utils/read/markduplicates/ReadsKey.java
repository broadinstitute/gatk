package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Map;
import java.util.Objects;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicatesGATK.
 *
 * This class was changed to primarily operate on key hashing instead of generating long string keys as it was discovered
 * that it had performance implications for serialization in MarkDuplicatesSpark. 
 */
public class ReadsKey {

    /**
     * Makes a hash key for the fragment.
     */
    public static int hashKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, String library) {

        int key = library != null ? library.hashCode() : 1;
        key = key * 31 + referenceIndex;
        key = key * 31 + strandedUnclippedStart;
        return key * 31 + (reverseStrand ? 0 : 1);
    }
//
//    /**
//     * Makes a hash key for the paired reads.
//     */
//    public static int hashKeyForPair(final SAMFileHeader header, final GATKRead first, final GATKRead second) {
//        int key = hashKeyForFragment(ReadUtils.getStrandedUnclippedStart(first), first.isReverseStrand(),
//                                     ReadUtils.getReferenceIndex(first, header), ReadUtils.getLibrary(first, header));
//        if (second == null) {
//            return key;
//        }
//
//        key = 31 * key + ReadUtils.getReferenceIndex(second, header);
//        key = 31 * key + ReadUtils.getStrandedUnclippedStart(second);
//        return 31 * key + (second.isReverseStrand() ? 0 : 1);
//    }

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final GATKRead read) {
        return read.getName();
    }

    /**
     * Makes a hash key for the read.
     */
    public static ReadsKey hashKeyForRead(final GATKRead read) {
        return new keyForFragment(read.getName().hashCode()) ;
    }

    public static ReadsKey getKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, byte library) {
        return new keyForFragment(longKeyForFragment(strandedUnclippedStart, reverseStrand, referenceIndex, library));
    }

    public static ReadsKey getKeyForPair(final SAMFileHeader header, final GATKRead first, final GATKRead second, final Map<String, Byte> libraryKeyMap) {
        return new keyForPair(longKeyForFragment(ReadUtils.getStrandedUnclippedStart(first),
                                first.isReverseStrand(),
                                ReadUtils.getReferenceIndex(first, header),
                                libraryKeyMap.get(ReadUtils.getLibrary(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY))),
                longKeyForPair(ReadUtils.getStrandedUnclippedStart(second),
                        second.isReverseStrand(),
                        ReadUtils.getReferenceIndex(second, header))
                );
    }

    /**
     *
     */
    public static class keyForFragment extends ReadsKey {
        final long keyValue;

        keyForFragment(final long key) {
            this.keyValue = key;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            keyForFragment that = (keyForFragment) o;
            return keyValue == that.keyValue;
        }

        @Override
        public int hashCode() {

            return Objects.hash(keyValue);
        }

        @Override
        public String toString() {
            return Long.toString(keyValue);
        }
    }

    /**
     *
     */
    public static class keyForPair extends ReadsKey {
        final long fragmentValue;
        final long keyValue;

        keyForPair(final long fragmentValue, final long pairValue) {
            this.fragmentValue = fragmentValue;
            this.keyValue = pairValue;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            keyForPair that = (keyForPair) o;
            return fragmentValue == that.fragmentValue &&
                    keyValue == that.keyValue;
        }

        @Override
        public int hashCode() {

            return Objects.hash(fragmentValue, keyValue);
        }

        @Override
        public String toString() {
            return fragmentValue + " " + keyValue;
        }
    }


    private static long longKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, byte library) {
        return (((long)strandedUnclippedStart) << 32) |
                        (referenceIndex << 16) |
                        (library << 8) |
                        (reverseStrand ? 1 : 0);
    }
    private static long longKeyForPair(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex) {
        return (((long)strandedUnclippedStart) << 32) |
                (referenceIndex << 16) |
                (reverseStrand ? 1 : 0);
    }
}
