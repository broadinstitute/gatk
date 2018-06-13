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
     * Key class for representing relevant duplicate marking identifiers into a single long key for fragment data.
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
     * Key class for representing relevant duplicate marking identifiers into a two long key values for pair data data.
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

    // Helper methods for generating summary longs
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
