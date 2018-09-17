package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
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
public abstract class ReadsKey {

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final GATKRead read) {
        return read.getName();
    }

    /**
     * Makes a hash key for the read.
     */
    public static ReadsKey hashKeyForPassthroughRead(final GATKRead read) {
        return new KeyForFragment(read.getName().hashCode()) ;
    }

    public static ReadsKey getKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, byte library) {
        return new KeyForFragment(longKeyForFragment(strandedUnclippedStart, reverseStrand, referenceIndex, library));
    }

    public static ReadsKey getKeyForPair(final SAMFileHeader header, final GATKRead first, final GATKRead second, final Map<String, Byte> libraryKeyMap) {
        return new KeyForPair(longKeyForFragment(ReadUtils.getStrandedUnclippedStart(first),
                                first.isReverseStrand(),
                                ReadUtils.getReferenceIndex(first, header),
                                libraryKeyMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY))),
                longKeyForPair(ReadUtils.getStrandedUnclippedStart(second),
                        second.isReverseStrand(),
                        ReadUtils.getReferenceIndex(second, header))
                );
    }

    /**
     * Key class for representing relevant duplicate marking identifiers into a single long key for fragment data.
     *
     * Note: This class is intended for internal MarkDuplicatesSpark key purposes, it is only exposed so it can
     *       be accessed by {@link org.broadinstitute.hellbender.engine.spark.GATKRegistrator} for kryo serialization
     */
    public static class KeyForFragment extends ReadsKey {
        final long keyValue;

        KeyForFragment(final long key) {
            this.keyValue = key;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            KeyForFragment that = (KeyForFragment) o;
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
     *
     * Note: This class is intended for internal MarkDuplicatesSpark key purposes, it is only exposed so it can
     *       be accessed by {@link org.broadinstitute.hellbender.engine.spark.GATKRegistrator} for kryo serialization
     */
    public static class KeyForPair extends ReadsKey {
        final long firstReadKeyValue;
        final long secondReadKeyValue;

        KeyForPair(final long firstReadKeyValue, final long secondReadKeyValue) {
            this.firstReadKeyValue = firstReadKeyValue;
            this.secondReadKeyValue = secondReadKeyValue;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            KeyForPair that = (KeyForPair) o;
            return firstReadKeyValue == that.firstReadKeyValue &&
                    secondReadKeyValue == that.secondReadKeyValue;
        }

        @Override
        public int hashCode() {
            return Objects.hash(firstReadKeyValue, secondReadKeyValue);
        }

        @Override
        public String toString() {
            return firstReadKeyValue + " " + secondReadKeyValue;
        }
    }

    // Helper methods for generating summary longs
    private static long longKeyForFragment(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex, byte library) {
        return (((long)strandedUnclippedStart) << 32) |
                        (referenceIndex << 16 & (0xFFFF0000)) | // Note, the bitmasks are being used here because upcasting a negative int to a long in java results in the top bits being filled with 1s, which will ruin the rest of the key. So we mask it for saftey.
                        ((library << 8) & (0x0000FF00)) |
                        (reverseStrand ? 1 : 0);
    }

    private static long longKeyForPair(int strandedUnclippedStart, boolean reverseStrand, int referenceIndex) {
        return (((long)strandedUnclippedStart) << 32) |
                ((referenceIndex << 16) & (0xFFFF0000)) | // Note, the bitmasks are being used here because upcasting a negative int to a long in java results in the top bits being filled with 1s, which will ruin the rest of the key. So we mask it for saftey.
                (reverseStrand ? 1 : 0);
    }
}
