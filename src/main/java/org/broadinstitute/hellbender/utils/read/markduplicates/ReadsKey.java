package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicatesGATK.
 */
public final class ReadsKey {

    public static final String FRAGMENT_PREFIX = "f";

    public static final String DELIMETER = "|";

    private static final String PAIRED_ENDS_PREFIX = "p";

    /**
     * Makes a unique key for the fragment.
     */
    public static int keyForFragment(final SAMFileHeader header, final GATKRead read, final int strandedUnclippedStart) {
        //return FRAGMENT_PREFIX + DELIMETER + subkeyForFragment(header, read);
        return subkeyForFragment(header, read, strandedUnclippedStart);
    }

    /**
     * Makes a unique key for the fragment (excluding the prefix).
     */
    private static int subkeyForFragment(final SAMFileHeader header, final GATKRead read, int strandedUnclippedStart) {
        final String library = ReadUtils.getLibrary(read, header);
//
//        return new StringBuilder().append(library != null ? library : "-").append(DELIMETER)
//                .append(ReadUtils.getReferenceIndex(read, header)).append(DELIMETER)
//                .append(ReadUtils.getStrandedUnclippedStart(read)).append(DELIMETER)
//                .append(read.isReverseStrand() ? "r" : "f").toString();

        int key = library.hashCode();
        key = key * 31 + ReadUtils.getReferenceIndex(read, header);
        key = key * 31 + ReadUtils.getStrandedUnclippedStart(read);
        return key * 31 + (read.isReverseStrand() ? 0 : 1);

    }

    /**
     * Makes a unique key for the paired reads.
     */
    public static int keyForPairedEnds(final SAMFileHeader header, final GATKRead first, final GATKRead second, int strandedUnclippedStart) {
        int key = subkeyForFragment(header, first, strandedUnclippedStart);
        if (second == null) {
            return key;
        }

//        return new StringBuilder().append(PAIRED_ENDS_PREFIX).append(key).append(DELIMETER)
//                .append(ReadUtils.getReferenceIndex(second, header)).append(DELIMETER)
//                .append(ReadUtils.getStrandedUnclippedStart(second)).append(DELIMETER)
//                .append(second.isReverseStrand() ? "r" : "f").toString();

        key = 31 * key + ReadUtils.getReferenceIndex(second, header);
        key = 31 * key + ReadUtils.getStrandedUnclippedStart(second);
        return 31 * key + (second.isReverseStrand() ? 0 : 1);

    }

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final SAMFileHeader header, final GATKRead read) {
        return read.getName();
        //return new StringBuilder().append(read.getReadGroup()).append(DELIMETER).append(read.getName()).toString();
    }

    /**
     * Returns true if the key is a fragment key.
     */
    public static boolean isFragment(String key) {
        return key.charAt(0)=='f';
    }
}
