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
    public static String keyForFragment(final SAMFileHeader header, final GATKRead read) {
        return FRAGMENT_PREFIX + DELIMETER + subkeyForFragment(header, read);
    }

    /**
     * Makes a unique key for the fragment (excluding the prefix).
     */
    private static String subkeyForFragment(final SAMFileHeader header, final GATKRead read) {
        final String library = ReadUtils.getLibrary(read, header);

        return new StringBuilder().append(library != null ? library : "-").append(DELIMETER)
                .append(ReadUtils.getReferenceIndex(read, header)).append(DELIMETER)
                .append(ReadUtils.getStrandedUnclippedStart(read)).append(DELIMETER)
                .append(read.isReverseStrand() ? "r" : "f").toString();
    }

    /**
     * Makes a unique key for the paired reads.
     */
    public static String keyForPairedEnds(final SAMFileHeader header, final GATKRead first, final GATKRead second) {
        final String key = subkeyForFragment(header, first);
        if (second == null) {
            return PAIRED_ENDS_PREFIX + DELIMETER + key;
        }

        return new StringBuilder().append(PAIRED_ENDS_PREFIX).append(key).append(DELIMETER)
                .append(ReadUtils.getReferenceIndex(second, header)).append(DELIMETER)
                .append(ReadUtils.getStrandedUnclippedStart(second)).append(DELIMETER)
                .append(second.isReverseStrand() ? "r" : "f").toString();
    }

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final SAMFileHeader header, final GATKRead read) {
        return new StringBuilder().append(read.getReadGroup()).append(DELIMETER).append(read.getName()).toString();
    }

    /**
     * Returns true if the key is a fragment key.
     */
    public static boolean isFragment(String key) {
        return key.charAt(0)=='f';
    }
}
