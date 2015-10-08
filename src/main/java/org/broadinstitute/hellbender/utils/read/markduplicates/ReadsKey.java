package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicates.
 */
public final class ReadsKey {

    public static final String FRAGMENT_PREFIX = "f|";

    private static final String PAIRED_ENDS_PREFIX = "p|";

    /**
     * Makes a unique key for the fragment.
     */
    public static String keyForFragment(final SAMFileHeader header, final GATKRead read) {
        return FRAGMENT_PREFIX + subkeyForFragment(header, read);
    }

    /**
     * Makes a unique key for the fragment (excluding the prefix).
     */
    private static String subkeyForFragment(final SAMFileHeader header, final GATKRead read) {
        final String library = ReadUtils.getLibrary(read, header);

        return String.format(
                "%s|%d|%d|%s",
                library != null ? library : "-",
                ReadUtils.getReferenceIndex(read, header),
                ReadUtils.getStrandedUnclippedStart(read),
                read.isReverseStrand() ? "r" : "f");
    }

    /**
     * Makes a unique key for the paired reads.
     */
    public static String keyForPairedEnds(final SAMFileHeader header, final GATKRead first, final GATKRead second) {
        final String key = subkeyForFragment(header, first);
        if (second == null) {
            return PAIRED_ENDS_PREFIX + key;
        }

        return String.format(
                PAIRED_ENDS_PREFIX + "%s|%d|%d|%s",
                key,
                ReadUtils.getReferenceIndex(second, header),
                ReadUtils.getStrandedUnclippedStart(second),
                second.isReverseStrand() ? "r" : "f");
    }

    /**
     * Makes a unique key for the read.
     */
    public static String keyForRead(final SAMFileHeader header, final GATKRead read) {
        return String.format(
                "%s|%s",
                read.getReadGroup(),
                read.getName());
    }

    /**
     * Returns true if the key is a fragment key.
     */
    public static boolean isFragment(String key) {
        return key.startsWith(FRAGMENT_PREFIX);
    }
}
