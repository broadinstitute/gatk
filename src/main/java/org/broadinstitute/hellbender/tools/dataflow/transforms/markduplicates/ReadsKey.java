package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Encodes a unique key for read, read pairs and fragments. Used to identify duplicates for MarkDuplicates.
 */
final class ReadsKey {

    /**
     * Makes a unique key for the fragment.
     */
    public static String keyForFragment(final SAMFileHeader header, final GATKRead read) {
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
        final String key = keyForFragment(header, first);
        if (second == null) {
            return key;
        }

        return String.format(
                "%s|%d|%d|%s",
                key,
                ReadUtils.getReferenceIndex(second, header),
                ReadUtils.getStrandedUnclippedStart(second),
                second.isReverseStrand() ? "r" : "f");
    }

    /**
     * Makes a unique key for the pair.
     */
    public static String keyForPair(final SAMFileHeader header, final GATKRead read) {
        return String.format(
                "%s|%s",
                read.getReadGroup(),
                read.getName());
    }
}

