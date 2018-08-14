package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Common functions for PathSeq
 */
public final class PSUtils {

    private static final Logger logger = LogManager.getLogger(PSUtils.class);

    public static JavaRDD<GATKRead> primaryReads(final JavaRDD<GATKRead> reads) {
        return reads.filter(read -> !(read.isSecondaryAlignment() || read.isSupplementaryAlignment()));
    }

    public static String[] parseCommaDelimitedArgList(final String arg) {
        if (arg == null || arg.isEmpty()) {
            return new String[0];
        }
        return arg.split(",");
    }

    /**
     * Parses command-line option for specifying kmer spacing masks
     */
    public static byte[] parseMask(final String maskArg, final int kSize) {

        final String[] kmerMaskString = parseCommaDelimitedArgList(maskArg);
        final byte[] kmerMask = new byte[kmerMaskString.length];
        for (int i = 0; i < kmerMaskString.length; i++) {
            kmerMask[i] = (byte) Integer.parseInt(kmerMaskString[i]);
            Utils.validateArg(kmerMask[i] >= 0 && kmerMask[i] < kSize, "Invalid kmer mask index: " + kmerMaskString[i]);
        }
        return kmerMask;
    }

    /**
     * Prints warning message followed by a list of relevant items
     */
    public static void logItemizedWarning(final Logger logger, final Collection<?> items, final String warning) {
        if (!items.isEmpty()) {
            final String str =  items.stream().map(String::valueOf).collect(Collectors.joining(", "));
            logger.warn(warning + " : " + str);
        }
    }

    /**
     * Same as GATKSparkTool's getRecommendedNumReducers(), but can specify input BAM path (for when --input is not used)
     */
    public static int pathseqGetRecommendedNumReducers(final String inputPath, final int numReducers,
                                                       final int targetPartitionSize) {
        if (numReducers != 0) {
            return numReducers;
        }
        return 1 + (int) (BucketUtils.dirSize(inputPath) / targetPartitionSize);
    }

    /**
     * Returns a deep copy of the input header with an empty sequence dictionary, and logs warnings if the input may
     * be aligned but --is-host-aligned was not set to true (or vice versa).
     */
    public static SAMFileHeader checkAndClearHeaderSequences(final SAMFileHeader inputHeader, final PSFilterArgumentCollection filterArgs, final Logger logger) {

        Utils.nonNull(inputHeader, "Cannot check and clear null input header");
        Utils.nonNull(filterArgs, "Cannot check header against null filter arguments");
        Utils.nonNull(logger, "Cannot check header using null logger");

        //Deep copy of header, otherwise aligned reads will be filtered out by WellformedReadFilter because the sequence dictionary is cleared
        final SAMFileHeader header = inputHeader.clone();

        if (filterArgs.alignedInput && (header.getSequenceDictionary() == null || header.getSequenceDictionary().isEmpty())) {
            logger.warn("--is-host-aligned is true but the BAM header contains no sequences");
        }
        if (!filterArgs.alignedInput && header.getSequenceDictionary() != null && !header.getSequenceDictionary().isEmpty()) {
            logger.warn("--is-host-aligned is false but there are one or more sequences in the BAM header");
        }

        //Clear header sequences
        header.setSequenceDictionary(new SAMSequenceDictionary());

        return header;
    }

    public static int getMatchesLessDeletions(final Cigar cigar, final int nmTagValue) {
        Utils.nonNull(cigar, "Cannot get match score for null cigar");
        Utils.validateArg(nmTagValue >= 0, "NM tag value cannot be negative");
        int numMatchesOrMismatches = 0;
        int numDeletions = 0;
        int numInsertions = 0;
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        for (final CigarElement e : cigarElements) {
            if (e.getOperator().isAlignment()) {
                numMatchesOrMismatches += e.getLength();
            } else if (e.getOperator().equals(CigarOperator.DELETION)) {
                numDeletions += e.getLength();
            } else if (e.getOperator().equals(CigarOperator.INSERTION)) {
                numInsertions += e.getLength();
            }
        }
        final int numSubstitutions = nmTagValue - numDeletions - numInsertions;
        final int numMatches = numMatchesOrMismatches - numSubstitutions;
        if (numSubstitutions < 0) {
            logger.warn("Invalid arguments passed to getMatchesLessDeletions(): NM tag value was less than the number of insertions and deletions combined. Returning 0.");
            return 0;
        }
        if (numMatches < 0) {
            logger.warn("Invalid arguments passed to getMatchesLessDeletions(): Combined number of matches and mismatches was less than the number of substitutions. Returning 0.");
            return 0;
        }
        return numMatches - numDeletions;
    }
}
