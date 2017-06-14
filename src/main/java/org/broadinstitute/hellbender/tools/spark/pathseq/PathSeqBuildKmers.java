package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongBloomFilter;

import java.util.Collection;

/**
 * SparkTool to build kmer hash set or Bloom filter from a given host reference. The output file is required by PathSeqFilterSpark.
 */
@CommandLineProgramProperties(summary = "Builds host kmer set used in the PathSeq subtraction phase.",
        oneLineSummary = "Builds PathSeq kmer library",
        programGroup = SparkProgramGroup.class)
public final class PathSeqBuildKmers extends CommandLineProgram {

    @Argument(doc = "File for kmer library output. Extension will be automatically added if not present ("
            + PSKmerUtils.HOPSCOTCH_SET_EXTENSION + " for hash set or "
            + PSKmerUtils.BLOOM_FILTER_EXTENSION + " for Bloom filter)",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String outputFile;

    @Argument(doc = "Reference file path on local disk",
            fullName = "referencePath")
    public String referencePath;

    @Argument(doc = "If non-zero, creates a Bloom filter with this false positive probability",
            fullName = "bloomFalsePositiveProbability",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double bloomFpp = 0;

    @Argument(doc = "Kmer size, must be odd and less than 32",
            fullName = "kSize",
            minValue = 1,
            maxValue = 31,
            optional = true)
    public int kmerSize = 31;

    @Argument(doc = "Comma-delimited list of base indices (starting with 0) to mask in each kmer",
            fullName = "kmerMask",
            optional = true)
    public String kmerMaskString = "";

    @Argument(doc = "Spacing between successive kmers",
            fullName = "kmerSpacing",
            minValue = 1,
            optional = true)
    public int kmerSpacing = 1;

    /**
     * Get the list of distinct kmers in the reference, and write them to a file as a HopScotch set or Bloom filter.
     */
    @Override
    protected Object doWork() {

        final ReferenceFileSource reference = new ReferenceFileSource(referencePath);

        final byte[] maskBytes = PSUtils.parseMask(kmerMaskString, kmerSize);
        final SVKmerShort kmerMask = SVKmerShort.getMask(maskBytes, kmerSize);

        logger.info("Loading reference kmers...");
        final Collection<long[]> maskedKmerCollection = PSKmerUtils.getMaskedKmersFromLocalReference(reference, kmerSize, kmerSpacing, kmerMask);
        final long numLongs = PSKmerUtils.longArrayCollectionSize(maskedKmerCollection);
        if (bloomFpp > 0) {
            logger.info("Building Bloom filter with false positive probability " + bloomFpp + "...");
            final LongBloomFilter bloomFilter = PSKmerUtils.longArrayCollectionToBloomFilter(maskedKmerCollection, numLongs, bloomFpp);
            final PSKmerBloomFilter kmerBloomFilter = new PSKmerBloomFilter(bloomFilter, kmerSize, kmerMask, numLongs);
            logger.info("Theoretical Bloom filter false positive probability: " + kmerBloomFilter.getFalsePositiveProbability());
            PSKmerUtils.writeKmerBloomFilter(outputFile, kmerBloomFilter);
        } else {
            logger.info("Building kmer hash set...");
            final LargeLongHopscotchSet kmerHopscotchSet = PSKmerUtils.longArrayCollectionToSet(maskedKmerCollection, numLongs);
            final PSKmerSet kmerSet = new PSKmerSet(kmerHopscotchSet, kmerSize, kmerMask);
            PSKmerUtils.writeKmerSet(outputFile, kmerSet);
        }
        return null;
    }
}
