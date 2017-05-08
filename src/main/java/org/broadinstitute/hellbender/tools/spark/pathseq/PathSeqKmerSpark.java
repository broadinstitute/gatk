package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongBloomFilter;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;

import java.util.ArrayList;
import java.util.List;

/**
 * SparkTool to build kmer Bloom filter from a given host reference. The output file is required by PathSeqFilterSpark.
 */
@CommandLineProgramProperties(summary = "Builds kmer library used in the PathSeq subtraction phase.",
        oneLineSummary = "Builds kmer library",
        programGroup = SparkProgramGroup.class)
public final class PathSeqKmerSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "File for kmer library output. Extension will be automatically added if not present ("
            + PSKmerLibUtils.HOPSCOTCH_SET_EXTENSION + " for kmer set or "
            + PSKmerLibUtils.BLOOM_FILTER_EXTENSION + " for Bloom filter)",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String OUTPUT_FILE;

    @Argument(doc = "If given, uses this comma-delimited list of .hss file URIs instead of a reference to construct the kmer library. "
            + "Sets will be merged together. kSize must be set to that of the input libraries, i.e. the full length minus the number of masked bases.",
            fullName = "kmerFiles",
            optional = true)
    public String KMER_URIS = null;

    @Argument(doc = "If set to a non-zero value, creates a Bloom filter with this false positive probability",
            fullName = "bloomFalsePositiveProbability",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double BLOOM_FPP = 0;

    @Argument(doc = "Kmer size, must be odd and less than 32",
            fullName = "kSize",
            minValue = 1,
            maxValue = 31,
            optional = true)
    public int KMER_SIZE = 31;

    @Argument(doc = "Comma-delimited list of (0-based) indices to mask in each kmer, must have an even number of entries"
            + " and be in ascending order. If using --kmerFile, then the indices correspond to those of the unmasked "
            + "bases in the input kmers",
            fullName = "kmerMask",
            optional = true)
    public String KMER_MASK = "";

    @Argument(doc = "Maximum size of the kmer library partitions, in megabytes",
            fullName = "maxPartitionSize",
            minValue = 1,
            maxValue = 2048,
            optional = true)
    public int PARTITION_SIZE_MB = 2048;

    @Argument(doc = "Spacing between sampled kmers. Ignored if --kmerFile is specified.",
            fullName = "kmerSpacing",
            minValue = 1,
            optional = true)
    public int KMER_DIST = 1;

    @Argument(doc = "Downsampling probability. Increasing kmerDistance is more efficient than downsampling and should "
            + "be used instead when possible.",
            fullName = "downsampleProbability",
            minValue = 0.0,
            maxValue = 1.0,
            optional = true)
    public double DOWNSAMPLE_PROBABILITY = 1.0;

    @Argument(doc = "Downsampling seed",
            fullName = "downsampleSeed",
            optional = true)
    public long DOWNSAMPLE_SEED = 1L;

    @Override
    public boolean requiresReference() {
        return false;
    }

    /**
     * Get the list of distinct kmers in the reference, and write them to a file as a HopScotch set or Bloom filter.
     */
    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final String[] kmerFiles = PSUtils.parseCommaDelimitedArgList(KMER_URIS);
        final byte[] kmerMask = PSUtils.parseMask(KMER_MASK, KMER_SIZE);

        final long maxPartitionBytes = PARTITION_SIZE_MB * 1024L * 1024L;
        final PipelineOptions gcsOptions = getAuthenticatedGCSOptions();
        LargeLongHopscotchSet kmerSet;

        if (kmerFiles.length == 0) {
            final SAMSequenceDictionary dict = getReferenceSequenceDictionary();
            final ReferenceMultiSource referenceMultiSource = getReference();
            kmerSet = PSKmerLibUtils.getKmersFromReference(ctx, maxPartitionBytes, KMER_SIZE, KMER_DIST, kmerMask, referenceMultiSource, gcsOptions, dict);
        } else {
            if (kmerFiles.length == 1) {
                logger.debug("Loading kmers from file " + kmerFiles[0] + "...");
                kmerSet = PSKmerLibUtils.readLargeLongHopscotchSet(kmerFiles[0], gcsOptions);
            } else {
                final List<LargeLongHopscotchSet> loadedSets = new ArrayList<>(kmerFiles.length);
                long totalSize = 0;
                for (final String kmerFilePath : kmerFiles) {
                    logger.debug("Loading kmers from file " + kmerFilePath + "...");
                    final LargeLongHopscotchSet set = PSKmerLibUtils.readLargeLongHopscotchSet(kmerFilePath, gcsOptions);
                    loadedSets.add(set);
                    totalSize += set.size();
                }
                kmerSet = new LargeLongHopscotchSet(maxPartitionBytes, totalSize);
                for (final LargeLongHopscotchSet set : loadedSets) {
                    final LongIterator itr = set.iterator();
                    while (itr.hasNext()) {
                        kmerSet.add(itr.next());
                    }
                }
            }
        }

        logger.debug("Loaded " + kmerSet.size() + " kmers in a set with " + kmerSet.getSets().size() + " partitions.");

        if (kmerMask.length > 0) {
            logger.debug("Applying kmer mask...");
            kmerSet = PSKmerLibUtils.maskKmers(kmerSet, KMER_SIZE, kmerMask, maxPartitionBytes);
            logger.debug("Masked set contains " + kmerSet.size() + " kmers.");
        }

        if (DOWNSAMPLE_PROBABILITY < 1) {
            logger.debug("Downsampling kmer set...");
            kmerSet = PSKmerLibUtils.downsampleKmerSet(kmerSet, maxPartitionBytes, DOWNSAMPLE_PROBABILITY, DOWNSAMPLE_SEED);
            logger.debug("Downsampled set contains " + kmerSet.size() + " kmers.");
        }

        if (BLOOM_FPP > 0) {
            logger.debug("Building Bloom filter with false positive probability " + BLOOM_FPP + "...");
            final LargeLongBloomFilter bloomFilter = PSKmerLibUtils.createBloomFilterFromSet(kmerSet, maxPartitionBytes, BLOOM_FPP);
            logger.debug("Built Bloom filter with " + bloomFilter.getSets().size() + " partition(s).");

            final double fpp = LargeLongBloomFilter.estimateFalsePositiveProbability(bloomFilter, (long) (1000 / BLOOM_FPP), 28399432L);
            logger.info("Estimated Bloom filter false positive probability: " + fpp);

            logger.debug("Writing Bloom filter file...");
            PSKmerLibUtils.writeLargeLongBloomFilter(bloomFilter, OUTPUT_FILE, gcsOptions);
        } else {
            logger.debug("Writing kmer set file...");
            PSKmerLibUtils.writeLargeLongHopscotchSet(kmerSet, OUTPUT_FILE, gcsOptions);
        }
    }
}
