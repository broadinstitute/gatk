package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.MetagenomicsProgramGroup;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceFileSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.utils.LargeLongHopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.LongBloomFilter;

import java.util.Collection;

/**
 * Produce a set of k-mers from the given host reference. The output file from this tool is required to run the PathSeq pipeline.
 *
 * <p>The tool works by scanning the reference one position at a time. It takes the k-mer (k-base subsequence) starting at each
 * consecutive position and adds it to a set. By default, the set is stored as a
 * <a href="https://en.wikipedia.org/wiki/Hash_table">hash table</a>.</p>
 *
 * <p>Users also have the option to represent the k-mers set using a <a href="https://en.wikipedia.org/wiki/Bloom_filter">Bloom
 * filter</a> by specifying a non-zero value for the --bloom-false-positive-probability parameter. This uses less memory
 * than the default hash set but also can produce false positives. In other words, when
 * asked whether a non-host k-mer exists in the set, it will incorrectly say yes with a probability, p. The user can
 * specify p so that the probability of incorrectly subtracting a non-host read is negligibly
 * small. For p = 0.0001 and read length of 151 bases, the probability of the PathSeq incorrectly subtracting a non-host
 * read is < 1.5%, but the amount of memory used is reduced 4-fold compared to a hash table. For this reason, Bloom
 * filters are generally recommended.</p>
 *
 * <p>Note that the file formats used for storing these k-mer data structures are only readable by the PathSeq tools.</p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>An indexed host reference in FASTA format</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A set of the k-mers in the reference</li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Builds a hash table of every k-mer (k = 31) in the reference. Each k-mer is masked at the 16th position.</h4>
 * <pre>
 * gatk PathSeqBuildKmers  \
 *   --reference host_reference.fasta \
 *   --output host_reference.hss \
 *   --kmer-mask 16 \
 *   --kmer-size 31
 * </pre>
 *
 * <h4>Builds a Bloom filter with false positive probability p < 0.001.</h4>
 * <pre>
 * gatk PathSeqBuildKmers  \
 *   --reference host_reference.fasta \
 *   --output host_reference.hss \
 *   --bloom-false-positive-probability 0.001 \
 *   --kmer-mask 16 \
 *   --kmer-size 31
 * </pre>
 *
 * <h3>Notes</h3>
 *
 * <p>For most references, the Java VM will run out of memory with the default settings. The Java heap size limit should
 * be set at least 20x the size of the reference (less if building a Bloom filter). For example, for a 3 GB reference set
 * the limit to 60 GB by adding --java-options "-Xmx60g" to the command.</p>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Produce a set of k-mers from the given host reference. The output file from this tool is required to run the PathSeq pipeline.",
        oneLineSummary = "Builds set of host reference k-mers",
        programGroup = MetagenomicsProgramGroup.class)
public final class PathSeqBuildKmers extends CommandLineProgram {

    public static final String REFERENCE_LONG_NAME = StandardArgumentDefinitions.REFERENCE_LONG_NAME;
    public static final String REFERENCE_SHORT_NAME = StandardArgumentDefinitions.REFERENCE_SHORT_NAME;
    public static final String BLOOM_FILTER_FALSE_POSITIVE_P_LONG_NAME = "bloom-false-positive-probability";
    public static final String BLOOM_FILTER_FALSE_POSITIVE_P_SHORT_NAME = "P";
    public static final String KMER_SIZE_LONG_NAME = "kmer-size";
    public static final String KMER_SIZE_SHORT_NAME = "SZ";
    public static final String KMER_MASK_LONG_NAME = "kmer-mask";
    public static final String KMER_MASK_SHORT_NAME = "M";
    public static final String KMER_SPACING_LONG_NAME = "kmer-spacing";
    public static final String KMER_SPACING_SHORT_NAME = "SP";

    @Argument(doc = "File for k-mer set output. Extension will be automatically added if not present ("
            + PSKmerUtils.HOPSCOTCH_SET_EXTENSION + " for hash set or "
            + PSKmerUtils.BLOOM_FILTER_EXTENSION + " for Bloom filter)",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String outputFile;

    @Argument(doc = "Reference FASTA file path on local disk",
            fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME)
    public String reference;

    /**
     * <p>Note that the provided argument is used as an upper limit on the probability, and the actual false positive
     * probability may be less.</p>
     */
    @Argument(doc = "If non-zero, creates a Bloom filter with this false positive probability",
            fullName = BLOOM_FILTER_FALSE_POSITIVE_P_LONG_NAME,
            shortName = BLOOM_FILTER_FALSE_POSITIVE_P_SHORT_NAME,
            minValue = 0.0,
            maxValue = 1.0,
            maxRecommendedValue = 0.001,
            optional = true)
    public double bloomFpp = 0;

    /**
     * Reducing the k-mer length will increase the number of host reads subtracted in the
     * filtering phase of the pipeline, but it may also increase the number of non-host (i.e. microbial)
     * reads that are incorrectly subtracted. Note that changing the length of the k-mer does not affect memory usage.
     */
    @Argument(doc = "K-mer size, must be odd and less than 32",
            fullName = KMER_SIZE_LONG_NAME,
            shortName = KMER_SIZE_SHORT_NAME,
            minValue = 1,
            maxValue = 31,
            optional = true)
    public int kmerSize = 31;

    /**
     * K-mer masking allows mismatches to occur at one or more specified positions. Masking the middle base is recommended
     * to enhance host read detection.
     */
    @Argument(doc = "Comma-delimited list of base indices (starting with 0) to mask in each k-mer",
            fullName = KMER_MASK_LONG_NAME,
            shortName = KMER_MASK_SHORT_NAME,
            optional = true)
    public String kmerMaskString = "";

    /**
     * The k-mer set size can be reduced by only storing k-mers starting at every n bases in the reference. By default
     * every k-mer, starting at consecutive bases in the reference, is stored.
     */
    @Argument(doc = "Spacing between successive k-mers",
            fullName = KMER_SPACING_LONG_NAME,
            shortName = KMER_SPACING_SHORT_NAME,
            minValue = 1,
            optional = true)
    public int kmerSpacing = 1;

    /**
     * Get the list of distinct kmers in the reference, and write them to a file as a HopScotch set or Bloom filter.
     */
    @Override
    protected Object doWork() {

        final ReferenceFileSparkSource reference = new ReferenceFileSparkSource(this.reference);

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
