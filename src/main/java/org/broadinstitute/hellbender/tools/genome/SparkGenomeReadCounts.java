package org.broadinstitute.hellbender.tools.genome;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SampleCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetWriter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Collects read counts on whole genome sequencing (WGS) alignments using Spark.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * <h3>Examples</h3>
 *
 * <p>
 *     The command encompasses empirically determined parameters for TCGA project data.
 *     You may obtain better results with different parameters.
 * </p>
 *
 * <p>For whole genome sequencing (WGS) data: </p>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" SparkGenomeReadCounts \
 *   --input sample.bam \
 *   --disableReadFilter NotDuplicateReadFilter \
 *   --outputFile sample.readCounts.tsv
 * </pre>
 *
 * <p>
 *     For whole exome sequencing (WGS) data, use CalculateTargetCoverage instead.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Collect read counts on a WGS bam file using Spark.  This creates a set of intervals that span " +
                "the entire genome.  Use the 'binLength' parameter to specify the size of each interval.  By default, any " +
                "contigs X, Y, M, and MT are excluded.\n" +
        "Please see the " + SparkGenomeReadCounts.DROP_NON_AUTOSOMES_LONG_NAME + " option if using this tool on a non-human genome.",
        oneLineSummary = "Collect read counts on a WGS bam file using Spark.",
        programGroup = CopyNumberProgramGroup.class)
@DocumentedFeature
public class SparkGenomeReadCounts extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(SparkGenomeReadCounts.class);

    private static final Set<String> NON_AUTOSOMAL_CONTIGS = new HashSet<>(Arrays.asList("X", "Y", "MT", "M", "x", "y",
            "m", "chrX", "chrY", "chrMT", "chrM", "chrm"));

    public static final String DROP_NON_AUTOSOMES_LONG_NAME = "keepXYMT";
    public static final String DROP_NON_AUTOSOMES_SHORT_NAME = "keepXYMT";

    public static final String BIN_LENGTH_LONG_NAME = "binLength";
    public static final String BIN_LENGTH_SHORT_NAME = "binLength";

    public static final String WRITE_HDF5_LONG_NAME = "writeHdf5";
    public static final String WRITE_HDF5_SHORT_NAME = "hdf5";

    public static final String WRITE_PROPORTIONAL_COVERAGE_LONG_NAME = "writeProportionalCoverage";
    public static final String WRITE_PROPORTIONAL_COVERAGE_SHORT_NAME = "pcov";

    static final String TSV_EXT = ".tsv";
    static final String HDF5_EXT = ".hdf5";
    static final String PROPORTIONAL_COVERAGE_EXT = ".pcov" + TSV_EXT;
    static final String INTERVALS_EXT = ".intervals" + TSV_EXT;

    @Argument(doc = "Keep X, Y, GL*, NC_*, and MT regions.  If this option is not specified, these regions will be dropped, regardless of intervals specified.  Use -L (or -XL) and enable this option for exact specification of intervals.  This option may be removed in the future.",
            fullName = DROP_NON_AUTOSOMES_LONG_NAME,
            shortName = DROP_NON_AUTOSOMES_SHORT_NAME,
            optional = true
    )
    private boolean keepNonAutosomes = false;

    @Argument(doc = "The length of bins (in bp) for each interval specified.  E.g. chr2:1-100 --> chr2:1-50, chr2:51-100 if binLength = 50.",
            fullName = BIN_LENGTH_LONG_NAME,
            shortName = BIN_LENGTH_SHORT_NAME,
            optional = true)
    private int binLength = 1000;

    @Argument(doc = "Whether we should write an additional read-counts HDF5 file with extension " + HDF5_EXT,
            fullName = WRITE_HDF5_LONG_NAME,
            shortName = WRITE_HDF5_SHORT_NAME,
            optional = true)
    private boolean isWritingHdf5 = true;

    @Argument(doc = "Whether we should write an additional proportional-coverage TSV file with extension " + PROPORTIONAL_COVERAGE_EXT,
            fullName = WRITE_PROPORTIONAL_COVERAGE_LONG_NAME,
            shortName = WRITE_PROPORTIONAL_COVERAGE_SHORT_NAME,
            optional = true)
    private boolean isWritingProportionalCoverage = true;

    @Argument(doc = "Output read-counts TSV file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    /**
     * Determine the intervals to consider for read-count collection.  Honors the keepAutosome parameter.
     *
     * <p>Developer's note:  This CLI will always set the attribute, intervals, to a non-null value.</p>
     *
     * @param rawIntervals Specified by the user.  If null, converts to SimpleIntervals specifying the entire
     *                     reference genome.  If keepNonAutosomes is NOT specified, it will prune these intervals (or the
     *                     ones specified by the user), to remove the contigs that are listed in the docs for
     *                     {@link SparkGenomeReadCounts#keepNonAutosomes}
     * @return Never {@code null}  Specified list of intervals.  These will be treated as if user had specified on the
     *  CLI.
     */
    @Override
    protected List<SimpleInterval> editIntervals(final List<SimpleInterval> rawIntervals) {
        List<SimpleInterval> modifiedIntervals = rawIntervals;
        if (rawIntervals == null) {
            modifiedIntervals = IntervalUtils.getAllIntervalsForReference(getReferenceSequenceDictionary());
        }

        if (keepNonAutosomes) {
            return modifiedIntervals;
        }

        // Enforce the elimination of certain contigs when proper option is set.
        logger.info("Dropping non-autosomes, as requested...");
        return modifiedIntervals.stream()
                .filter(s -> !NON_AUTOSOMAL_CONTIGS.contains(s.getContig()))
                .filter(s -> !(s.getContig().startsWith("GL")) && !(s.getContig().startsWith("NC_")))
                .collect(Collectors.toList());
    }

    private void collectReads() {
        if (readArguments.getReadFilesNames().size() != 1) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }

        final SampleCollection sampleCollection = new SampleCollection(getHeaderForReads());
        if (sampleCollection.sampleCount() > 1){
            throw new UserException.BadInput("We do not support bams with more than one sample.");
        }
        final String sampleName = sampleCollection.sampleIds().get(0);
        final String[] commentsForReadCountsFile = {"##fileFormat  = tsv",
                "##commandLine = " + getCommandLine(),
                String.format("##title = Read counts in %d base bins for WGS", binLength)};

        final ReadFilter filter = makeGenomeReadFilter();
        final SAMSequenceDictionary sequenceDictionary = getReferenceSequenceDictionary();

        logger.info("Starting Spark read-count collection...");
        final long readCountCollectionStartTime = System.currentTimeMillis();
        final JavaRDD<GATKRead> rawReads = getReads();
        final JavaRDD<GATKRead> reads = rawReads.filter(filter::test);

        //Note: using a field inside a closure will pull in the whole enclosing object to serialization
        // (which leads to bad performance and can blow up if some objects in the fields are not
        // Serializable - closures always use java Serializable and not Kryo)
        //Solution here is to use a temp variable for binLength because it's just an int.
        final int binLength = this.binLength;
        final JavaRDD<SimpleInterval> readIntervals = reads
                .filter(read -> sequenceDictionary.getSequence(read.getContig()) != null)
                .map(read -> SparkGenomeReadCounts.createKey(read, sequenceDictionary, binLength));
        final Map<SimpleInterval, Long> byKey = readIntervals.countByValue();
        final Set<SimpleInterval> readIntervalKeySet = byKey.keySet();
        final long totalReads = byKey.values().stream().mapToLong(v -> v).sum();
        final long readCountCollectionEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished Spark read-count collection with %d intervals and %d reads. Elapse of %d seconds",
                readIntervalKeySet.size(), totalReads, (readCountCollectionEndTime - readCountCollectionStartTime) / 1000));

        logger.info("Creating full genome bins...");
        final long createGenomeBinsStartTime = System.currentTimeMillis();
        final List<SimpleInterval> fullGenomeBins = createFullGenomeBins(this.binLength);
        List<Target> fullGenomeTargetCollection = createTargetListFromSimpleInterval(fullGenomeBins);
        final String intervalsFile = outputFile.getAbsolutePath().endsWith(TSV_EXT)
                ? FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + INTERVALS_EXT    //if a TSV file, replace extension with .intervals.tsv extension
                : outputFile.getAbsolutePath() + INTERVALS_EXT;                                  //else just append HDF5 extension
        TargetWriter.writeTargetsToFile(new File(intervalsFile), fullGenomeTargetCollection);
        final long createGenomeBinsEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating genome bins. Elapse of %d seconds",
                (createGenomeBinsEndTime - createGenomeBinsStartTime) / 1000));

        logger.info("Creating missing genome bins...");
        final long createMissingGenomeBinsStartTime = System.currentTimeMillis();
        logger.info("Creating missing genome bins: Creating a mutable mapping...");
        final Map<SimpleInterval, Long> byKeyMutable = new HashMap<>();
        byKeyMutable.putAll(byKey);

        logger.info("Creating missing genome bins: Populating mutable mapping with zero counts for empty regions...");
        fullGenomeBins.forEach(b -> byKeyMutable.putIfAbsent(b, 0l));

        final long createMissingGenomeBinsEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating missing genome bins. Elapse of %d seconds",
                (createMissingGenomeBinsEndTime - createMissingGenomeBinsStartTime) / 1000));

        logger.info("Creating final map...");
        final long createFinalMapStartTime = System.currentTimeMillis();
        final SortedMap<SimpleInterval, Long> byKeySorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        byKeySorted.putAll(byKeyMutable);
        final long createFinalMapEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished creating final map. Elapse of %d seconds",
                (createFinalMapEndTime - createFinalMapStartTime) / 1000));

        logger.info("Writing read-counts TSV file...");
        final long writingCovFileStartTime = System.currentTimeMillis();
        ReadCountCollectionUtils.writeReadCountsFromSimpleInterval(new File(outputFile.getAbsolutePath()), sampleName, byKeySorted, commentsForReadCountsFile);
        final long writingCovFileEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished writing read-counts TSV file. Elapse of %d seconds",
                (writingCovFileEndTime - writingCovFileStartTime) / 1000));

        if (isWritingHdf5) {
            logger.info("Writing read-counts HDF5 file...");
            final String hdf5File = outputFile.getAbsolutePath().endsWith(TSV_EXT)
                    ? FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + HDF5_EXT    //if a TSV file, replace extension with HDF5 extension
                    : outputFile.getAbsolutePath() + HDF5_EXT;                                  //else just append HDF5 extension
            final long writingHdf5CovFileStartTime = System.currentTimeMillis();
            writeReadCountsFromSimpleIntervalToHdf5(new File(hdf5File), sampleName,
                    byKeySorted);
            final long writingHdf5CovFileEndTime = System.currentTimeMillis();
            logger.info(String.format("Finished writing read-counts HDF5 file. Elapse of %d seconds",
                    (writingHdf5CovFileEndTime - writingHdf5CovFileStartTime) / 1000));
        }

        if (isWritingProportionalCoverage) {
            logger.info("creating proportional-coverage map...");
            final String[] commentsForProportionalCoverage = {commentsForReadCountsFile[0], commentsForReadCountsFile[1],
                    String.format("##title = Proportional coverage in %d base bins for WGS (total reads: %d)", this.binLength, totalReads)};
            final long pCovFileStartTime = System.currentTimeMillis();
            final SortedMap<SimpleInterval, Double> byKeyProportionalSorted = new TreeMap<>(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
            byKeySorted.forEach((key, value) -> byKeyProportionalSorted.put(key, (double) value / totalReads));
            final long pCovFileEndTime = System.currentTimeMillis();
            logger.info(String.format("Finished creating proportional-coverage map. Elapse of %d seconds",
                    (pCovFileEndTime - pCovFileStartTime) / 1000));
            logger.info("Writing proportional-coverage file ...");
            final String pCovFile = outputFile.getAbsolutePath().endsWith(TSV_EXT)
                    ? FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + PROPORTIONAL_COVERAGE_EXT   //if a TSV file, replace extension with pcov extension
                    : outputFile.getAbsolutePath() + PROPORTIONAL_COVERAGE_EXT;                                 //else just append pcov extension
            final long writingPCovFileStartTime = System.currentTimeMillis();
            ReadCountCollectionUtils.writeReadCountsFromSimpleInterval(new File(pCovFile), sampleName, byKeyProportionalSorted,
                    commentsForProportionalCoverage);
            final long writingPCovFileEndTime = System.currentTimeMillis();
            logger.info(String.format("Finished writing proportional-coverage file. Elapse of %d seconds",
                    -                (writingPCovFileEndTime - writingPCovFileStartTime) / 1000));

        }
    }

    private List<SimpleInterval> createFullGenomeBins(final int binLength){
        return IntervalUtils.cutToShards(getIntervals(), binLength);
    }

    private static SimpleInterval createKey(final GATKRead read, final SAMSequenceDictionary sequenceDictionary, final int binLength) {
        final String contig = read.getContig();
        final int contigLength = sequenceDictionary.getSequence(contig).getSequenceLength();
        final int newStart = Math.min((read.getStart() / binLength) * binLength + 1, contigLength);
        final int newEnd = Math.min(newStart + binLength - 1, contigLength);

        return new SimpleInterval(contig, newStart, newEnd);
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        collectReads();
    }

    /**
     * Filter reads based on marked parameters
     * @return never {@code null}
     */
    private ReadFilter makeGenomeReadFilter() {
        return new WellformedReadFilter(getHeaderForReads())
                .and(ReadFilterLibrary.MAPPED)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT)
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
    }

    private List<Target> createTargetListFromSimpleInterval(final List<SimpleInterval> intervalList){
        return intervalList.stream().map(Target::new).collect(Collectors.toList());
    }

    /**
     *  Write a simple interval to number mapping as a HDF5 file.
     *
     *  When this file is read, all numbers will have been converted to double.  See {@link ReadCountCollectionUtils::parseHdf5AsDouble}
     *
     * @param outFile output file.  Cannot already exist.  Not {@code null}
     * @param sampleName Name for the values per interval.  Not {@code null}
     * @param byKeySorted Mapping of interval to number.  Not {@code null}
     * @param <N> Class that extends Number.  Stored as a double.
     */
    private static <N extends Number> void writeReadCountsFromSimpleIntervalToHdf5(final File outFile, final String sampleName,
                                                                                   final SortedMap<SimpleInterval, N> byKeySorted) {
        Utils.nonNull(outFile, "Output file cannot be null.");
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        Utils.nonNull(byKeySorted, "Intervals cannot be null.");

        final List<SimpleCount> simpleCounts = new ArrayList<>(byKeySorted.size());
        byKeySorted.forEach((key, value) -> simpleCounts.add(new SimpleCount(key, value.intValue())));
        final SimpleCountCollection scc = new SimpleCountCollection(new SimpleSampleMetadata(sampleName), simpleCounts);
        scc.writeHDF5(outFile);
    }
}

