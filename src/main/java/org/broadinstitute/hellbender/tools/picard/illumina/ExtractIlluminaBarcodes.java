package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.IlluminaProgramGroup;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ClusterData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProvider;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProviderFactory;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadDescriptor;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.broadinstitute.hellbender.utils.illumina.IlluminaUtil;
import org.broadinstitute.hellbender.utils.text.parsers.TabbedTextFileWithHeaderParser;

import java.io.BufferedWriter;
import java.io.File;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static htsjdk.samtools.util.IOUtil.assertDirectoryIsWritable;
import static htsjdk.samtools.util.IOUtil.assertFileIsWritable;
import static htsjdk.samtools.util.IOUtil.openFileForBufferedWriting;
import static htsjdk.samtools.util.Log.getInstance;
import static htsjdk.samtools.util.SequenceUtil.basesEqual;
import static htsjdk.samtools.util.SequenceUtil.isNoCall;
import static htsjdk.samtools.util.StringUtil.bytesToString;
import static htsjdk.samtools.util.StringUtil.join;
import static htsjdk.samtools.util.StringUtil.repeatCharNTimes;
import static htsjdk.samtools.util.StringUtil.stringToBytes;
import static java.lang.Math.min;
import static java.lang.Runtime.getRuntime;
import static java.lang.String.format;
import static java.lang.String.valueOf;
import static java.text.NumberFormat.getNumberInstance;
import static java.util.Arrays.asList;
import static java.util.concurrent.Executors.newFixedThreadPool;
import static java.util.concurrent.TimeUnit.HOURS;
import static java.util.concurrent.TimeUnit.SECONDS;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.LANE_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME;
import static org.broadinstitute.hellbender.tools.picard.illumina.ExtractIlluminaBarcodes.BarcodeMetric.copy;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.BaseCalls;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.PF;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.QualityScores;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure.PARAMETER_DOC;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadType.Barcode;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.BARCODE_DELIMITER;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.barcodeSeqsToString;

/**
 * Determine the barcode for each read in an Illumina lane.
 * For each tile, a file is written to the basecalls directory of the form s_<lane>_<tile>_barcode.txt.
 * An output file contains a line for each read in the tile, aligned with the regular basecall output
 * The output file contains the following tab-separated columns:
 * - read subsequence at barcode position
 * - Y or N indicating if there was a barcode match
 * - matched barcode sequence (empty if read did not match one of the barcodes).  If there is no match
 * but we're close to the threshold of calling it a match we output the barcode that would have been
 * matched but in lower case
 *
 * @author jburke@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Determine the barcode for each read in an Illumina lane.\n" +
                "For each tile, a file is written to the basecalls directory of the form s_<lane>_<tile>_barcode.txt. " +
                "An output file contains a line for each read in the tile, aligned with the regular basecall output. \n" +
                "The output file contains the following tab-separated columns: \n" +
                "    * read subsequence at barcode position\n" +
                "    * Y or N indicating if there was a barcode match\n" +
                "    * matched barcode sequence\n" +
                "Note that the order of specification of barcodes can cause arbitrary differences in output for poorly matching barcodes.\n\n",
        usageShort = "Tool to determine the barcode for each read in an Illumina lane",
        programGroup = IlluminaProgramGroup.class
)
public class ExtractIlluminaBarcodes extends PicardCommandLineProgram {

    @Argument(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "Where to write _barcode.txt files.  By default, these are written to BASECALLS_DIR.", optional = true)
    public File OUTPUT_DIR;

    @Argument(doc = "Lane number. ", shortName = LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(doc = PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(doc = "Barcode sequence.  These must be unique, and all the same length.  This cannot be used with reads that " +
            "have more than one barcode; use BARCODE_FILE in that case. ", mutex = {"BARCODE_FILE"})
    public List<String> BARCODE = new ArrayList<>();

    @Argument(doc = "Tab-delimited file of barcode sequences, barcode name and, optionally, library name.  " +
            "Barcodes must be unique and all the same length.  Column headers must be 'barcode_sequence_1', " +
            "'barcode_sequence_2' (optional), 'barcode_name', and 'library_name'.", mutex = {"BARCODE"})
    public File BARCODE_FILE;

    @Argument(doc = "Per-barcode and per-lane metrics written to this file.", shortName = METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Argument(doc = "Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Argument(doc = "Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 1;

    @Argument(doc = "Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Argument(shortName = "Q", doc = "Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even in the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Argument(doc = "The minimum quality (after transforming 0s to 1s) expected from reads.  If qualities are lower than this value, an error is thrown." +
            "The default of 2 is what the Illumina's spec describes as the minimum, but in practice the value has been observed lower.")
    public int MINIMUM_QUALITY = ILLUMINA_ALLEGED_MINIMUM_QUALITY;

    @Argument(shortName = "GZIP", doc = "Compress output s_l_t_barcode.txt files using gzip and append a .gz extension to the file names.")
    public boolean COMPRESS_OUTPUTS = false;

    @Argument(doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.")
    public int NUM_PROCESSORS = 1;

    private static final Log LOG = getInstance(ExtractIlluminaBarcodes.class);

    /**
     * The read structure of the actual Illumina Run, i.e. the readStructure of the input data
     */
    private ReadStructure readStructure;

    private IlluminaDataProviderFactory factory;

    private final Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();

    private final NumberFormat tileNumberFormatter = getNumberInstance();
    private BclQualityEvaluationStrategy bclQualityEvaluationStrategy;

    public ExtractIlluminaBarcodes() {
        tileNumberFormatter.setMinimumIntegerDigits(4);
        tileNumberFormatter.setGroupingUsed(false);
    }

    @Override
    protected Object doWork() {

        assertFileIsWritable(METRICS_FILE);
        if (OUTPUT_DIR == null) {
            OUTPUT_DIR = BASECALLS_DIR;
        }
        assertDirectoryIsWritable(OUTPUT_DIR);

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.barcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == Barcode) {
                noMatchBarcode[index++] = repeatCharNTimes('N', d.length);
            }
        }

        final BarcodeMetric noMatchMetric = new BarcodeMetric(null, null, barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

        final int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            numProcessors = getRuntime().availableProcessors() + NUM_PROCESSORS;
        } else {
            numProcessors = NUM_PROCESSORS;
        }

        LOG.info("Processing with " + numProcessors + " PerTileBarcodeExtractor(s).");
        final ExecutorService pool = newFixedThreadPool(numProcessors);

        // TODO: This is terribly inefficient; we're opening a huge number of files via the extractor constructor and we never close them.
        final List<PerTileBarcodeExtractor> extractors = new ArrayList<>(factory.getAvailableTiles().size());
        for (final int tile : factory.getAvailableTiles()) {
            final PerTileBarcodeExtractor extractor = new PerTileBarcodeExtractor(
                    tile,
                    getBarcodeFile(tile),
                    barcodeToMetrics,
                    noMatchMetric,
                    factory,
                    MINIMUM_BASE_QUALITY,
                    MAX_NO_CALLS,
                    MAX_MISMATCHES,
                    MIN_MISMATCH_DELTA
            );
            extractors.add(extractor);
        }
        try {
            for (final PerTileBarcodeExtractor extractor : extractors) {
                pool.submit(extractor);
            }
            pool.shutdown();
            // Wait a while for existing tasks to terminate
            if (!pool.awaitTermination(6, HOURS)) {
                pool.shutdownNow(); // Cancel any still-executing tasks
                // Wait a while for tasks to respond to being cancelled
                if (!pool.awaitTermination(60, SECONDS))
                    LOG.error("Pool did not terminate");
                return 1;
            }
        } catch (final Throwable e) {
            // (Re-)Cancel if current thread also interrupted
            LOG.error(e, "Parent thread encountered problem submitting extractors to thread pool or awaiting shutdown of threadpool.  Attempting to kill threadpool.");
            pool.shutdownNow();
            return 2;
        }

        LOG.info("Processed " + extractors.size() + " tiles.");
        for (final PerTileBarcodeExtractor extractor : extractors) {
            for (final String key : barcodeToMetrics.keySet()) {
                barcodeToMetrics.get(key).merge(extractor.getMetrics().get(key));
            }
            noMatchMetric.merge(extractor.getNoMatchMetric());
            if (extractor.getException() != null) {
                LOG.error("Abandoning metrics calculation because one or more PerTileBarcodeExtractors failed.");
                return 4;
            }
        }

        // Finish metrics tallying.
        int totalReads = noMatchMetric.READS;
        int totalPfReads = noMatchMetric.PF_READS;
        int totalPfReadsAssigned = 0;
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            totalReads += barcodeMetric.READS;
            totalPfReads += barcodeMetric.PF_READS;
            totalPfReadsAssigned += barcodeMetric.PF_READS;
        }

        if (totalReads > 0) {
            noMatchMetric.PCT_MATCHES = noMatchMetric.READS / (double) totalReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PCT_MATCHES = barcodeMetric.READS / (double) totalReads;
                if (barcodeMetric.PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PCT_MATCHES / bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PCT_MATCHES / bestPctOfAllBarcodeMatches;
                }
            }
        }

        if (totalPfReads > 0) {
            noMatchMetric.PF_PCT_MATCHES = noMatchMetric.PF_READS / (double) totalPfReads;
            double bestPctOfAllBarcodeMatches = 0;
            for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                barcodeMetric.PF_PCT_MATCHES = barcodeMetric.PF_READS / (double) totalPfReads;
                if (barcodeMetric.PF_PCT_MATCHES > bestPctOfAllBarcodeMatches) {
                    bestPctOfAllBarcodeMatches = barcodeMetric.PF_PCT_MATCHES;
                }
            }
            if (bestPctOfAllBarcodeMatches > 0) {
                noMatchMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                        noMatchMetric.PF_PCT_MATCHES / bestPctOfAllBarcodeMatches;
                for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
                    barcodeMetric.PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT =
                            barcodeMetric.PF_PCT_MATCHES / bestPctOfAllBarcodeMatches;
                }
            }
        }

        // Warn about minimum qualities and assert that we've achieved the minimum.
        for (Map.Entry<Byte, Integer> entry : bclQualityEvaluationStrategy.getPoorQualityFrequencies().entrySet()) {
            LOG.warn(format("Observed low quality of %s %s times.", entry.getKey(), entry.getValue()));
        }
        bclQualityEvaluationStrategy.assertMinimumQualities();

        // Calculate the normalized matches
        if (totalPfReadsAssigned > 0) {
            final double mean = (double) totalPfReadsAssigned / (double) barcodeToMetrics.values().size();
            for (final BarcodeMetric m : barcodeToMetrics.values()) {
                m.PF_NORMALIZED_MATCHES = m.PF_READS / mean;
            }
        }

        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : barcodeToMetrics.values()) {
            metrics.addMetric(barcodeMetric);
        }
        metrics.addMetric(noMatchMetric);
        metrics.write(METRICS_FILE);
        return 0;
    }

    /**
     * Create a barcode filename corresponding to the given tile qseq file.
     */
    private File getBarcodeFile(final int tile) {
        return new File(OUTPUT_DIR,
                "s_" + LANE + "_" + tileNumberFormatter.format(tile) + "_barcode.txt" + (COMPRESS_OUTPUTS ? ".gz" : ""));
    }

    /**
     * Validate that POSITION >= 1, and that all BARCODEs are the same length and unique
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     * to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> messages = new ArrayList<>();

        this.bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(MINIMUM_QUALITY);

        /**
         * In extract illumina barcodes we NEVER want to look at the template reads, therefore replace them with skips because
         * IlluminaDataProvider and its factory will not open these nor produce ClusterData with the template reads in them, thus reducing
         * the file IO and value copying done by the data provider
         */
        readStructure = new ReadStructure(READ_STRUCTURE.replaceAll("T", "S"));
        final IlluminaDataType[] datatypes = (MINIMUM_BASE_QUALITY > 0) ?
                new IlluminaDataType[]{BaseCalls, PF, QualityScores} :
                new IlluminaDataType[]{BaseCalls, PF};
        factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy, datatypes);

        if (BARCODE_FILE != null) {
            parseBarcodeFile(messages);
        } else {
            final Set<String> barcodes = new HashSet<>();
            for (final String barcode : BARCODE) {
                if (barcodes.contains(barcode)) {
                    messages.add("Barcode " + barcode + " specified more than once.");
                }
                barcodes.add(barcode);
                final BarcodeMetric metric = new BarcodeMetric(null, null, barcode, new String[]{barcode});
                barcodeToMetrics.put(barcode, metric);
            }
        }
        if (barcodeToMetrics.keySet().size() == 0) {
            messages.add("No barcodes have been specified.");
        }
        if (messages.size() == 0) {
            return null;
        }
        return messages.toArray(new String[messages.size()]);
    }

    private static final String BARCODE_SEQUENCE_COLUMN = "barcode_sequence";
    private static final String BARCODE_SEQUENCE_1_COLUMN = "barcode_sequence_1";
    private static final String BARCODE_NAME_COLUMN = "barcode_name";
    private static final String LIBRARY_NAME_COLUMN = "library_name";

    private void parseBarcodeFile(final ArrayList<String> messages) {
        final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(BARCODE_FILE);
        final String sequenceColumn = barcodesParser.hasColumn(BARCODE_SEQUENCE_COLUMN)
                ? BARCODE_SEQUENCE_COLUMN : barcodesParser.hasColumn(BARCODE_SEQUENCE_1_COLUMN)
                ? BARCODE_SEQUENCE_1_COLUMN : null;
        if (sequenceColumn == null) {
            messages.add(BARCODE_FILE + " does not have " + BARCODE_SEQUENCE_COLUMN + " or " +
                    BARCODE_SEQUENCE_1_COLUMN + " column header");
            return;
        }
        final boolean hasBarcodeName = barcodesParser.hasColumn(BARCODE_NAME_COLUMN);
        final boolean hasLibraryName = barcodesParser.hasColumn(LIBRARY_NAME_COLUMN);
        final int numBarcodes = readStructure.barcodes.length();
        final Set<String> barcodes = new HashSet<>();
        for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
            final String bcStrings[] = new String[numBarcodes];
            int barcodeNum = 1;
            for (final ReadDescriptor rd : readStructure.descriptors) {
                if (rd.type != Barcode) continue;
                final String header = barcodeNum == 1 ? sequenceColumn : "barcode_sequence_" + valueOf(barcodeNum);
                bcStrings[barcodeNum - 1] = row.getField(header);
                barcodeNum++;
            }
            final String bcStr = barcodeSeqsToString(bcStrings);
            if (barcodes.contains(bcStr)) {
                messages.add("Barcode " + bcStr + " specified more than once in " + BARCODE_FILE);
            }
            barcodes.add(bcStr);
            final String barcodeName = (hasBarcodeName ? row.getField(BARCODE_NAME_COLUMN) : "");
            final String libraryName = (hasLibraryName ? row.getField(LIBRARY_NAME_COLUMN) : "");
            final BarcodeMetric metric = new BarcodeMetric(barcodeName, libraryName, bcStr, bcStrings);
            barcodeToMetrics.put(join("", bcStrings), metric);
        }
        barcodesParser.close();
    }

    /**
     * Metrics produced by the ExtractIlluminaBarcodes program that is used to parse data in
     * the basecalls directory and determine to which barcode each read should be assigned.
     */
    public static class BarcodeMetric extends MetricBase {
        /**
         * The barcode (from the set of expected barcodes) for which the following metrics apply.
         * Note that the "symbolic" barcode of NNNNNN is used to report metrics for all reads that
         * do not match a barcode.
         */
        public String BARCODE;
        public String BARCODE_NAME = "";
        public String LIBRARY_NAME = "";
        /**
         * The total number of reads matching the barcode.
         */
        public int READS = 0;
        /**
         * The number of PF reads matching this barcode (always less than or equal to READS).
         */
        public int PF_READS = 0;
        /**
         * The number of all reads matching this barcode that matched with 0 errors or no-calls.
         */
        public int PERFECT_MATCHES = 0;
        /**
         * The number of PF reads matching this barcode that matched with 0 errors or no-calls.
         */
        public int PF_PERFECT_MATCHES = 0;
        /**
         * The number of all reads matching this barcode that matched with 1 error or no-call.
         */
        public int ONE_MISMATCH_MATCHES = 0;
        /**
         * The number of PF reads matching this barcode that matched with 1 error or no-call.
         */
        public int PF_ONE_MISMATCH_MATCHES = 0;
        /**
         * The percentage of all reads in the lane that matched to this barcode.
         */
        public double PCT_MATCHES = 0d;
        /**
         * The rate of all reads matching this barcode to all reads matching the most prevelant barcode. For the
         * most prevelant barcode this will be 1, for all others it will be less than 1 (except for the possible
         * exception of when there are more orphan reads than for any other barcode, in which case the value
         * may be arbitrarily large).  One over the lowest number in this column gives you the fold-difference
         * in representation between barcodes.
         */
        public double RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;
        /**
         * The percentage of PF reads in the lane that matched to this barcode.
         */
        public double PF_PCT_MATCHES = 0d;

        /**
         * The rate of PF reads matching this barcode to PF reads matching the most prevelant barcode. For the
         * most prevelant barcode this will be 1, for all others it will be less than 1 (except for the possible
         * exception of when there are more orphan reads than for any other barcode, in which case the value
         * may be arbitrarily large).  One over the lowest number in this column gives you the fold-difference
         * in representation of PF reads between barcodes.
         */
        public double PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT = 0d;

        /**
         * The "normalized" matches to each barcode. This is calculated as the number of pf reads matching
         * this barcode over the sum of all pf reads matching any barcode (excluding orphans). If all barcodes
         * are represented equally this will be 1.
         */
        public double PF_NORMALIZED_MATCHES;

        protected byte[][] barcodeBytes;

        public BarcodeMetric(final String barcodeName, final String libraryName,
                             final String barcodeDisplay, final String[] barcodeSeqs) {

            this.BARCODE = barcodeDisplay;
            this.BARCODE_NAME = barcodeName;
            this.LIBRARY_NAME = libraryName;
            this.barcodeBytes = new byte[barcodeSeqs.length][];
            for (int i = 0; i < barcodeSeqs.length; i++) {
                barcodeBytes[i] = stringToBytes(barcodeSeqs[i]);
            }
        }

        /**
         * This ctor is necessary for when reading metrics from file
         */
        public BarcodeMetric() {
            barcodeBytes = null;
        }

        /**
         * Creates a copy of metric initialized with only non-accumulated and non-calculated values set
         */
        public static BarcodeMetric copy(final BarcodeMetric metric) {
            final BarcodeMetric result = new BarcodeMetric();
            result.BARCODE = metric.BARCODE;
            result.BARCODE_NAME = metric.BARCODE_NAME;
            result.LIBRARY_NAME = metric.LIBRARY_NAME;
            result.barcodeBytes = metric.barcodeBytes;
            return result;
        }

        /**
         * Adds the non-calculated
         *
         * @param metric
         */
        public void merge(final BarcodeMetric metric) {
            this.READS += metric.READS;
            this.PF_READS += metric.PF_READS;
            this.PERFECT_MATCHES += metric.PERFECT_MATCHES;
            this.PF_PERFECT_MATCHES += metric.PF_PERFECT_MATCHES;
            this.ONE_MISMATCH_MATCHES += metric.ONE_MISMATCH_MATCHES;
            this.PF_ONE_MISMATCH_MATCHES += metric.PF_ONE_MISMATCH_MATCHES;
        }

    }

    /**
     * Extracts barcodes and accumulates metrics for an entire tile.
     */
    private static class PerTileBarcodeExtractor implements Runnable {
        private final int tile;
        private final File barcodeFile;
        private final Map<String, BarcodeMetric> metrics;
        private final BarcodeMetric noMatch;
        private Exception exception = null;
        private final boolean usingQualityScores;
        private final IlluminaDataProvider provider;
        private final ReadStructure outputReadStructure;
        private final int maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality;

        /**
         * Utility class to hang onto data about the best match for a given barcode
         */
        class BarcodeMatch {
            boolean matched;
            String barcode;
            int mismatches;
            int mismatchesToSecondBest;
        }

        /**
         * Constructor
         *
         * @param tile             The number of the tile being processed; used for logging only.
         * @param barcodeFile      The file to write the barcodes to
         * @param noMatchMetric    A "template" metric that is cloned and the clone is stored internally for accumulating data
         * @param barcodeToMetrics A "template" metric map whose metrics are cloned, and the clones are stored internally for accumulating data
         */
        public PerTileBarcodeExtractor(
                final int tile,
                final File barcodeFile,
                final Map<String, BarcodeMetric> barcodeToMetrics,
                final BarcodeMetric noMatchMetric,
                final IlluminaDataProviderFactory factory,
                final int minimumBaseQuality,
                final int maxNoCalls,
                final int maxMismatches,
                final int minMismatchDelta
        ) {
            this.tile = tile;
            this.barcodeFile = barcodeFile;
            this.usingQualityScores = minimumBaseQuality > 0;
            this.maxNoCalls = maxNoCalls;
            this.maxMismatches = maxMismatches;
            this.minMismatchDelta = minMismatchDelta;
            this.minimumBaseQuality = minimumBaseQuality;
            this.metrics = new LinkedHashMap<>(barcodeToMetrics.size());
            for (final String key : barcodeToMetrics.keySet()) {
                this.metrics.put(key, copy(barcodeToMetrics.get(key)));
            }
            this.noMatch = copy(noMatchMetric);
            this.provider = factory.makeDataProvider(asList(tile));
            this.outputReadStructure = factory.getOutputReadStructure();

        }

        // These methods return the results of the extraction
        public synchronized Map<String, BarcodeMetric> getMetrics() {
            return this.metrics;
        }

        public synchronized BarcodeMetric getNoMatchMetric() {
            return this.noMatch;
        }

        public synchronized Exception getException() {
            return this.exception;
        }

        /**
         * run method which extracts barcodes and accumulates metrics for an entire tile
         */
        synchronized public void run() {
            try {
                LOG.info("Extracting barcodes for tile " + tile);

                //Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                //is non-overlapping sets of files so make the  data providers in the individual threads for PerTileBarcodeExtractors
                //so they are not all waiting for each others file operations

                //Most likely we have SKIPS in our read structure since we replace all template reads with skips in the input data structure
                //(see customCommnandLineValidation), therefore we must use the outputReadStructure to index into the output cluster data
                final int[] barcodeIndices = outputReadStructure.barcodes.getIndices();
                final BufferedWriter writer = openFileForBufferedWriting(barcodeFile);
                final byte barcodeSubsequences[][] = new byte[barcodeIndices.length][];
                final byte qualityScores[][] = usingQualityScores ? new byte[barcodeIndices.length][] : null;
                while (provider.hasNext()) {
                    // Extract the barcode from the cluster and write it to the file for the tile
                    final ClusterData cluster = provider.next();
                    for (int i = 0; i < barcodeIndices.length; i++) {
                        barcodeSubsequences[i] = cluster.getRead(barcodeIndices[i]).getBases();
                        if (usingQualityScores) qualityScores[i] = cluster.getRead(barcodeIndices[i]).getQualities();
                    }
                    final boolean passingFilter = cluster.isPf();
                    final BarcodeMatch match = findBestBarcodeAndUpdateMetrics(barcodeSubsequences, qualityScores, passingFilter, metrics, noMatch);

                    final String yOrN = (match.matched ? "Y" : "N");

                    for (final byte[] bc : barcodeSubsequences) {
                        writer.write(bytesToString(bc));
                    }
                    writer.write("\t" + yOrN + "\t" + match.barcode + "\t" + valueOf(match.mismatches) +
                            "\t" + valueOf(match.mismatchesToSecondBest));
                    writer.newLine();
                }
                writer.close();
            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                provider.close();
            }
        }

        /**
         * Find the best barcode match for the given read sequence, and accumulate metrics
         *
         * @param readSubsequences portion of read containing barcode
         * @param passingFilter    PF flag for the current read
         * @return perfect barcode string, if there was a match within tolerance, or null if not.
         */
        private BarcodeMatch findBestBarcodeAndUpdateMetrics(final byte[][] readSubsequences,
                                                             final byte[][] qualityScores,
                                                             final boolean passingFilter,
                                                             final Map<String, BarcodeMetric> metrics,
                                                             final BarcodeMetric noMatchBarcodeMetric) {
            BarcodeMetric bestBarcodeMetric = null;
            int totalBarcodeReadBases = 0;
            int numNoCalls = 0; // NoCalls are calculated for all the barcodes combined

            for (final byte[] bc : readSubsequences) {
                totalBarcodeReadBases += bc.length;
                for (final byte b : bc) if (isNoCall(b)) ++numNoCalls;
            }

            // PIC-506 When forcing all reads to match a single barcode, allow a read to match even if every
            // base is a mismatch.
            int numMismatchesInBestBarcode = totalBarcodeReadBases + 1;
            int numMismatchesInSecondBestBarcode = totalBarcodeReadBases + 1;

            for (final BarcodeMetric barcodeMetric : metrics.values()) {
                final int numMismatches = countMismatches(barcodeMetric.barcodeBytes, readSubsequences, qualityScores);
                if (numMismatches < numMismatchesInBestBarcode) {
                    if (bestBarcodeMetric != null) {
                        numMismatchesInSecondBestBarcode = numMismatchesInBestBarcode;
                    }
                    numMismatchesInBestBarcode = numMismatches;
                    bestBarcodeMetric = barcodeMetric;
                } else if (numMismatches < numMismatchesInSecondBestBarcode) {
                    numMismatchesInSecondBestBarcode = numMismatches;
                }
            }

            final boolean matched = bestBarcodeMetric != null &&
                    numNoCalls <= maxNoCalls &&
                    numMismatchesInBestBarcode <= maxMismatches &&
                    numMismatchesInSecondBestBarcode - numMismatchesInBestBarcode >= minMismatchDelta;

            final BarcodeMatch match = new BarcodeMatch();

            // If we have something that's not a "match" but matches one barcode
            // slightly, we output that matching barcode in lower case
            if (numNoCalls + numMismatchesInBestBarcode < totalBarcodeReadBases) {
                match.mismatches = numMismatchesInBestBarcode;
                match.mismatchesToSecondBest = numMismatchesInSecondBestBarcode;
                match.barcode = bestBarcodeMetric.BARCODE.toLowerCase().replaceAll(BARCODE_DELIMITER, "");
            } else {
                match.mismatches = totalBarcodeReadBases;
                match.barcode = "";
            }

            if (matched) {
                ++bestBarcodeMetric.READS;
                if (passingFilter) {
                    ++bestBarcodeMetric.PF_READS;
                }
                if (numMismatchesInBestBarcode == 0) {
                    ++bestBarcodeMetric.PERFECT_MATCHES;
                    if (passingFilter) {
                        ++bestBarcodeMetric.PF_PERFECT_MATCHES;
                    }
                } else if (numMismatchesInBestBarcode == 1) {
                    ++bestBarcodeMetric.ONE_MISMATCH_MATCHES;
                    if (passingFilter) {
                        ++bestBarcodeMetric.PF_ONE_MISMATCH_MATCHES;
                    }
                }

                match.matched = true;
                match.barcode = bestBarcodeMetric.BARCODE.replaceAll(BARCODE_DELIMITER, "");
            } else {
                ++noMatchBarcodeMetric.READS;
                if (passingFilter) {
                    ++noMatchBarcodeMetric.PF_READS;
                }
            }

            return match;
        }

        /**
         * Compare barcode sequence to bases from read
         *
         * @return how many bases did not match
         */
        private int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence, final byte[][] qualities) {
            int numMismatches = 0;
            // Read sequence and barcode length may not be equal, so we just use the shorter of the two
            for (int j = 0; j < barcodeBytes.length; j++) {
                final int basesToCheck = min(barcodeBytes[j].length, readSubsequence[j].length);
                for (int i = 0; i < basesToCheck; ++i) {
                    if (!isNoCall(readSubsequence[j][i])) {
                        if (!basesEqual(barcodeBytes[j][i], readSubsequence[j][i])) ++numMismatches;
                        else if (qualities != null && qualities[j][i] < minimumBaseQuality) ++numMismatches;
                    }
                }
            }
            return numMismatches;
        }
    }
}
