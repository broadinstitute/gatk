package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ClusterData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProvider;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProviderFactory;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static htsjdk.samtools.util.IOUtil.assertFileIsWritable;
import static htsjdk.samtools.util.Log.getInstance;
import static java.lang.Long.MAX_VALUE;
import static java.lang.Runtime.getRuntime;
import static java.util.Arrays.asList;
import static java.util.concurrent.Executors.newFixedThreadPool;
import static java.util.concurrent.TimeUnit.DAYS;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.LANE_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.CollectHiSeqXPfFailMetrics.ReadClassifier.PfFailReason.EMPTY;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.CollectHiSeqXPfFailMetrics.ReadClassifier.PfFailReason.MISALIGNED;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.CollectHiSeqXPfFailMetrics.ReadClassifier.PfFailReason.POLYCLONAL;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.CollectHiSeqXPfFailMetrics.ReadClassifier.PfFailReason.UNKNOWN;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.BaseCalls;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.PF;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Position;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.QualityScores;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

/**
 * Collect metrics regarding the reason for reads (sequenced by HiSeqX) not passing the Illumina PF Filter. (BETA)
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Classify PF-Failing reads in a HiSeqX Illumina Basecalling directory into various categories. The classification is based on a heuristic that was derived by looking at a few titration experiments.",
        usageShort = "Classify PF-Failing reads in a HiSeqX Illumina Basecalling directory into various categories.",
        programGroup = QCProgramGroup.class
)

public class CollectHiSeqXPfFailMetrics extends PicardCommandLineProgram {
    @Argument(doc = "The Illumina basecalls directory. ", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(shortName = OUTPUT_SHORT_NAME, doc = "Basename for metrics file. Resulting file will be" +
            " <OUTPUT>" + summaryMetricsExtension, optional = false)
    public File OUTPUT;

    @Argument(shortName = "P", doc = "The fraction of (non-PF) reads for which to output explicit classification. Output file will be <OUTPUT>" + detailedMetricsExtension + " (if PROB_EXPLICIT_READS != 0)", optional = true)
    public double PROB_EXPLICIT_READS = 0;

    @Argument(doc = "Lane number.", shortName = LANE_SHORT_NAME)
    public Integer LANE;

    @Argument(shortName = "NP", doc = "Run this many PerTileBarcodeExtractors in parallel.  If NUM_PROCESSORS = 0, number of cores is automatically set to " +
            "the number of cores available on the machine. If NUM_PROCESSORS < 0 then the number of cores used will be " +
            "the number available on the machine less NUM_PROCESSORS.", optional = true)
    public int NUM_PROCESSORS = 1;

    @Argument(doc = "Number of cycles to look at. At time of writing PF status gets determined at cycle 24 so numbers greater than this will yield strange results. " +
            "In addition, PF status is currently determined at cycle 24, so running this with any other value is neither tested nor recommended.", optional = true)
    public int N_CYCLES = 24;

    private static final Log LOG = getInstance(CollectHiSeqXPfFailMetrics.class);

    private final Map<Integer, PFFailSummaryMetric> tileToSummaryMetrics = new LinkedHashMap<>();
    private final Map<Integer, List<PFFailDetailedMetric>> tileToDetailedMetrics = new LinkedHashMap<>();

    //Add "T" to the number of cycles to create a "TemplateRead" of the desired length.
    private final ReadStructure READ_STRUCTURE = new ReadStructure(N_CYCLES + "T");

    public final static String detailedMetricsExtension = ".pffail_detailed_metrics";
    public final static String summaryMetricsExtension = ".pffail_summary_metrics";

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();

        if (N_CYCLES < 0) {
            errors.add("Number of Cycles to look at must be greater than 0");
        }

        if (PROB_EXPLICIT_READS > 1 || PROB_EXPLICIT_READS < 0) {
            errors.add("PROB_EXPLICIT_READS must be a probability, i.e., 0 <= PROB_EXPLICIT_READS <= 1");
        }

        if (errors.size() > 0) {
            return errors.toArray(new String[errors.size()]);
        } else {
            return super.customCommandLineValidation();
        }
    }

    @Override
    protected Object doWork() {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, READ_STRUCTURE,
                new BclQualityEvaluationStrategy(ILLUMINA_ALLEGED_MINIMUM_QUALITY),
                BaseCalls,
                PF,
                QualityScores,
                Position);

        final File summaryMetricsFileName = new File(OUTPUT + summaryMetricsExtension);
        final File detailedMetricsFileName = new File(OUTPUT + detailedMetricsExtension);

        assertFileIsWritable(summaryMetricsFileName);
        if (PROB_EXPLICIT_READS != 0) {
            assertFileIsWritable(detailedMetricsFileName);
        }

        final int numProcessors;
        if (NUM_PROCESSORS == 0) {
            numProcessors = getRuntime().availableProcessors();
        } else if (NUM_PROCESSORS < 0) {
            numProcessors = getRuntime().availableProcessors() + NUM_PROCESSORS;
        } else {
            numProcessors = NUM_PROCESSORS;
        }

        // Create thread-pool submit jobs and what for their completion
        LOG.info("Processing with " + numProcessors + " PerTilePFMetricsExtractor(s).");
        final ExecutorService pool = newFixedThreadPool(numProcessors);

        final List<PerTilePFMetricsExtractor> extractors = new ArrayList<>(factory.getAvailableTiles().size());
        for (final int tile : factory.getAvailableTiles()) {
            tileToSummaryMetrics.put(tile, new PFFailSummaryMetric(Integer.toString(tile)));
            tileToDetailedMetrics.put(tile, new ArrayList<>());

            final PerTilePFMetricsExtractor extractor = new PerTilePFMetricsExtractor(
                    tile,
                    tileToSummaryMetrics.get(tile),
                    tileToDetailedMetrics.get(tile),
                    factory,
                    PROB_EXPLICIT_READS
            );
            extractors.add(extractor);
        }
        try {
            for (final PerTilePFMetricsExtractor extractor : extractors) {
                pool.submit(extractor);
            }
            pool.shutdown();
            // Wait forever for tasks to terminate
            pool.awaitTermination(MAX_VALUE, DAYS);
        } catch (final Throwable e) {
            // Cancel if current thread also interrupted
            LOG.error(e, "Parent thread encountered problem submitting extractors to thread pool or awaiting shutdown of threadpool.  Attempting to kill threadpool.");
            pool.shutdownNow();
            throw new GATKException(e.getMessage());
        }

        LOG.info("Processed " + extractors.size() + " tiles.");

        // Check for exceptions from extractors
        for (final PerTilePFMetricsExtractor extractor : extractors) {
            if (extractor.getException() != null) {
                LOG.error("Abandoning metrics calculation because one or more PerTilePFMetricsExtractors failed.");
                throw new GATKException(extractor.getException().getMessage());
            }
        }

        // Add detailed metrics to file
        final MetricsFile<PFFailDetailedMetric, ?> detailedMetrics = getMetricsFile();
        for (final Collection<PFFailDetailedMetric> detailedMetricCollection : tileToDetailedMetrics.values()) {
            for (final PFFailDetailedMetric metric : detailedMetricCollection) {
                detailedMetrics.addMetric(metric);
            }
        }

        // If detailed metrics were requested, write them now.
        if (PROB_EXPLICIT_READS > 0) {
            detailedMetrics.write(detailedMetricsFileName);
        }

        // Finish metrics tallying. Looping twice so that the "All" metrics will come out on top.
        final PFFailSummaryMetric totalMetric = new PFFailSummaryMetric("All"); // a "fake" tile that will contain the total tally
        for (final PFFailSummaryMetric summaryMetric : tileToSummaryMetrics.values()) {
            totalMetric.merge(summaryMetric);
        }

        // Derive fields for total metric and add to file
        totalMetric.calculateDerivedFields();
        final MetricsFile<PFFailSummaryMetric, ?> summaryMetricsFile = getMetricsFile();
        summaryMetricsFile.addMetric(totalMetric);

        // Prepare each tile's derived fields and add it to the file
        for (final PFFailSummaryMetric summaryMetric : tileToSummaryMetrics.values()) {
            summaryMetric.calculateDerivedFields();
            summaryMetricsFile.addMetric(summaryMetric);
        }

        // Actually write the summary metrics to their file.
        summaryMetricsFile.write(summaryMetricsFileName);

        return null;
    }

    /**
     * Extracts metrics from a HiSeqX tile.
     */
    private static class PerTilePFMetricsExtractor implements Runnable {

        private final int tile;
        private final PFFailSummaryMetric summaryMetric;
        final Collection<PFFailDetailedMetric> detailedMetrics;
        private Exception exception = null;
        private final IlluminaDataProvider provider;
        final private double pWriteDetailed;
        final private Random random = new Random();

        /**
         * Constructor
         *
         * @param tile            The number of the tile being processed.
         * @param summaryMetric   A summaryMetric for collecting the tile data in.
         * @param detailedMetrics A set of metrics for collecting the classification data in.
         * @param factory         A dataprovider for IlluminaData
         */
        public PerTilePFMetricsExtractor(
                final int tile,
                final PFFailSummaryMetric summaryMetric,
                final Collection<PFFailDetailedMetric> detailedMetrics,
                final IlluminaDataProviderFactory factory,
                final double pWriteDetailed
        ) {
            this.tile = tile;
            this.summaryMetric = summaryMetric;
            this.detailedMetrics = detailedMetrics;
            this.pWriteDetailed = pWriteDetailed;
            this.provider = factory.makeDataProvider(asList(tile));
        }

        public Exception getException() {
            return this.exception;
        }

        /**
         * run method which extracts accumulates metrics for a tile
         */
        public void run() {
            try {
                LOG.info("Extracting PF metrics for tile " + tile);

                /**
                 *   Sometimes makeDataProvider takes a while waiting for slow file IO, for each tile the needed set of files
                 *   is non-overlapping sets of files so make the data providers in the individual threads for Extractors
                 *   so they are not all waiting for each others file operations
                 */
                while (provider.hasNext()) {
                    // Extract the PF status and infer reason if FAIL from the cluster and update the summaryMetric for the tile
                    final ClusterData cluster = provider.next();
                    this.summaryMetric.READS++;
                    if (!cluster.isPf()) {
                        this.summaryMetric.PF_FAIL_READS++;

                        final ReadClassifier readClassifier = new ReadClassifier(cluster.getRead(0));

                        if (random.nextDouble() < pWriteDetailed) {
                            detailedMetrics.add(new PFFailDetailedMetric(tile, cluster.getX(), cluster.getY(), readClassifier.numNs, readClassifier.numQGtTwo, readClassifier.failClass));
                        }
                        switch (readClassifier.failClass) {
                            case EMPTY:
                                this.summaryMetric.PF_FAIL_EMPTY++;
                                break;
                            case MISALIGNED:
                                this.summaryMetric.PF_FAIL_MISALIGNED++;
                                break;
                            case POLYCLONAL:
                                this.summaryMetric.PF_FAIL_POLYCLONAL++;
                                break;
                            case UNKNOWN:
                                this.summaryMetric.PF_FAIL_UNKNOWN++;
                                break;
                            default:
                                LOG.error("Got unexpected fail Reason");
                        }
                    }
                }
            } catch (final Exception e) {
                LOG.error(e, "Error processing tile ", this.tile);
                this.exception = e;
            } finally {
                provider.close();
            }
        }
    }

    protected static class ReadClassifier {
        public enum PfFailReason {
            EMPTY,
            POLYCLONAL,
            MISALIGNED,
            UNKNOWN
        }

        private final int numNs; // The number of Ns in the base calls
        private final int numQGtTwo; // The number of quality scores greater than 2
        private PfFailReason failClass = null; // The classification of the failure mode

        /**
         * Heart of CLP.
         * This class actually classifies ReadData into the reason why it failed PF
         * classification is based on a small set of titrated flowcells sequenced at the Broad Institute by the Genomics Platform.
         * Three cluster were observed:
         * - numNs~24 and was found only near the boundaries of tiles. it didn't seem to depend on concentration. For this reason it
         * was classified as MISALIGNED
         * <p>
         * - numNs~0 and numQGtTwo<=8 these were found throughout the tiles and _decreased_ in number as the concentration of the library increased
         * Thus it was concluded that these correspond to the EMPTY wells
         * <p>
         * - numNs~0 and numQGtTwo>=12 there were found throughout the tiles and _increased_ in number as the concentration of the library increased
         * Thus it was concluded that these correspond to the POLYCLONAL wells
         * <p>
         * - the remaining reads were few in number the classification for them wasn't clear. Thus they are left as UNKNOWN.
         * <p>
         * We use the length of the read as a parameter and scale the 8 and the 12 accordingly as length/3 and length/2, but in reality this has only
         * been tested on length=24.
         *
         * @param read The read to classify.
         */
        public ReadClassifier(final ReadData read) {

            final int length = read.getBases().length;

            numNs = countEquals(read.getBases(), (byte) '.'); // Ns are returned as periods from Illumina
            numQGtTwo = countGreaterThan(read.getQualities(), (byte) 2);

            failClass = UNKNOWN; //for cases not covered below
            if (numNs >= (length - 1)) {
                failClass = MISALIGNED;
            } else if (numNs <= 1) {
                if (numQGtTwo <= length / 3) {
                    failClass = EMPTY;
                } else if (numQGtTwo >= length / 2) {
                    failClass = POLYCLONAL;
                }
            }
        }
    }

    /**
     * a metric class for describing FP failing reads from an Illumina HiSeqX lane *
     */
    public static class PFFailDetailedMetric extends MetricBase {
        // The Tile that is described by this metric.
        public Integer TILE;

        //The X coordinate of the read within the tile
        public int X;

        //The Y coordinate of the read within the tile
        public int Y;

        //The number of Ns found in this read.
        public int NUM_N;

        //The number of Quality scores greater than 2 found in this read
        public int NUM_Q_GT_TWO;

        /**
         * The classification of this read: {EMPTY, POLYCLONAL, MISALIGNED, UNKNOWN}
         * (See PFFailSummaryMetric for explanation regarding the possible classification.)
         */
        public ReadClassifier.PfFailReason CLASSIFICATION;

        public PFFailDetailedMetric(final Integer TILE, final int x, final int y, final int NUM_N, final int NUM_Q_GT_TWO, final ReadClassifier.PfFailReason CLASSIFICATION) {
            this.TILE = TILE;
            X = x;
            Y = y;
            this.NUM_N = NUM_N;
            this.NUM_Q_GT_TWO = NUM_Q_GT_TWO;
            this.CLASSIFICATION = CLASSIFICATION;
        }

        /**
         * This ctor is necessary for when reading metrics from file
         */
        public PFFailDetailedMetric() {
        }
    }

    /**
     * Metrics produced by the GetHiSeqXPFFailMetrics program. Used to diagnose lanes from HiSeqX
     * Sequencing, providing the number and fraction of each of the reasons that reads could have not passed PF.
     * Possible reasons are EMPTY (reads from empty wells with no template strand), POLYCLONAL (reads from wells that had more than one strand
     * cloned in them), MISALIGNED (reads from wells that are near the edge of the tile), UNKNOWN (reads that didn't pass PF but couldn't be diagnosed)
     */
    public static class PFFailSummaryMetric extends MetricBase {
        /**
         * The Tile that is described by this metric. Can be a string (like "All") to mean some marginal over tiles. *
         */
        public String TILE = null;

        /**
         * The total number of reads examined
         */
        public int READS = 0;

        /**
         * The number of non-PF reads in this tile.
         */
        public int PF_FAIL_READS = 0;

        /**
         * The fraction of PF_READS
         */
        public double PCT_PF_FAIL_READS = 0.0;

        /**
         * The number of non-PF reads in this tile that are deemed empty.
         */
        public int PF_FAIL_EMPTY = 0;

        /**
         * The fraction of non-PF reads in this tile that are deemed empty (as fraction of all non-PF reads).
         */
        public double PCT_PF_FAIL_EMPTY = 0.0;

        /**
         * The number of non-PF reads in this tile that are deemed multiclonal.
         */
        public int PF_FAIL_POLYCLONAL = 0;

        /**
         * The fraction of non-PF reads in this tile that are deemed multiclonal (as fraction of all non-PF reads).
         */
        public double PCT_PF_FAIL_POLYCLONAL = 0.0;

        /**
         * The number of non-PF reads in this tile that are deemed "misaligned".
         */
        public int PF_FAIL_MISALIGNED = 0;

        /**
         * The fraction of non-PF reads in this tile that are deemed "misaligned" (as fraction of all non-PF reads).
         */
        public double PCT_PF_FAIL_MISALIGNED = 0.0;

        /**
         * The number of non-PF reads in this tile that have not been classified.
         */
        public int PF_FAIL_UNKNOWN = 0;

        /**
         * The fraction of non-PF reads in this tile that have not been classified (as fraction of all non-PF reads).
         */
        public double PCT_PF_FAIL_UNKNOWN = 0.0;

        // constructor takes a String for tile since we want to have one instance with tile="All". This tile will contain the summary of all the tiles
        public PFFailSummaryMetric(final String tile) {
            TILE = tile;
        }

        /**
         * This ctor is necessary for when reading metrics from file
         */
        public PFFailSummaryMetric() {
        }

        /**
         * Adds the non-calculated fields from the other metric to this one.
         *
         * @param metric
         */
        public void merge(final PFFailSummaryMetric metric) {
            this.READS += metric.READS;
            this.PF_FAIL_READS += metric.PF_FAIL_READS;
            this.PF_FAIL_EMPTY += metric.PF_FAIL_EMPTY;
            this.PF_FAIL_MISALIGNED += metric.PF_FAIL_MISALIGNED;
            this.PF_FAIL_POLYCLONAL += metric.PF_FAIL_POLYCLONAL;
            this.PF_FAIL_UNKNOWN += metric.PF_FAIL_UNKNOWN;
        }

        public void calculateDerivedFields() {
            //protect against divide by zero
            if (this.READS != 0) {
                this.PCT_PF_FAIL_READS = (double) this.PF_FAIL_READS / this.READS;
                this.PCT_PF_FAIL_EMPTY = (double) this.PF_FAIL_EMPTY / this.READS;
                this.PCT_PF_FAIL_MISALIGNED = (double) this.PF_FAIL_MISALIGNED / this.READS;
                this.PCT_PF_FAIL_POLYCLONAL = (double) this.PF_FAIL_POLYCLONAL / this.READS;
                this.PCT_PF_FAIL_UNKNOWN = (double) this.PF_FAIL_UNKNOWN / this.READS;
            }
        }
    }

    /**
     * a simple function that counts how many elements in array are equal to 'toCount'
     *
     * @param array
     * @param toCount
     * @return number of elements in array that == 'toCount'
     */
    static private int countEquals(final byte[] array, final byte toCount) {
        int count = 0;
        for (final byte t : array) {
            if (t == toCount) count++;
        }
        return count;
    }

    /**
     * a simple function that counts how many elements in array are greater-than-or-equal-to 'value'
     *
     * @param array
     * @param value
     * @return number of elements in array that >= 'value'
     */
    static private int countGreaterThan(final byte[] array, final byte value) {
        int count = 0;
        for (final int t : array) {
            if (t > value) count++;
        }
        return count;
    }
}
