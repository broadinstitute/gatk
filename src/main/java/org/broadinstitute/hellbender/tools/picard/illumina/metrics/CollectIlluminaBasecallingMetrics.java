package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ClusterData;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProvider;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataProviderFactory;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.broadinstitute.hellbender.utils.text.parsers.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.text.DecimalFormat;
import java.util.*;

import static htsjdk.samtools.util.IOUtil.assertDirectoryIsReadable;
import static htsjdk.samtools.util.IOUtil.assertFileIsReadable;
import static htsjdk.samtools.util.IOUtil.assertFileIsWritable;
import static htsjdk.samtools.util.StringUtil.repeatCharNTimes;
import static java.lang.Double.isNaN;
import static java.lang.Double.valueOf;
import static java.lang.Math.round;
import static java.lang.String.format;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.INPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.LANE_SHORT_NAME;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Barcodes;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.PF;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.Position;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure.PARAMETER_DOC;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY;

/**
 * A Command line tool to collect Illumina Basecalling metrics for a sequencing run
 * Requires a Lane and an input file of Barcodes to expect.
 * Outputs metrics:
 * *  Mean Clusters Per Tile
 * *  Standard Deviation of Clusters Per Tile
 * *  Mean Pf Clusters Per Tile
 * *  Standard Deviation of Pf Clusters Per Tile
 * *  Mean Percentage of Pf Clusters Per Tile
 * *  Standard Deviation of Percentage of Pf Clusters Per Tile
 */
@CommandLineProgramProperties(
        usage = CollectIlluminaBasecallingMetrics.USAGE,
        usageShort = CollectIlluminaBasecallingMetrics.USAGE,
        programGroup = QCProgramGroup.class
)
public class CollectIlluminaBasecallingMetrics extends PicardCommandLineProgram {
    //Command Line Arguments
    static final String USAGE = "Given an Illumina basecalling and a lane, produces per-lane-barcode basecalling metrics";

    @Argument(doc = "The Illumina basecalls output directory from which data are read", shortName = "B")
    public File BASECALLS_DIR;

    @Argument(doc = "The lane whose data will be read", shortName = LANE_SHORT_NAME)
    public Integer LANE;

    // TODO: No longer optional after old workflows are through
    @Argument(doc = "The file containing barcodes to expect from the run - barcodeData.#", shortName = INPUT_SHORT_NAME, optional = true)
    public File INPUT;

    @Argument(doc = PARAMETER_DOC, shortName = "RS")
    public String READ_STRUCTURE;

    @Argument(doc = "The file to which the collected metrics are written", shortName = OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    private int barcodeLength = 0;
    private String unmatched_barcode;
    private final SortedMap<String, IlluminaMetricCounts> barcodeToMetricCounts;

    private static final String BARCODE_NAME_COLUMN = "barcode_name";
    private static final String BARCODE_SEQUENCE_COLUMN_NAME_STUB = "barcode_sequence_";

    public CollectIlluminaBasecallingMetrics() {
        this.barcodeToMetricCounts = new TreeMap<String, IlluminaMetricCounts>();
    }

    @Override
    protected Object doWork() {
        // File and Directory Validation
        assertDirectoryIsReadable(BASECALLS_DIR);
        if (OUTPUT == null) OUTPUT = new File(BASECALLS_DIR, format("LANE%s_basecalling_metrics", LANE));
        assertFileIsWritable(OUTPUT);

        final IlluminaDataProviderFactory factory;
        final ReadStructure readStructure = new ReadStructure(READ_STRUCTURE);
        final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(ILLUMINA_ALLEGED_MINIMUM_QUALITY);

        if (INPUT == null) {
            // TODO: Legacy support. Remove when INPUT is required, after all old workflows are through
            factory = new IlluminaDataProviderFactory(BASECALLS_DIR, LANE, readStructure, bclQualityEvaluationStrategy,
                    PF, Position);
        } else {
            // Grab expected barcode data from barcodeData.<LANE>
            assertFileIsReadable(INPUT);
            final TabbedTextFileWithHeaderParser barcodesParser = new TabbedTextFileWithHeaderParser(INPUT);
            for (final TabbedTextFileWithHeaderParser.Row row : barcodesParser) {
                final String barcodeName = row.getField(BARCODE_NAME_COLUMN);
                final StringBuilder barcode = new StringBuilder();
                for (int i = 1; i <= readStructure.barcodes.length(); i++) {
                    barcode.append(row.getField(BARCODE_SEQUENCE_COLUMN_NAME_STUB + i));
                    if (barcodeLength == 0) barcodeLength = barcode.length();
                }

                // Only add the barcode to the hash if it has sequences. For libraries
                // that don't have barcodes this won't be set in the file.
                if (barcode.length() > 0) {
                    barcodeToMetricCounts.put(barcode.toString(), new IlluminaMetricCounts(barcode.toString(), barcodeName, LANE));
                }
            }

            factory = barcodeToMetricCounts.isEmpty()
                    ? new IlluminaDataProviderFactory(
                    BASECALLS_DIR,
                    LANE,
                    readStructure,
                    bclQualityEvaluationStrategy,
                    PF,
                    Position)
                    : new IlluminaDataProviderFactory(
                    BASECALLS_DIR,
                    LANE,
                    readStructure,
                    bclQualityEvaluationStrategy,
                    PF,
                    Position,
                    Barcodes);
        }

        unmatched_barcode = repeatCharNTimes('N', barcodeLength);

        //Initialize data provider, iterate over clusters, and collect statistics
        final IlluminaDataProvider provider = factory.makeDataProvider();

        while (provider.hasNext()) {
            final ClusterData cluster = provider.next();
            addCluster(cluster);
        }

        onComplete();
        return null;
    }

    /**
     * Process new cluster of Illumina data - increment a running counter of data
     */
    private void addCluster(final ClusterData cluster) {
        //compute hash of Barcode and Lane for key
        String barcode = cluster.getMatchedBarcode();
        if (barcode == null) barcode = unmatched_barcode;

        //increment counts
        IlluminaMetricCounts counters = barcodeToMetricCounts.get(barcode);
        if (counters == null) {
            counters = new IlluminaMetricCounts(barcode, null, LANE);
            barcodeToMetricCounts.put(barcode, counters);
        }
        final int tileNumber = cluster.getTile();
        counters.incrementClusterCount(tileNumber, cluster.isPf());
    }

    /**
     * Handles completion of metric collection. Metrics are computed from counts of data and written out to a file.
     */
    private void onComplete() {
        try {
            final MetricsFile<IlluminaBasecallingMetrics, Comparable<?>> file = getMetricsFile();
            final IlluminaMetricCounts allLaneCounts = new IlluminaMetricCounts(null, null, LANE);
            for (final String s : barcodeToMetricCounts.keySet()) {
                final IlluminaMetricCounts counts = barcodeToMetricCounts.get(s);
                counts.addMetricsToFile(file);
                allLaneCounts.addIlluminaMetricCounts(counts);
            }
            if (!barcodeToMetricCounts.keySet().contains(""))
                allLaneCounts.addMetricsToFile(file);  // detect non-indexed case
            file.write(OUTPUT);
        } catch (final Exception ex) {
            throw new RuntimeIOException("Error writing output file " + OUTPUT.getPath(), ex);
        }
    }

    /**
     * This class manages counts of Illumina Basecalling data on a Per Barcode Per Lane basis.  Cluster and PFCluster
     * counts are stored per tile number.
     */
    private class IlluminaMetricCounts {
        /**
         * Stores counts of clusters found for a specific Barcode-Lane combination across all tiles.  Key = Tile Number, Value = count of clusters**
         */
        private final Histogram<Integer> tileToClusterHistogram;
        /**
         * Stores counts of pf clusters found for a specific Barcode-Lane combination across all tiles.  Key = Tile Number, Value = count of clusters**
         */
        private final Histogram<Integer> tileToPfClusterHistogram;
        final IlluminaBasecallingMetrics metrics;

        public IlluminaMetricCounts(final String barcode, final String barcodeName, final Integer laneNumber) {
            this.tileToClusterHistogram = new Histogram<Integer>();
            this.tileToPfClusterHistogram = new Histogram<Integer>();
            this.metrics = new IlluminaBasecallingMetrics();
            this.metrics.MOLECULAR_BARCODE_SEQUENCE_1 = barcode;
            this.metrics.MOLECULAR_BARCODE_NAME = barcodeName;
            this.metrics.LANE = Integer.toString(laneNumber);
        }

        /*  Increments cluster count by 1 for a given tile number */
        public void incrementClusterCount(final int tileNumber, final boolean isPf) {
            incrementClusterCount(tileNumber, 1d, isPf);
        }

        /*  Increments cluster count by an amount for a given tile number */
        public void incrementClusterCount(final int tileNumber, final double incrementAmount, final boolean isPf) {
            incrementClusterCount(tileNumber, incrementAmount, (isPf ? 1d : 0d));
        }

        /*  Increments cluster count by an amount for a given tile number */
        public void incrementClusterCount(final Integer tileNumber, final double incrementAmount, final double pfIncrementAmount) {
            tileToClusterHistogram.increment(tileNumber, incrementAmount);
            tileToPfClusterHistogram.increment(tileNumber, pfIncrementAmount);
        }

        /* Handles calculating final metrics and updating the metric object */
        private void onComplete() {
            final double meanClustersPerTile = tileToClusterHistogram.getMeanBinSize();
            metrics.MEAN_CLUSTERS_PER_TILE = round(meanClustersPerTile);
            metrics.SD_CLUSTERS_PER_TILE = round(tileToClusterHistogram.getStandardDeviationBinSize(meanClustersPerTile));

            final double meanPfClustersPerTile = tileToPfClusterHistogram.getMeanBinSize();
            metrics.MEAN_PF_CLUSTERS_PER_TILE = round(meanPfClustersPerTile);
            metrics.SD_PF_CLUSTERS_PER_TILE = round(tileToPfClusterHistogram.getStandardDeviationBinSize(meanPfClustersPerTile));

            final DecimalFormat decFormat = new DecimalFormat("#.##");
            final Histogram<Integer> laneToPctPfClusterHistogram = tileToPfClusterHistogram.divideByHistogram(tileToClusterHistogram);
            final double meanPctPfClustersPerTile = laneToPctPfClusterHistogram.getMeanBinSize();
            metrics.MEAN_PCT_PF_CLUSTERS_PER_TILE = (isNaN(meanPctPfClustersPerTile) ? 0 : valueOf(decFormat.format(meanPctPfClustersPerTile * 100)));
            metrics.SD_PCT_PF_CLUSTERS_PER_TILE = valueOf(decFormat.format(laneToPctPfClusterHistogram.getStandardDeviationBinSize(meanPctPfClustersPerTile) * 100));

            metrics.TOTAL_CLUSTERS = (long) this.tileToClusterHistogram.getSumOfValues();
            metrics.PF_CLUSTERS = (long) this.tileToPfClusterHistogram.getSumOfValues();

            final ReadStructure readStructure = new ReadStructure(READ_STRUCTURE);
            int templateBaseCountPerCluster = 0;
            for (int i = 0; i < readStructure.templates.length(); i++)
                templateBaseCountPerCluster += readStructure.templates.get(i).length;
            metrics.TOTAL_READS = metrics.TOTAL_CLUSTERS * readStructure.templates.length();
            metrics.PF_READS = metrics.PF_CLUSTERS * readStructure.templates.length();
            metrics.TOTAL_BASES = metrics.TOTAL_CLUSTERS * templateBaseCountPerCluster;
            metrics.PF_BASES = metrics.PF_CLUSTERS * templateBaseCountPerCluster;

        }

        /* Computes final metric based on data counts and writes to output metric file */
        public void addMetricsToFile(final MetricsFile<IlluminaBasecallingMetrics, Comparable<?>> file) {
            onComplete();
            file.addMetric(metrics);
        }

        /*  Merges data from another IlluminaMetricCount object into current one.*/
        public void addIlluminaMetricCounts(final IlluminaMetricCounts counts) {
            this.tileToClusterHistogram.addHistogram(counts.tileToClusterHistogram);
            this.tileToPfClusterHistogram.addHistogram(counts.tileToPfClusterHistogram);
        }
    }
}
