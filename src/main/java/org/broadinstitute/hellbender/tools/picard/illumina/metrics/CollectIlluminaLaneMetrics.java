package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.Tile;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.TileMetricsUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.partition;
import static htsjdk.samtools.util.Log.getInstance;
import static java.lang.String.format;
import static org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions.OUTPUT_SHORT_NAME;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.CollectIlluminaLaneMetrics.IlluminaLaneMetricsCollector.collectLaneMetrics;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.IlluminaPhasingMetrics.getExtension;
import static org.broadinstitute.hellbender.tools.picard.illumina.metrics.IlluminaPhasingMetrics.getPhasingMetricsForTiles;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure.PARAMETER_DOC;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.TileMetricsUtil.parseTileMetrics;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.TileMetricsUtil.renderTileMetricsFileFromBasecallingDirectory;

/**
 * @author mccowan
 */
@CommandLineProgramProperties(
        usage = CollectIlluminaLaneMetrics.USAGE,
        usageShort = CollectIlluminaLaneMetrics.USAGE,
        programGroup = QCProgramGroup.class
)
public class CollectIlluminaLaneMetrics extends PicardCommandLineProgram {
    static final String USAGE = "Collects Illumina lane metrics for the given basecalling analysis directory";

    @Argument(doc = "The Illumina run directory of the run for which the lane metrics are to be generated")
    public File RUN_DIRECTORY;

    @Argument(doc = "The directory to which the output file will be written")
    public File OUTPUT_DIRECTORY;

    @Argument(doc = "The prefix to be prepended to the file name of the output file; an appropriate suffix will be applied",
            shortName = OUTPUT_SHORT_NAME)
    public String OUTPUT_PREFIX;

    @Argument(doc = PARAMETER_DOC, shortName = "RS")
    public ReadStructure READ_STRUCTURE;

    @Override
    protected Object doWork() {
        final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile = this.getMetricsFile();
        final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile = this.getMetricsFile();
        collectLaneMetrics(RUN_DIRECTORY, OUTPUT_DIRECTORY, OUTPUT_PREFIX, laneMetricsFile, phasingMetricsFile, READ_STRUCTURE);
        return null;
    }

    /**
     * Utility for collating Tile records from the Illumina TileMetrics file into lane-level and phasing-level metrics.
     */
    public static class IlluminaLaneMetricsCollector {

        private final static Log LOG = getInstance(IlluminaLaneMetricsCollector.class);

        /**
         * Returns a partitioned collection of lane number to Tile objects from the provided basecall directory.
         */
        public static Map<Integer, Collection<Tile>> readLaneTiles(final File illuminaRunDirectory, final ReadStructure readStructure) {
            final Collection<Tile> tiles;
            try {
                tiles = parseTileMetrics(renderTileMetricsFileFromBasecallingDirectory(illuminaRunDirectory), readStructure);
            } catch (final FileNotFoundException e) {
                throw new RuntimeIOException("Unable to open laneMetrics file.", e);
            }

            return partition(tiles,
                    new CollectionUtil.Partitioner<Tile, Integer>() {
                        @Override
                        public Integer getPartition(final Tile tile) {
                            return tile.getLaneNumber();
                        }
                    });
        }

        /**
         * Parses the tile data from the basecall directory and writes to both the lane and phasing metrics files
         */
        public static void collectLaneMetrics(final File runDirectory, final File outputDirectory, final String outputPrefix,
                                              final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile,
                                              final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile,
                                              final ReadStructure readStructure) {
            final Map<Integer, Collection<Tile>> laneTiles = readLaneTiles(runDirectory, readStructure);
            writeLaneMetrics(laneTiles, outputDirectory, outputPrefix, laneMetricsFile);
            writePhasingMetrics(laneTiles, outputDirectory, outputPrefix, phasingMetricsFile);
        }

        public static File writePhasingMetrics(final Map<Integer, Collection<Tile>> laneTiles, final File outputDirectory,
                                               final String outputPrefix, final MetricsFile<MetricBase, Comparable<?>> phasingMetricsFile) {
            for (final Integer lane : laneTiles.keySet()) {
                for (final IlluminaPhasingMetrics phasingMetric : getPhasingMetricsForTiles(lane.longValue(), laneTiles.get(lane))) {
                    phasingMetricsFile.addMetric(phasingMetric);
                }
            }

            return writeMetrics(phasingMetricsFile, outputDirectory, outputPrefix, getExtension());
        }

        public static File writeLaneMetrics(final Map<Integer, Collection<Tile>> laneTiles, final File outputDirectory,
                                            final String outputPrefix, final MetricsFile<MetricBase, Comparable<?>> laneMetricsFile) {
            for (final Map.Entry<Integer, Collection<Tile>> entry : laneTiles.entrySet()) {
                final IlluminaLaneMetrics laneMetric = new IlluminaLaneMetrics();
                laneMetric.LANE = entry.getKey().longValue();
                laneMetric.CLUSTER_DENSITY = calculateLaneDensityFromTiles(entry.getValue());
                laneMetricsFile.addMetric(laneMetric);
            }

            return writeMetrics(laneMetricsFile, outputDirectory, outputPrefix, IlluminaLaneMetrics.getExtension());
        }

        private static File writeMetrics(final MetricsFile<MetricBase, Comparable<?>> metricsFile, final File outputDirectory,
                                         final String outputPrefix, final String outputExtension) {
            final File outputFile = new File(outputDirectory, format("%s.%s", outputPrefix, outputExtension));
            LOG.info(format("Writing %s lane metrics to %s ...", metricsFile.getMetrics().size(), outputFile));
            metricsFile.write(outputFile);
            return outputFile;
        }

        private static double calculateLaneDensityFromTiles(final Collection<Tile> tiles) {
            double area = 0;
            double clusters = 0;
            for (final Tile tile : tiles) {
                area += (tile.getClusterCount() / tile.getClusterDensity());
                clusters += tile.getClusterCount();
            }
            return clusters / area;
        }
    }
}
