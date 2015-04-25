package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IterableAdapter;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.TileMetricsOutReader;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.TileMetricsOutReader.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.getSoleElement;
import static htsjdk.samtools.util.CollectionUtil.partition;
import static java.lang.String.format;
import static java.util.Collections.unmodifiableCollection;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaMetricsCode.*;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadType.Template;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.TileTemplateRead.*;

/**
 * Utility for reading the tile data from an Illumina run directory's TileMetricsOut.bin file
 *
 * @author mccowan
 */
public class TileMetricsUtil {
    /**
     * The path to the directory containing the tile metrics file relative to the basecalling directory.
     */
    public static String INTEROP_SUBDIRECTORY_NAME = "InterOp";

    /**
     * The expected name of the tile metrics output file.
     */
    public static String TILE_METRICS_OUT_FILE_NAME = "TileMetricsOut.bin";

    /**
     * Returns the path to the TileMetrics file given the basecalling directory.
     */
    public static File renderTileMetricsFileFromBasecallingDirectory(final File illuminaRunDirectory) {
        return new File(new File(illuminaRunDirectory, INTEROP_SUBDIRECTORY_NAME), TILE_METRICS_OUT_FILE_NAME);
    }

    /**
     * Returns an unmodifiable collection of tile data read from the provided file. For each tile we will extract:
     * - lane number
     * - tile number
     * - density
     * - cluster ID
     * - Phasing & Prephasing for first template read (if available)
     * - Phasing & Prephasing for second template read (if available)
     */
    public static Collection<Tile> parseTileMetrics(final File tileMetricsOutFile, final ReadStructure readStructure) throws FileNotFoundException {
        // Get the tile metrics lines from TileMetricsOut, keeping only the last value for any Lane/Tile/Code combination
        final Collection<IlluminaTileMetrics> tileMetrics = determineLastValueForLaneTileMetricsCode(new TileMetricsOutReader
                (tileMetricsOutFile));

        // Collect the tiles by lane & tile, and then collect the metrics by lane
        final Map<String, Collection<IlluminaTileMetrics>> locationToMetricsMap = partitionTileMetricsByLocation(tileMetrics);
        final Collection<Tile> tiles = new LinkedList<>();
        for (final Map.Entry<String, Collection<IlluminaTileMetrics>> entry : locationToMetricsMap.entrySet()) {
            final Collection<IlluminaTileMetrics> tileRecords = entry.getValue();

            // Get a mapping from metric code number to the corresponding IlluminaTileMetrics
            final Map<Integer, Collection<IlluminaTileMetrics>> codeMetricsMap = partitionTileMetricsByCode(tileRecords);

            final Set<Integer> observedCodes = codeMetricsMap.keySet();
            if (!(observedCodes.contains(DENSITY_ID.getMetricsCode()) && observedCodes.contains(CLUSTER_ID.getMetricsCode())))
                throw new IlluminaParserException(format("Expected to find cluster and density record codes (%s and %s) in records read for tile location %s (lane:tile), but found only %s.",
                        CLUSTER_ID.getMetricsCode(), DENSITY_ID.getMetricsCode(), entry.getKey(), observedCodes));

            final IlluminaTileMetrics densityRecord = getSoleElement(codeMetricsMap.get(DENSITY_ID.getMetricsCode()));
            final IlluminaTileMetrics clusterRecord = getSoleElement(codeMetricsMap.get(CLUSTER_ID.getMetricsCode()));

            // Snag the phasing data for each read in the read structure. For both types of phasing values, this is the median of all of the individual values seen
            final Collection<TilePhasingValue> tilePhasingValues = getTilePhasingValues(codeMetricsMap, readStructure);

            tiles.add(new Tile(densityRecord.getLaneNumber(), densityRecord.getTileNumber(), densityRecord.getMetricValue(), clusterRecord.getMetricValue(),
                    tilePhasingValues.toArray(new TilePhasingValue[tilePhasingValues.size()])));
        }

        return unmodifiableCollection(tiles);
    }

    /**
     * Pulls out the phasing & prephasing value for the template reads and returns a collection of TilePhasingValues representing these
     */
    private static Collection<TilePhasingValue> getTilePhasingValues(final Map<Integer, Collection<IlluminaTileMetrics>> codeMetricsMap, final ReadStructure readStructure) {
        boolean isFirstRead = true;
        final Collection<TilePhasingValue> tilePhasingValues = new ArrayList<>();
        for (int descriptorIndex = 0; descriptorIndex < readStructure.descriptors.size(); descriptorIndex++) {
            if (readStructure.descriptors.get(descriptorIndex).type == Template) {
                final TileTemplateRead tileTemplateRead = isFirstRead ? FIRST : SECOND;
                // For both phasing & prephasing, pull out the value and create a TilePhasingValue for further processing
                final int phasingCode = getPhasingCode(descriptorIndex, PHASING_BASE);
                final int prePhasingCode = getPhasingCode(descriptorIndex, PREPHASING_BASE);

                if (!(codeMetricsMap.containsKey(phasingCode) && codeMetricsMap.containsKey(prePhasingCode))) {
                    throw new IlluminaParserException("Don't have both phasing and prephasing values for tile");
                }

                tilePhasingValues.add(new TilePhasingValue(tileTemplateRead,
                        getSoleElement(codeMetricsMap.get(phasingCode)).getMetricValue(),
                        getSoleElement(codeMetricsMap.get(prePhasingCode)).getMetricValue()));
                isFirstRead = false;
            }
        }

        return tilePhasingValues;
    }

    /**
     * According to Illumina, for every lane/tile/code combination they will only use the last value. Filter out the previous values
     */
    private static Collection<IlluminaTileMetrics> determineLastValueForLaneTileMetricsCode(final Iterator<IlluminaTileMetrics>
                                                                                                    tileMetricsIterator) {
        final Map<IlluminaLaneTileCode, IlluminaTileMetrics> filteredTileMetrics = new HashMap<>();
        for (final IlluminaTileMetrics illuminaTileMetrics : new IterableAdapter<>(tileMetricsIterator)) {
            filteredTileMetrics.put(illuminaTileMetrics.getLaneTileCode(), illuminaTileMetrics);
        }

        return filteredTileMetrics.values();
    }

    private static String renderMetricLocationKey(final IlluminaTileMetrics metric) {
        return format("%s:%s", metric.getLaneNumber(), metric.getTileNumber());
    }

    // Wrapper around CollectionUtil.Partitioner, purely to de-bulk the actual methods
    private static Map<Integer, Collection<IlluminaTileMetrics>> partitionTileMetricsByCode(final Collection<IlluminaTileMetrics> tileMetrics) {
        return partition(tileMetrics, new CollectionUtil.Partitioner<IlluminaTileMetrics, Integer>() {
            @Override
            public Integer getPartition(final IlluminaTileMetrics metric) {
                return metric.getMetricCode();
            }
        });
    }

    // Wrapper around CollectionUtil.Partitioner, purely to de-bulk the actual methods
    private static Map<String, Collection<IlluminaTileMetrics>> partitionTileMetricsByLocation(final Collection<IlluminaTileMetrics> tileMetrics) {
        return partition(tileMetrics, new CollectionUtil.Partitioner<IlluminaTileMetrics, String>() {
            @Override
            public String getPartition(final IlluminaTileMetrics metric) {
                return renderMetricLocationKey(metric);
            }
        });
    }
}
