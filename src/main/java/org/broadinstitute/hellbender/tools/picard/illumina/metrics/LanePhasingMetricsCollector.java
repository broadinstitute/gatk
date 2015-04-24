package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.util.CollectionUtil;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.Tile;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.TileTemplateRead;
import org.broadinstitute.hellbender.utils.Median;

import java.util.*;

import static java.util.Collections.unmodifiableMap;

/**
 * Helper class used to transform tile data for a lane into a collection of IlluminaPhasingMetrics
 */
public final class LanePhasingMetricsCollector {
    private final Map<TileTemplateRead, Float> medianPhasingMap;
    private final Map<TileTemplateRead, Float> medianPrePhasingMap;

    /**
     * Constructor takes a lane's collection of Tiles and calculates the median phasing/prephasing for the
     * first and second (if available) reads
     */
    public LanePhasingMetricsCollector(final Collection<Tile> laneTiles) {
        final Map<TileTemplateRead, Float> medianPhasingMap = new TreeMap<TileTemplateRead, Float>();
        final Map<TileTemplateRead, Float> medianPrePhasingMap = new TreeMap<TileTemplateRead, Float>();

        final CollectionUtil.MultiMap<TileTemplateRead, Float> phasingValues = new CollectionUtil.MultiMap<TileTemplateRead, Float>();
        final CollectionUtil.MultiMap<TileTemplateRead, Float> prePhasingValues = new CollectionUtil.MultiMap<TileTemplateRead, Float>();

        // Collect the phasing/prephasing values from all of the tiles, sorted by template read #
        for (final Tile tile : laneTiles) {
            for (final TileTemplateRead tileTemplateRead : tile.getPhasingMap().keySet()) {
                phasingValues.append(tileTemplateRead, tile.getPhasingMap().get(tileTemplateRead));
                prePhasingValues.append(tileTemplateRead, tile.getPrePhasingMap().get(tileTemplateRead));
            }
        }

        // Calculate the medians for the collected data
        for (final TileTemplateRead tileTemplateRead : phasingValues.keySet()) {
            medianPhasingMap.put(tileTemplateRead, medianPercentage(phasingValues.get(tileTemplateRead)));
            medianPrePhasingMap.put(tileTemplateRead, medianPercentage(prePhasingValues.get(tileTemplateRead)));
        }

        this.medianPhasingMap = unmodifiableMap(medianPhasingMap);
        this.medianPrePhasingMap = unmodifiableMap(medianPrePhasingMap);
    }

    public Map<TileTemplateRead, Float> getMedianPhasingMap() {
        return medianPhasingMap;
    }

    public Map<TileTemplateRead, Float> getMedianPrePhasingMap() {
        return medianPrePhasingMap;
    }

    private static float medianPercentage(final Collection<Float> phaseValues) {
        final Median<Double> median = new Median<>();
        for (final Float f : phaseValues) median.add((double) f);
        return median.getMedian().floatValue() * 100;
    }
}
