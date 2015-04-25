package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.Tile;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.TileTemplateRead;

import java.util.*;

/**
 * Metrics for Illumina Basecalling that stores median phasing and prephasing percentages on a per-template-read, per-lane basis.
 * For each lane/template read # (i.e. FIRST, SECOND) combination we will store the median values of both the phasing and prephasing
 * values for every tile in that lane/template read pair.
 *
 * @author jgentry
 */
public class IlluminaPhasingMetrics extends MetricBase {
    public long LANE;
    public String TYPE_NAME;
    public double PHASING_APPLIED;
    public double PREPHASING_APPLIED;

    /**
     * Calculate the median phasing & prephasing values for a lane's tiles and create the appropriate IlluminaPhasingMetrics for them
     */
    public static Collection<IlluminaPhasingMetrics> getPhasingMetricsForTiles(final long lane, final Collection<Tile> tilesForLane) {
        final LanePhasingMetricsCollector lanePhasingMetricsCollector = new LanePhasingMetricsCollector(tilesForLane);
        final Collection<IlluminaPhasingMetrics> phasingMetrics = new ArrayList<>();
        for (final TileTemplateRead tileTemplateRead : lanePhasingMetricsCollector.getMedianPhasingMap().keySet()) {
            final IlluminaPhasingMetrics phasingMetric = new IlluminaPhasingMetrics();
            phasingMetric.LANE = lane;
            phasingMetric.TYPE_NAME = tileTemplateRead.toString();
            phasingMetric.PHASING_APPLIED = lanePhasingMetricsCollector.getMedianPhasingMap().get(tileTemplateRead);
            phasingMetric.PREPHASING_APPLIED = lanePhasingMetricsCollector.getMedianPrePhasingMap().get(tileTemplateRead);
            phasingMetrics.add(phasingMetric);
        }

        return phasingMetrics;
    }

    /**
     * This property is not exposed in a field to avoid complications with MetricBase's dependency on reflection.
     */
    public static String getExtension() {
        return "illumina_phasing_metrics";
    }
}
