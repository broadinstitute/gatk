package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Embodies characteristics that describe a lane.  
 * @author mccowan
 */
public final class IlluminaLaneMetrics extends MetricBase {
    /** The number of clusters per unit area on the this lane expressed in units of [cluster / mm^2]. */
    public Double CLUSTER_DENSITY;
    
    /** This lane's number. */
    public Long LANE;
    
    /** This property is not exposed in a field to avoid complications with MetricBase's dependency on reflection. */
    public static String getExtension() {
        return "illumina_lane_metrics";
    }
}
