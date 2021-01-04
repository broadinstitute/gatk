package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

/**
 * Arguments for use with {@link SVClusterEngine}.
 */
public class SVClusterEngineArgumentsCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    protected static final String BASE_INTERVAL_OVERLAP_FRACTION_NAME = "-overlap-fraction";
    protected static final String BASE_BREAKEND_WINDOW_NAME = "-breakend-window";
    protected static final String BASE_SAMPLE_OVERLAP_FRACTION_NAME = "-sample-overlap";

    protected static final String BASE_INTERVAL_OVERLAP_FRACTION_DOC = " interval reciprocal overlap fraction";
    protected static final String BASE_BREAKEND_WINDOW_DOC = " window size for breakend proximity";
    protected static final String BASE_SAMPLE_OVERLAP_FRACTION_DOC = " shared sample overlap fraction";

    public static final String DEPTH_INTERVAL_OVERLAP_FRACTION_NAME = "depth" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;
    public static final String MIXED_INTERVAL_OVERLAP_FRACTION_NAME = "mixed" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;
    public static final String PESR_INTERVAL_OVERLAP_FRACTION_NAME = "pesr" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;

    public static final String DEPTH_BREAKEND_WINDOW_NAME = "depth" + BASE_BREAKEND_WINDOW_NAME;
    public static final String MIXED_BREAKEND_WINDOW_NAME = "mixed" + BASE_BREAKEND_WINDOW_NAME;
    public static final String PESR_BREAKEND_WINDOW_NAME = "pesr" + BASE_BREAKEND_WINDOW_NAME;

    public static final String DEPTH_SAMPLE_OVERLAP_FRACTION_NAME = "depth" + BASE_SAMPLE_OVERLAP_FRACTION_NAME;
    public static final String MIXED_SAMPLE_OVERLAP_FRACTION_NAME = "mixed" + BASE_SAMPLE_OVERLAP_FRACTION_NAME;
    public static final String PESR_SAMPLE_OVERLAP_FRACTION_NAME = "pesr" + BASE_SAMPLE_OVERLAP_FRACTION_NAME;

    /**
     * Minimum interval reciprocal overlap fraction to cluster depth-only/depth-only variant pairs.
     */
    @Argument(fullName = DEPTH_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="Depth/Depth" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double depthOverlapFraction = SVClusterEngine.DEFAULT_DEPTH_ONLY_PARAMS.getReciprocalOverlap();

    /**
     * Minimum interval reciprocal overlap fraction to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="PESR/Depth" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double mixedOverlapFraction = SVClusterEngine.DEFAULT_MIXED_PARAMS.getReciprocalOverlap();

    /**
     * Minimum interval reciprocal overlap fraction to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double pesrOverlapFraction = SVClusterEngine.DEFAULT_EVIDENCE_PARAMS.getReciprocalOverlap();

    /**
     * Maximum allowed distance between endpoints (in bp) to cluster depth-only/depth-only variant pairs.
     */
    @Argument(fullName = DEPTH_BREAKEND_WINDOW_NAME,
            doc="Depth/Depth" + BASE_BREAKEND_WINDOW_DOC, optional=true, minValue = 0)
    public int depthBreakendWindow = SVClusterEngine.DEFAULT_DEPTH_ONLY_PARAMS.getWindow();

    /**
     * Maximum allowed distance between endpoints (in bp) to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_BREAKEND_WINDOW_NAME,
            doc="Depth/PESR" + BASE_BREAKEND_WINDOW_DOC, optional=true, minValue = 0)
    public int mixedBreakendWindow = SVClusterEngine.DEFAULT_MIXED_PARAMS.getWindow();

    /**
     * Maximum allowed distance between endpoints (in bp) to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_BREAKEND_WINDOW_NAME,
            doc="PESR/PESR" + BASE_BREAKEND_WINDOW_DOC, optional=true, minValue = 0)
    public int pesrBreakendWindow = SVClusterEngine.DEFAULT_EVIDENCE_PARAMS.getWindow();

    /**
     * Minimum carrier sample reciprocal overlap fraction to cluster depth-only/depth-only variant pairs.
     */
    @Argument(fullName = DEPTH_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="Depth/Depth" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double depthSampleOverlapFraction = SVClusterEngine.DEFAULT_DEPTH_ONLY_PARAMS.getSampleOverlap();

    /**
     * Minimum carrier sample reciprocal overlap fraction to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="Depth/PESR" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double mixedSampleOverlapFraction = SVClusterEngine.DEFAULT_MIXED_PARAMS.getSampleOverlap();

    /**
     * Minimum carrier sample reciprocal overlap fraction to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double pesrSampleOverlapFraction = SVClusterEngine.DEFAULT_EVIDENCE_PARAMS.getSampleOverlap();

    public final SVClusterEngine.DepthClusteringParameters getDepthParameters() {
        return new SVClusterEngine.DepthClusteringParameters(depthOverlapFraction, depthBreakendWindow, depthSampleOverlapFraction);
    }

    public final SVClusterEngine.MixedClusteringParameters getMixedParameters() {
        return new SVClusterEngine.MixedClusteringParameters(mixedOverlapFraction, mixedBreakendWindow, mixedSampleOverlapFraction);
    }

    public final SVClusterEngine.EvidenceClusteringParameters getPESRParameters() {
        return new SVClusterEngine.EvidenceClusteringParameters(pesrOverlapFraction, pesrBreakendWindow, pesrSampleOverlapFraction);
    }
}
