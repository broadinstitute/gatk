package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

/**
 * Arguments for use with {@link SVClusterEngine}.
 */
public class SVClusterEngineArgumentsCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    protected static final String BASE_INTERVAL_OVERLAP_FRACTION_NAME = "-interval-overlap";
    protected static final String BASE_SIZE_SIMILARITY_FRACTION_NAME = "-size-similarity";
    protected static final String BASE_BREAKEND_WINDOW_NAME = "-breakend-window";
    protected static final String BASE_SAMPLE_OVERLAP_FRACTION_NAME = "-sample-overlap";

    protected static final String BASE_INTERVAL_OVERLAP_FRACTION_DOC = " interval reciprocal overlap fraction";
    protected static final String BASE_SIZE_SIMILARITY_DOC = " size similarity";
    protected static final String BASE_BREAKEND_WINDOW_DOC = " window size for breakend proximity";
    protected static final String BASE_SAMPLE_OVERLAP_FRACTION_DOC = " shared sample overlap fraction";

    public static final String DEPTH_INTERVAL_OVERLAP_FRACTION_NAME = "depth" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;
    public static final String MIXED_INTERVAL_OVERLAP_FRACTION_NAME = "mixed" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;
    public static final String PESR_INTERVAL_OVERLAP_FRACTION_NAME = "pesr" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;

    public static final String DEPTH_SIZE_SIMILARITY_NAME = "depth" + BASE_SIZE_SIMILARITY_FRACTION_NAME;
    public static final String MIXED_SIZE_SIMILARITY_NAME = "mixed" + BASE_SIZE_SIMILARITY_FRACTION_NAME;
    public static final String PESR_SIZE_SIMILARITY_NAME = "pesr" + BASE_SIZE_SIMILARITY_FRACTION_NAME;

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
    public double depthOverlapFraction = CanonicalSVLinkage.DEFAULT_RECIPROCAL_OVERLAP_DEPTH_ONLY;

    /**
     * Minimum interval reciprocal overlap fraction to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="PESR/Depth" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double mixedOverlapFraction = CanonicalSVLinkage.DEFAULT_RECIPROCAL_OVERLAP_MIXED;

    /**
     * Minimum interval reciprocal overlap fraction to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double pesrOverlapFraction = CanonicalSVLinkage.DEFAULT_RECIPROCAL_OVERLAP_PESR;

    /**
     * Minimum size similarity to cluster depth-only/depth-only variant pairs.
     */
    @Argument(fullName = DEPTH_SIZE_SIMILARITY_NAME,
            doc="Depth/Depth" + BASE_SIZE_SIMILARITY_DOC, optional=true, minValue = 0, maxValue = 1)
    public double depthSizeSimilarity = CanonicalSVLinkage.DEFAULT_SIZE_SIMILARITY_DEPTH_ONLY;

    /**
     * Minimum size similarity to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_SIZE_SIMILARITY_NAME,
            doc="Depth/PESR" + BASE_SIZE_SIMILARITY_DOC, optional=true, minValue = 0, maxValue = 1)
    public double mixedSizeSimilarity = CanonicalSVLinkage.DEFAULT_SIZE_SIMILARITY_MIXED;

    /**
     * Minimum size similarity to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_SIZE_SIMILARITY_NAME,
            doc="PESR/PESR" + BASE_SIZE_SIMILARITY_DOC, optional=true, minValue = 0, maxValue = 1)
    public double pesrSizeSimilarity = CanonicalSVLinkage.DEFAULT_SIZE_SIMILARITY_PESR;

    /**
     * Maximum allowed distance between endpoints (in bp) to cluster depth-only/depth-only variant pairs.
     */
    @Argument(fullName = DEPTH_BREAKEND_WINDOW_NAME,
            doc="Depth/Depth" + BASE_BREAKEND_WINDOW_DOC, optional=true, minValue = 0)
    public int depthBreakendWindow = CanonicalSVLinkage.DEFAULT_WINDOW_DEPTH_ONLY;

    /**
     * Maximum allowed distance between endpoints (in bp) to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_BREAKEND_WINDOW_NAME,
            doc="Depth/PESR" + BASE_BREAKEND_WINDOW_DOC, optional=true, minValue = 0)
    public int mixedBreakendWindow = CanonicalSVLinkage.DEFAULT_WINDOW_MIXED;

    /**
     * Maximum allowed distance between endpoints (in bp) to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_BREAKEND_WINDOW_NAME,
            doc="PESR/PESR" + BASE_BREAKEND_WINDOW_DOC, optional=true, minValue = 0)
    public int pesrBreakendWindow = CanonicalSVLinkage.DEFAULT_WINDOW_PESR;

    /**
     * Minimum carrier sample reciprocal overlap fraction to cluster depth-only/depth-only variant pairs.
     */
    @Argument(fullName = DEPTH_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="Depth/Depth" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double depthSampleOverlapFraction = CanonicalSVLinkage.DEFAULT_SAMPLE_OVERLAP_DEPTH_ONLY;

    /**
     * Minimum carrier sample reciprocal overlap fraction to cluster depth-only/PESR variant pairs.
     */
    @Argument(fullName = MIXED_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="Depth/PESR" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double mixedSampleOverlapFraction = CanonicalSVLinkage.DEFAULT_SAMPLE_OVERLAP_MIXED;

    /**
     * Minimum carrier sample reciprocal overlap fraction to cluster PESR/PESR variant pairs.
     */
    @Argument(fullName = PESR_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true, minValue = 0, maxValue = 1)
    public double pesrSampleOverlapFraction = CanonicalSVLinkage.DEFAULT_SAMPLE_OVERLAP_PESR;

    public final ClusteringParameters getDepthParameters() {
        return ClusteringParameters.createDepthParameters(depthOverlapFraction, depthSizeSimilarity, depthBreakendWindow, depthSampleOverlapFraction);
    }

    public final ClusteringParameters getMixedParameters() {
        return ClusteringParameters.createMixedParameters(mixedOverlapFraction, mixedSizeSimilarity, mixedBreakendWindow, mixedSampleOverlapFraction);
    }

    public final ClusteringParameters getPESRParameters() {
        return ClusteringParameters.createPesrParameters(pesrOverlapFraction, pesrSizeSimilarity, pesrBreakendWindow, pesrSampleOverlapFraction);
    }
}
