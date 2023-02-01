package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.List;

/**
 * Some useful functions for creating different kinds of {@link SVClusterEngine}.
 */
public class SVClusterEngineFactory {

    public static SVClusterEngine createCanonical(final SVClusterEngine.CLUSTERING_TYPE type,
                                                  final CanonicalSVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy,
                                                  final CanonicalSVCollapser.AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                                  final SAMSequenceDictionary dictionary,
                                                  final ReferenceSequenceFile reference,
                                                  final boolean enableCNV,
                                                  final ClusteringParameters depthParameters,
                                                  final ClusteringParameters mixedParameters,
                                                  final ClusteringParameters pesrParameters) {
        final CanonicalSVLinkage<SVCallRecord> linkage = new CanonicalSVLinkage<SVCallRecord>(dictionary, enableCNV);
        linkage.setDepthOnlyParams(depthParameters);
        linkage.setMixedParams(mixedParameters);
        linkage.setEvidenceParams(pesrParameters);
        final CanonicalSVCollapser collapser = new CanonicalSVCollapser(reference, altAlleleSummaryStrategy, breakpointSummaryStrategy);
        return new SVClusterEngine(type, collapser::collapse, linkage, dictionary);
    }

    public static SVClusterEngine createCNVDefragmenter(final SAMSequenceDictionary dictionary,
                                                        final CanonicalSVCollapser.AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                                        final ReferenceSequenceFile reference,
                                                        final double paddingFraction,
                                                        final double minSampleOverlap) {
        final SVClusterLinkage<SVCallRecord> linkage = new CNVLinkage(dictionary, paddingFraction, minSampleOverlap);
        final CanonicalSVCollapser collapser = new CanonicalSVCollapser(reference, altAlleleSummaryStrategy, CanonicalSVCollapser.BreakpointSummaryStrategy.MIN_START_MAX_END);
        return new SVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, collapser::collapse, linkage, dictionary);
    }

    public static SVClusterEngine createBinnedCNVDefragmenter(final SAMSequenceDictionary dictionary,
                                                              final CanonicalSVCollapser.AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                                              final ReferenceSequenceFile reference,
                                                              final double paddingFraction,
                                                              final double minSampleOverlap,
                                                              final List<GenomeLoc> coverageIntervals) {
        final SVClusterLinkage<SVCallRecord> linkage = new BinnedCNVLinkage(dictionary, paddingFraction, minSampleOverlap, coverageIntervals);
        final CanonicalSVCollapser collapser = new CanonicalSVCollapser(reference, altAlleleSummaryStrategy, CanonicalSVCollapser.BreakpointSummaryStrategy.MIN_START_MAX_END);
        return new SVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, collapser::collapse, linkage, dictionary);
    }
}
