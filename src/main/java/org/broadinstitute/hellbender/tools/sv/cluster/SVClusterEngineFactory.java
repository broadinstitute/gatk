package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.List;

/**
 * Some useful functions for creating different kinds of {@link SVClusterEngine}.
 */
public class SVClusterEngineFactory {

    public static <T extends SVCallRecord> CanonicalSVClusterEngine<T> createCanonical(final SVClusterEngine.CLUSTERING_TYPE type,
                                                                              final SAMSequenceDictionary dictionary,
                                                                              final boolean enableCNV,
                                                                              final ClusteringParameters depthParameters,
                                                                              final ClusteringParameters mixedParameters,
                                                                              final ClusteringParameters pesrParameters) {
        final CanonicalSVLinkage<T> linkage = new CanonicalSVLinkage<>(dictionary, enableCNV);
        linkage.setDepthOnlyParams(depthParameters);
        linkage.setMixedParams(mixedParameters);
        linkage.setEvidenceParams(pesrParameters);
        return new CanonicalSVClusterEngine<T>(type, linkage, dictionary);
    }

    public static <T extends SVCallRecord> CanonicalSVClusterEngine<T> createCNVDefragmenter(final SAMSequenceDictionary dictionary,
                                                                                    final double paddingFraction,
                                                                                    final double minSampleOverlap) {
        final SVClusterLinkage<T> linkage = new CNVLinkage<>(dictionary, paddingFraction, minSampleOverlap);
        return new CanonicalSVClusterEngine<T>(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, linkage, dictionary);
    }

    public static <T extends SVCallRecord> CanonicalSVClusterEngine<T> createBinnedCNVDefragmenter(final SAMSequenceDictionary dictionary,
                                                                                          final double paddingFraction,
                                                                                          final double minSampleOverlap,
                                                                                          final List<GenomeLoc> coverageIntervals) {
        final SVClusterLinkage<T> linkage = new BinnedCNVLinkage<>(dictionary, paddingFraction, minSampleOverlap, coverageIntervals);
        return new CanonicalSVClusterEngine<>(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, linkage, dictionary);
    }
}
