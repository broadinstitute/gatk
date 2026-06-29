package org.broadinstitute.hellbender.tools.sv.cluster;

/**
 * Implemented by clustering records that carry a memory-compact {@link CnvSampleCopyState} instead of
 * per-sample CNV genotypes. {@link SVClusterLinkage#computeSampleOverlap} uses the compact state (when
 * non-null) to compute CNV sample overlap via primitive array scans rather than iterating genotype
 * objects, bounding pass-1 memory at large sample counts. Records that do not provide compact state
 * (single-pass, and every other tool that shares the linkage) return {@code null} and take the
 * unchanged genotype-based path.
 */
public interface CnvCopyStateProvider {
    /** Compact per-sample CNV copy state, or {@code null} if this record carries genotypes instead. */
    CnvSampleCopyState getCnvCopyState();
}
