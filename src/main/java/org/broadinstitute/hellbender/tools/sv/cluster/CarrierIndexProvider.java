package org.broadinstitute.hellbender.tools.sv.cluster;

/**
 * Implemented by clustering records that carry their non-CNV carrier sample set as a memory-compact,
 * sorted array of sample indices (into a shared {@link CnvSampleCopyState.Dictionary}) instead of
 * carrier {@link htsjdk.variant.variantcontext.Genotype} objects. Non-CNV sample overlap needs only the
 * carrier sample set, so at large sample counts a common variant (100k+ carriers) can retain ~carrier-count
 * genotype objects per pass-1 item; the {@code int[]} form is ~50x smaller and lets
 * {@link SVClusterLinkage#computeSampleOverlap} intersect via a sorted merge.
 *
 * <p>Returns {@code null} for records that carry genotypes instead (single-pass and every other linkage
 * caller), which take the unchanged {@code getCarrierSampleSet}-based path.</p>
 */
public interface CarrierIndexProvider {
    /** Sorted, de-duplicated carrier sample indices, or {@code null} if this record carries genotypes. */
    int[] getCarrierSampleIndices();
}
