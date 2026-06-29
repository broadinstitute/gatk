package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class SVClusterLinkage<T extends SVLocatable> {

    /**
     * Returns whether two given items cluster.
     * @param a first item
     * @param b second item
     */
    public abstract LinkageResult areClusterable(final T a, final T b);

    /**
     * Returns the maximum feasible starting position of any other item with the given item. That is, given item A and
     * X = getMaxClusterableStartingPosition(A), then for any item B on the current contig,
     * Y = start(B) > X => clusterTogether(A, B) == false. Note that this is an upper-bound, but tighter estimates
     * can greatly improve performance.
     * @param item item in question
     * @return max feasible clusterable start coordinate on the current contig
     */
    public abstract int getMaxClusterableStartingPosition(final T item);

    /**
     * Compute max feasible starting position of any other item for all items in the given collection. Note the items
     * must all have the same starting contig.
     */
    public int getMaxClusterableStartingPosition(final Collection<T> items) {
        final List<String> contigA = items.stream().map(T::getContigA).distinct().collect(Collectors.toList());
        if (contigA.size() > 1) {
            throw new IllegalArgumentException("Items start on multiple contigs");
        }
        return items.stream().mapToInt(item -> getMaxClusterableStartingPosition(item)).max().getAsInt();
    }

    /**
     * Returns number of overlapping items
     */
    protected static double getSampleSetOverlap(final Collection<String> a, final Set<String> b) {
        final double denom = Math.max(a.size(), b.size());
        if (denom == 0) {
            return 1;
        }
        return a.stream().filter(b::contains).count() / denom;
    }

    /**
     * Returns true if the overlap is null or exceeds threshold.
     */
    protected static boolean testSampleOverlap(final Double sampleOverlap, final double threshold) {
        return sampleOverlap == null || sampleOverlap >= threshold;
    }


    /**
     * Returns fractional carrier sample overlap in the two records. For CNVs, returns fraction of copy number states that match.
     */
    protected static Double computeSampleOverlap(final SVCallRecord a, final SVCallRecord b) {
        if (a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV || b.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            // CNV sample overlap. Memory-compact pass-1 CNV items expose their per-sample copy state as
            // primitive arrays (CnvSampleCopyState) instead of retaining a Genotype per sample; use those
            // when present, falling back to the genotype-based computation otherwise. Both paths produce
            // the identical fraction (locked by SVClusterLinkageCompactCopyStateTest).
            final CnvSampleCopyState compactA = a instanceof CnvCopyStateProvider ? ((CnvCopyStateProvider) a).getCnvCopyState() : null;
            final CnvSampleCopyState compactB = b instanceof CnvCopyStateProvider ? ((CnvCopyStateProvider) b).getCnvCopyState() : null;
            if (compactA != null && compactB != null) {
                return computeCompactCnvOverlap(compactA, compactB);
            } else if (compactA != null) {
                return computeMixedCnvOverlap(compactA, b.getGenotypes());
            } else if (compactB != null) {
                return computeMixedCnvOverlap(compactB, a.getGenotypes());
            }
            final GenotypesContext genotypesA = a.getGenotypes();
            final GenotypesContext genotypesB = b.getGenotypes();
            final Set<String> samples = new HashSet<>(SVUtils.hashMapCapacity(genotypesA.size() + genotypesB.size()));
            samples.addAll(genotypesA.getSampleNames());
            samples.addAll(genotypesB.getSampleNames());
            if (samples.isEmpty()) {
                return null;
            }
            int numMatches = 0;
            for (final String sample : samples) {
                final Genotype genotypeA = genotypesA.get(sample);
                final Genotype genotypeB = genotypesB.get(sample);
                // If one sample doesn't exist in the other set, assume reference copy state
                final int cnA = getCopyState(genotypeA, genotypeB);
                final int cnB = getCopyState(genotypeB, genotypeA);
                if (cnA == cnB) {
                    numMatches++;
                }
            }
            final int numSamples = samples.size();
            return (numMatches / (double) numSamples);
        } else {
            // Non-CNV
            final Set<String> samplesA = a.getCarrierSampleSet();
            final Set<String> samplesB = b.getCarrierSampleSet();
            return getSampleSetOverlap(samplesA, samplesB);
        }
    }

    /**
     * CNV sample overlap when both records carry compact copy state (share one {@link CnvSampleCopyState.Dictionary}).
     * Replicates {@link #computeSampleOverlap}'s CNV loop over the union of genotyped samples: present side
     * uses {@code copyState}, missing side uses the other record's {@code expectedCopyState} (ECN).
     */
    private static Double computeCompactCnvOverlap(final CnvSampleCopyState a, final CnvSampleCopyState b) {
        final int n = a.getDictionary().size();
        int numMatches = 0;
        int numSamples = 0;
        for (int i = 0; i < n; i++) {
            final boolean pa = a.isPresent(i);
            final boolean pb = b.isPresent(i);
            if (!pa && !pb) {
                continue; // sample genotyped in neither record -> not in the union
            }
            numSamples++;
            final int cnA = pa ? a.copyStateAt(i) : b.expectedCopyStateAt(i);
            final int cnB = pb ? b.copyStateAt(i) : a.expectedCopyStateAt(i);
            if (cnA == cnB) {
                numMatches++;
            }
        }
        return numSamples == 0 ? null : numMatches / (double) numSamples;
    }

    /**
     * CNV sample overlap when one record is compact ({@code a}) and the other carries genotypes ({@code bGenotypes},
     * e.g. a carrier-stripped non-CNV record). Union = a's genotyped samples plus any sample present only in
     * {@code bGenotypes}; per-sample compare matches {@link #computeSampleOverlap}'s CNV loop exactly.
     */
    private static Double computeMixedCnvOverlap(final CnvSampleCopyState a, final GenotypesContext bGenotypes) {
        final CnvSampleCopyState.Dictionary dict = a.getDictionary();
        int numMatches = 0;
        int numSamples = 0;
        // a's genotyped samples (present side a); b may or may not have each.
        for (int i = a.nextPresent(0); i >= 0; i = a.nextPresent(i + 1)) {
            final Genotype gB = bGenotypes.get(dict.sampleAt(i));
            final int cnA = a.copyStateAt(i);
            final int cnB = gB != null ? CnvSampleCopyState.copyStateOf(gB) : a.expectedCopyStateAt(i);
            if (cnA == cnB) {
                numMatches++;
            }
            numSamples++;
        }
        // Samples present only in b (not genotyped in a): a is the missing side -> a uses b's ECN.
        for (final Genotype gB : bGenotypes) {
            final int idx = dict.indexOf(gB.getSampleName());
            if (idx >= 0 && a.isPresent(idx)) {
                continue; // already counted in the a-present loop
            }
            final int cnA = VariantContextGetters.getAttributeAsInt(gB, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1);
            final int cnB = CnvSampleCopyState.copyStateOf(gB);
            if (cnA == cnB) {
                numMatches++;
            }
            numSamples++;
        }
        return numSamples == 0 ? null : numMatches / (double) numSamples;
    }

    /**
     * Tries to get the best copy state from the genotype. If the genotype is null, uses ploidy from a "backup"
     * genotype as the default. If we have no clue, just return -1 as a null default.
     */
    private static int getCopyState(final Genotype genotype, final Genotype matchedSampleGenotype) {
        if (genotype == null) {
            if (matchedSampleGenotype != null) {
                return VariantContextGetters.getAttributeAsInt(matchedSampleGenotype, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1);
            } else {
                throw new IllegalArgumentException("Both genotypes are null");
            }
        } else {
            return VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.COPY_NUMBER_FORMAT,
                    VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, -1));
        }
    }

    /**
     * Used for storing the result of a clustering check between two records and any additional metadata
     */
    public static class LinkageResult {
        private final boolean result;
        public LinkageResult(final boolean result) { this.result = result; }
        public boolean getResult() { return result; }
    }
}
