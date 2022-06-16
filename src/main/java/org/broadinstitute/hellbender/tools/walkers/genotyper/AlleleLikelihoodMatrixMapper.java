package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.List;

/**
 * Creates {@link org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix} mappers to be used when working with a subset of the original alleles.
 */
class AlleleLikelihoodMatrixMapper<A extends Allele> {

    private final AlleleListPermutation<A> permutation;
    /**
     * Constructs a new mapper given an allele-list permutation.
     * @param permutation the requested permutation.
     *
     * @throws IllegalArgumentException if {@code permutation} is {@code null}.
     */
    AlleleLikelihoodMatrixMapper(final AlleleListPermutation<A> permutation) {
        this.permutation = Utils.nonNull(permutation);
    }

    <EVIDENCE> LikelihoodMatrix<EVIDENCE, A> mapAlleles(final LikelihoodMatrix<EVIDENCE, A> original) {
        if (permutation.isNonPermuted()) {
            return original;
        }

        return new LikelihoodMatrix<EVIDENCE, A>() {

            @Override
            public List<EVIDENCE> evidence() {
                return original.evidence();
            }

            @Override
            public List<A> alleles() {
                return permutation.toList();
            }

            @Override
            public void set(final int alleleIndex, final int evidenceIndex, final double value) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                Utils.validateArg(evidenceIndex >= 0, "readIndex");
                original.set(permutation.fromIndex(alleleIndex), evidenceIndex, value);
            }

            @Override
            public double get(final int alleleIndex, final int evidenceIndex) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                Utils.validateArg(evidenceIndex >= 0, "readIndex");
                return original.get(permutation.fromIndex(alleleIndex), evidenceIndex);
            }

            @Override
            public int indexOfAllele(final Allele allele) {
                Utils.nonNull(allele);
                return permutation.indexOfAllele(allele);
            }

            @Override
            public int indexOfEvidence(final EVIDENCE read) {
                Utils.nonNull(read);
                return original.indexOfEvidence(read);
            }

            @Override
            public int numberOfAlleles() {
                return permutation.toSize();
            }

            @Override
            public int evidenceCount() {
                return original.evidenceCount();
            }

            @Override
            public A getAllele(final int alleleIndex) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                return original.getAllele(permutation.fromIndex(alleleIndex));
            }

            @Override
            public EVIDENCE getEvidence(final int evidenceIndex) {
                Utils.validateArg(evidenceIndex >= 0, "readIndex");
                return original.getEvidence(evidenceIndex);
            }

            @Override
            public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                Utils.nonNull(dest);
                original.copyAlleleLikelihoods(permutation.fromIndex(alleleIndex), dest, offset);
            }
        };
    }
}
