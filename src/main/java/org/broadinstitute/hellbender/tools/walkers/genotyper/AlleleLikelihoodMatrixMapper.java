package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.function.UnaryOperator;

/**
 * Creates {@link org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix} mappers to be used when working with a subset of the original alleles.
 */
@FunctionalInterface
public interface AlleleLikelihoodMatrixMapper<A extends Allele> extends UnaryOperator<LikelihoodMatrix<A>>{

    /**
     * Instantiates a new mapper given an allele-list permutation.
     * @param permutation the requested permutation.
     * @param <A> the allele type.
     *
     * @throws IllegalArgumentException if {@code permutation} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static <A extends Allele> AlleleLikelihoodMatrixMapper<A> newInstance(final AlleleListPermutation<A> permutation) {
        Utils.nonNull(permutation, "the permutation must not be null");
        if (permutation.isNonPermuted()) {
            return read -> read;
        }
        return original -> new LikelihoodMatrix<A>() {

            @Override
            public List<GATKRead> reads() {
                return original.reads();
            }

            @Override
            public List<A> alleles() {
                return permutation.toList();
            }

            @Override
            public void set(final int alleleIndex, final int readIndex, final double value) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                Utils.validateArg(readIndex >= 0, "readIndex");
                original.set(permutation.fromIndex(alleleIndex), readIndex, value);
            }

            @Override
            public double get(final int alleleIndex, final int readIndex) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                Utils.validateArg(readIndex >= 0, "readIndex");
                return original.get(permutation.fromIndex(alleleIndex),readIndex);
            }

            @Override
            public int indexOfAllele(final A allele) {
                Utils.nonNull(allele);
                return permutation.indexOfAllele(allele);
            }

            @Override
            public int indexOfRead(final GATKRead read) {
                Utils.nonNull(read);
                return original.indexOfRead(read);
            }

            @Override
            public int numberOfAlleles() {
                return permutation.toSize();
            }

            @Override
            public int numberOfReads() {
                return original.numberOfReads();
            }

            @Override
            public A getAllele(final int alleleIndex) {
                Utils.validateArg(alleleIndex >= 0, "alleleIndex");
                return original.getAllele(permutation.fromIndex(alleleIndex));
            }

            @Override
            public GATKRead getRead(final int readIndex) {
                Utils.validateArg(readIndex >= 0, "readIndex");
                return original.getRead(readIndex);
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
