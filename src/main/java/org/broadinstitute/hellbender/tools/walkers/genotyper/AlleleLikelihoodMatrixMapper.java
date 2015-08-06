package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.function.UnaryOperator;

/**
 * Creates {@link org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods.Matrix} mappers to be used when working with a subset of the original alleles.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
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
        return new AlleleLikelihoodMatrixMapper<A>() {
            @Override
            public LikelihoodMatrix<A> apply(final LikelihoodMatrix<A> original) {
                return new LikelihoodMatrix<A>() {

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
                        original.set(permutation.fromIndex(alleleIndex),readIndex,value);
                    }

                    @Override
                    public double get(final int alleleIndex, final int readIndex) {
                        return original.get(permutation.fromIndex(alleleIndex),readIndex);
                    }

                    @Override
                    public int indexOfAllele(final A allele) {
                        return permutation.indexOfAllele(allele);
                    }

                    @Override
                    public int indexOfRead(final GATKRead read) {
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
                        return original.getAllele(permutation.fromIndex(alleleIndex));
                    }

                    @Override
                    public GATKRead getRead(final int readIndex) {
                        return original.getRead(readIndex);
                    }

                    @Override
                    public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
                        original.copyAlleleLikelihoods(permutation.fromIndex(alleleIndex),dest,offset);
                    }
                };
            }
        };
    }

}
