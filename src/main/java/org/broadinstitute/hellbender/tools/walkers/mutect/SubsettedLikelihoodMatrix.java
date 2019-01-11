package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Fast wrapper for a LikelihoodMatrix that uses only a subset of alleles.  Useful for model comparison of different
 * allele subsets without having to copy the underlying likelihoods.
 * Created by davidben on 1/26/17.
 */
//TODO: consider making the constructor a static method that returns an anonymous class instance in AlleleSubsettingUtils
public class SubsettedLikelihoodMatrix<A extends Allele> implements LikelihoodMatrix<A> {
    private final LikelihoodMatrix<A> matrix;
    private final List<A> alleles;
    private final Int2IntMap newToOldIndexMap;

    public SubsettedLikelihoodMatrix(final LikelihoodMatrix<A> matrix, final List<A> alleles) {
        this.matrix = Utils.nonNull(matrix);
        this.alleles = Utils.nonNull(alleles);
        final int[] newIndices = new IndexRange(0, alleles.size()).mapToInteger(n -> n);
        final int[] oldIndices = alleles.stream().mapToInt(matrix::indexOfAllele).toArray();
        Utils.validateArg(Arrays.stream(oldIndices).noneMatch(n -> n < 0), "All alleles must be found in likelihoods matrix");
        newToOldIndexMap = new Int2IntArrayMap(newIndices, oldIndices);
    }

    public static <A extends Allele> SubsettedLikelihoodMatrix<A> excludingAllele(final LikelihoodMatrix<A> matrix, final Allele excludedAllele) {
        final List<A> alleles = matrix.alleles().stream().filter(a -> !basesMatch(a,excludedAllele)).collect(Collectors.toList());
        Utils.validateArg(alleles.size() == matrix.numberOfAlleles() - 1, "More than one allele excluded.");
        return new SubsettedLikelihoodMatrix<A>(matrix, alleles);
    }

    //TODO: take this hack out
    public static boolean basesMatch(final Allele a, final Allele b) { return a.getBases() == b.getBases() || Arrays.equals(a.getBases(), b.getBases()); }

    @Override
    public List<GATKRead> reads() { return matrix.reads(); }

    @Override
    public List<A> alleles() { return alleles; }

    @Override
    public void set(final int alleleIndex, final int readIndex, final double value) {
        throw new UnsupportedOperationException("Subsetted likelihood matrices are immutable.");
    }

    @Override
    public double get(final int alleleIndex, final int readIndex) { return matrix.get(newToOldIndexMap.get(alleleIndex), readIndex); }

    @Override
    public int indexOfAllele(final A allele) { return alleles.indexOf(allele); }

    @Override
    public int indexOfRead(final GATKRead read) { return matrix.indexOfRead(read); }

    @Override
    public int numberOfAlleles() { return alleles.size(); }

    @Override
    public int numberOfReads() { return matrix.numberOfReads(); }

    @Override
    public A getAllele(final int alleleIndex) { return alleles.get(alleleIndex); }

    @Override
    public GATKRead getRead(final int readIndex) { return matrix.getRead(readIndex); }

    @Override
    public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
        throw new UnsupportedOperationException("Subsetted likelihood matrices are not meant to be copied.");
    }
}
