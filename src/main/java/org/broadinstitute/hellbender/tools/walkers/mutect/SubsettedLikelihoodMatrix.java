package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Fast wrapper for a LikelihoodMatrix that uses only a subset of alleles.  Useful for model comparison of different
 * allele subsets without having to copy the underlying likelihoods.
 * Created by davidben on 1/26/17.
 */
//TODO: consider making the constructor a static method that returns an anonymous class instance in AlleleSubsettingUtils
public class SubsettedLikelihoodMatrix<EVIDENCE extends Locatable, A extends Allele> implements LikelihoodMatrix<EVIDENCE, A> {
    private final LikelihoodMatrix<EVIDENCE, A> matrix;
    private final List<A> alleles;
    private final Int2IntMap newToOldIndexMap;

    public SubsettedLikelihoodMatrix(final LikelihoodMatrix<EVIDENCE, A> matrix, final List<A> alleles) {
        this.matrix = Utils.nonNull(matrix);
        this.alleles = Utils.nonNull(alleles);
        final int[] newIndices = new IndexRange(0, alleles.size()).mapToInteger(n -> n);
        final int[] oldIndices = alleles.stream().mapToInt(matrix::indexOfAllele).toArray();
        Utils.validateArg(Arrays.stream(oldIndices).noneMatch(n -> n < 0), "All alleles must be found in likelihoods matrix");
        newToOldIndexMap = new Int2IntArrayMap(newIndices, oldIndices);
    }

    public static <EVIDENCE extends Locatable, A extends Allele> SubsettedLikelihoodMatrix<EVIDENCE,A> excludingAllele(final LikelihoodMatrix<EVIDENCE,A> matrix, final Allele excludedAllele) {
        final List<A> alleles = matrix.alleles().stream().filter(a -> !basesMatch(a,excludedAllele)).collect(Collectors.toList());
        Utils.validateArg(alleles.size() == matrix.numberOfAlleles() - 1, "More than one allele excluded.");
        return new SubsettedLikelihoodMatrix<EVIDENCE,A>(matrix, alleles);
    }

    //TODO: take this hack out
    public static boolean basesMatch(final Allele a, final Allele b) { return a.getBases() == b.getBases() || Arrays.equals(a.getBases(), b.getBases()); }

    @Override
    public List<EVIDENCE> evidence() { return matrix.evidence(); }

    @Override
    public List<A> alleles() { return alleles; }

    @Override
    public void set(final int alleleIndex, final int evidenceIndex, final double value) {
        throw new UnsupportedOperationException("Subsetted likelihood matrices are immutable.");
    }

    @Override
    public double get(final int alleleIndex, final int evidenceIndex) { return matrix.get(newToOldIndexMap.get(alleleIndex), evidenceIndex); }

    @Override
    public int indexOfAllele(final Allele allele) { return alleles.indexOf(allele); }

    @Override
    public int indexOfEvidence(final EVIDENCE evidence) { return matrix.indexOfEvidence(evidence); }

    @Override
    public int numberOfAlleles() { return alleles.size(); }

    @Override
    public int evidenceCount() { return matrix.evidenceCount(); }

    @Override
    public A getAllele(final int alleleIndex) { return alleles.get(alleleIndex); }

    @Override
    public EVIDENCE getEvidence(final int evidenceIndex) { return matrix.getEvidence(evidenceIndex); }

    @Override
    public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
        throw new UnsupportedOperationException("Subsetted likelihood matrices are not meant to be copied.");
    }
}
