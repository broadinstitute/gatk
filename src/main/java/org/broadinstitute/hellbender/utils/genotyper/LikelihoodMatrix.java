package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;

import java.util.List;

/**
 * Likelihood matrix between a set of alleles and evidence.
 * @param <A> the allele-type.
 */
public interface LikelihoodMatrix<EVIDENCE,A extends Allele> extends AlleleList<A> {

    /**
     * List of evidence in the matrix sorted by their index therein.
     * @return never {@code null}.
     */
    List<EVIDENCE> evidence();

    /**
     * List of alleles in the matrix sorted by their index in the collection.
     * @return never {@code null}.
     */
    List<A> alleles();

    /**
     * Set the likelihood of a unit of evidence given an allele through their indices.
     *
     * @param alleleIndex the target allele index.
     * @param evidenceIndex the target evidence index.
     * @param value new likelihood value for the target evidence give the target allele.
     *
     * @throws IllegalArgumentException if {@code alleleIndex} or {@code evidenceIndex}
     *  are not valid allele and evidence indices respectively.
     */
    void set(final int alleleIndex, final int evidenceIndex, final double value);

    /**
     * Returns the likelihood of a unit of evidence given a haplotype.
     *
     * @param alleleIndex the index of the given haplotype.
     * @param evidenceIndex the index of the target evidence.
     *
     * @throws IllegalArgumentException if {@code alleleIndex} or {@code evidenceIndex} is not a
     * valid allele or evidence index respectively.
     *
     * @return the requested likelihood, whatever value was provided using {@link #set(int,int,double) set}
     *    or 0.0 if none was set.
     */
    double get(final int alleleIndex, final int evidenceIndex);

    /**
     * Queries the index of an allele in the matrix.
     *
     * @param allele the target allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     * @return -1 if such allele does not exist, otherwise its index which 0 or greater.
     */
    @Override
    int indexOfAllele(final Allele allele);

    /**
     * Queries the index of a unit of evidence in the matrix.
     *
     * @param evidence the target evidence.
     *
     * @throws IllegalArgumentException if {@code evidence} is {@code null}.
     *
     * @return -1 if there is not such a evidence in the matrix, otherwise its index
     *    which is 0 or greater.
     */
    int indexOfEvidence(final EVIDENCE evidence);

    /**
     * Number of allele in the matrix.
     * @return never negative.
     */
    @Override
    int numberOfAlleles();

    /**
     * Count of evidence in the matrix.
     * @return never negative.
     */
    int evidenceCount();

    /**
     * Returns the allele given its index.
     *
     * @param alleleIndex the target allele index.
     *
     * @throws IllegalArgumentException if {@code alleleIndex} is not a valid allele index.
     * @return never {@code null}.
     */
    @Override
    A getAllele(final int alleleIndex);

    /**
     * Returns the allele given its index.
     *
     * @param evidenceIndex the target allele index.
     *
     * @throws IllegalArgumentException if {@code evidenceIndex} is not a valid evidence index.
     * @return never {@code null}.
     */
    EVIDENCE getEvidence(final int evidenceIndex);


    /**
     * Copies the likelihood of all the evidence for a given allele into an array from a particular offset.
     * @param alleleIndex the targeted allele
     * @param dest the destination array.
     * @param offset the copy offset within the destination allele
     */
    void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset);
}
