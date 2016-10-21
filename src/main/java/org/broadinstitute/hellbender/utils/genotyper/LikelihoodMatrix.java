package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Likelihood matrix between a set of alleles and reads.
 * @param <A> the allele-type.
 */
public interface LikelihoodMatrix<A extends Allele> extends AlleleList<A> {

    /**
     * List of reads in the matrix sorted by their index therein.
     * @return never {@code null}.
     */
    public List<GATKRead> reads();

    /**
     * List of alleles in the matrix sorted by their index in the collection.
     * @return never {@code null}.
     */
    public List<A> alleles();

    /**
     * Set the likelihood of a read given an allele through their indices.
     *
     * @param alleleIndex the target allele index.
     * @param readIndex the target read index.
     * @param value new likelihood value for the target read give the target allele.
     *
     * @throws IllegalArgumentException if {@code alleleIndex} or {@code readIndex}
     *  are not valid allele and read indices respectively.
     */
    public void set(final int alleleIndex, final int readIndex, final double value);

    /**
     * Returns the likelihood of a read given a haplotype.
     *
     * @param alleleIndex the index of the given haplotype.
     * @param readIndex the index of the target read.
     *
     * @throws IllegalArgumentException if {@code alleleIndex} or {@code readIndex} is not a
     * valid allele or read index respectively.
     *
     * @return the requested likelihood, whatever value was provided using {@link #set(int,int,double) set}
     *    or 0.0 if none was set.
     */
    public double get(final int alleleIndex, final int readIndex);

    /**
     * Queries the index of an allele in the matrix.
     *
     * @param allele the target allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     * @return -1 if such allele does not exist, otherwise its index which 0 or greater.
     */
    @Override
    public int indexOfAllele(final A allele);

    /**
     * Queries the index of a read in the matrix.
     *
     * @param read the target read.
     *
     * @throws IllegalArgumentException if {@code read} is {@code null}.
     *
     * @return -1 if there is not such a read in the matrix, otherwise its index
     *    which is 0 or greater.
     */
    public int indexOfRead(final GATKRead read);

    /**
     * Number of allele in the matrix.
     * @return never negative.
     */
    @Override
    public int numberOfAlleles();

    /**
     * Number of reads in the matrix.
     * @return never negative.
     */
    public int numberOfReads();

    /**
     * Returns the allele given its index.
     *
     * @param alleleIndex the target allele index.
     *
     * @throws IllegalArgumentException if {@code alleleIndex} is not a valid allele index.
     * @return never {@code null}.
     */
    @Override
    public A getAllele(final int alleleIndex);

    /**
     * Returns the allele given its index.
     *
     * @param readIndex the target allele index.
     *
     * @throws IllegalArgumentException if {@code readIndex} is not a valid read index.
     * @return never {@code null}.
     */
    public GATKRead getRead(final int readIndex);


    /**
     * Copies the likelihood of all the reads for a given allele into an array from a particular offset.
     * @param alleleIndex the targeted allele
     * @param dest the destination array.
     * @param offset the copy offset within the destination allele
     */
    public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset);
}
