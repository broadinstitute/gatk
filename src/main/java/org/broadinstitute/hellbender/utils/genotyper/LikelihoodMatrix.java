package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.AbstractRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.ojalgo.access.Access2D;
import org.ojalgo.commons.math3.linear.Access2DWrapper;

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
    public List<EVIDENCE> evidence();

    /**
     * List of alleles in the matrix sorted by their index in the collection.
     * @return never {@code null}.
     */
    public List<A> alleles();

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
    public void set(final int alleleIndex, final int evidenceIndex, final double value);

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
    public double get(final int alleleIndex, final int evidenceIndex);

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
     * Queries the index of a unit of evidence in the matrix.
     *
     * @param evidence the target evidence.
     *
     * @throws IllegalArgumentException if {@code evidence} is {@code null}.
     *
     * @return -1 if there is not such a evidence in the matrix, otherwise its index
     *    which is 0 or greater.
     */
    public int indexOfEvidence(final EVIDENCE evidence);

    /**
     * Number of allele in the matrix.
     * @return never negative.
     */
    @Override
    public int numberOfAlleles();

    /**
     * Count of evidence in the matrix.
     * @return never negative.
     */
    public int evidenceCount();

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
     * @param evidenceIndex the target allele index.
     *
     * @throws IllegalArgumentException if {@code evidenceIndex} is not a valid evidence index.
     * @return never {@code null}.
     */
    public EVIDENCE getEvidence(final int evidenceIndex);


    /**
     * Copies the likelihood of all the evidence for a given allele into an array from a particular offset.
     * @param alleleIndex the targeted allele
     * @param dest the destination array.
     * @param offset the copy offset within the destination allele
     */
    public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset);

    /**
     * Returns this matrix as a {@link RealMatrix}.
     * <p>
     *     Changes int he return matrix will affect the content of this one.
     * </p>
     * @return never {@code null}.
     */
    default RealMatrix asRealMatrix() {
        return Access2DWrapper.of(new Access2D<Number>() {
            @Override
            public double doubleValue(long row, long col) {
                return LikelihoodMatrix.this.get((int) row, (int) col);
            }

            @Override
            public Number get(long row, long col) {
                return LikelihoodMatrix.this.get((int) row, (int) col);
            }

            @Override
            public long countColumns() {
                return evidenceCount();
            }

            @Override
            public long countRows() {
                return numberOfAlleles();
            }
        });
    }
}
