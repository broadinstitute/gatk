package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.annotation.Strand;

/**
 * Class to represent a strand-corrected {@link htsjdk.variant.variantcontext.Allele}.
 * Created by jonn on 10/24/18.
 */
public class StrandCorrectedAllele extends htsjdk.variant.variantcontext.Allele {
    public static final long serialVersionUID = 1L;

    /**
     * The strand on which this allele occurs.
     */
    private final Strand strand;

    private StrandCorrectedAllele(final String allele, final boolean isRef, final Strand strand) {
        super(allele, isRef);

        FuncotatorUtils.assertValidStrand(strand);
        this.strand = strand;
    }

    /**
     * Create a new {@link StrandCorrectedAllele} object with the given bases.
     * {@code bases} are assumed to be on the {@link Strand#POSITIVE} strand.
     * @param bases {@link String} of bases to use as the allele.
     * @return A new {@link StrandCorrectedAllele} object containing the given bases.
     */
    public static StrandCorrectedAllele create(final String bases) {
        return new StrandCorrectedAllele(bases, false, Strand.POSITIVE);
    }

    /**
     * Create a new {@link StrandCorrectedAllele} object with the given bases.
     * {@code bases} are assumed to be strand-corrected already.
     * @param bases {@link String} of bases to use as the allele.
     * @param strand The {@link Strand} on which the allele occurs.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return A new {@link StrandCorrectedAllele} object containing the given bases.
     */
    public static StrandCorrectedAllele create(final String bases, final Strand strand) {
        return new StrandCorrectedAllele(bases, false, strand);
    }

    /**
     * Create a new {@link StrandCorrectedAllele} object with the given bases.
     * {@code bases} are assumed to be on the {@link Strand#POSITIVE} strand.
     * @param bases {@link String} of bases to use as the allele.
     * @param isRef {@code true} iff the given allele is a reference allele.  {@code false} otherwise.
     * @return A new {@link StrandCorrectedAllele} object containing the given bases.
     */
    public static StrandCorrectedAllele create(final String bases, final boolean isRef) {
        return new StrandCorrectedAllele(bases, isRef, Strand.POSITIVE);
    }

    /**
     * Create a new {@link StrandCorrectedAllele} object with the given bases.
     * {@code bases} are assumed to be strand-corrected already.
     * @param bases {@link String} of bases to use as the allele.
     * @param isRef {@code true} iff the given allele is a reference allele.  {@code false} otherwise.
     * @param strand The {@link Strand} on which the allele occurs.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return A new {@link StrandCorrectedAllele} object containing the given bases.
     */
    public static StrandCorrectedAllele create(final String bases, final boolean isRef, final Strand strand) {
        return new StrandCorrectedAllele(bases, isRef, strand);
    }

    /**
     * @return The {@link Strand} on which this allele occurs.
     */
    public Strand getStrand() {
        return strand;
    }

    @Override
    public String toString() {
        return "StrandCorrectedAllele{" +
                "strand=" + strand +
                " " + super.toString() + " }";
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;
        if ( !super.equals(o) ) return false;

        final StrandCorrectedAllele that = (StrandCorrectedAllele) o;

        return strand == that.strand;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (strand != null ? strand.hashCode() : 0);
        return result;
    }
}
