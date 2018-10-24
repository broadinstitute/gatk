package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Simple container class to represent bases that have been corrected for strandedness already.
 * That is, if they would occur on {@link htsjdk.tribble.annotation.Strand#NEGATIVE}, they have been reverse-complemented.
 * If they would occur on {@link htsjdk.tribble.annotation.Strand#POSITIVE}, they have not been altered.
 * Created by jonn on 10/24/18.
 */
public class StrandCorrectedReferenceBases {

    private final Strand strand;
    private final String strandCorrectedReferenceBases;

    public StrandCorrectedReferenceBases(final String bases, final Strand strand) {
        strandCorrectedReferenceBases = bases;
        this.strand = strand;
    }

    public StrandCorrectedReferenceBases(final byte[] bases, final Strand strand) {
        strandCorrectedReferenceBases = new String(bases);
        this.strand = strand;
    }

    public Strand getStrand() {
        return strand;
    }

    /**
     * Get the string of bases in this {@link StrandCorrectedReferenceBases} object as if it were on the given {@code strand}.
     * @param strand {@link Strand} to use as reference for retrieving the base string.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return A {@link String} of bases corresponding to the bases in this {@link StrandCorrectedReferenceBases} object if they were on the given {@code strand}.
     */
    public String getBaseString(final Strand strand) {
        FuncotatorUtils.assertValidStrand(strand);

        if ( strand == this.strand ) {
            return getBaseString();
        }
        else {
            return ReadUtils.getBasesReverseComplement(getBases());
        }
    }

    public String getBaseString() {
        return strandCorrectedReferenceBases;
    }

    public byte[] getBases() {
        return strandCorrectedReferenceBases.getBytes();
    }

    @Override
    public String toString() {
        return "StrandCorrectedReferenceBases{" +
                "strand=" + strand +
                ", strandCorrectedReferenceBases='" + strandCorrectedReferenceBases + '\'' +
                '}';
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final StrandCorrectedReferenceBases that = (StrandCorrectedReferenceBases) o;

        if ( strand != that.strand ) return false;
        return strandCorrectedReferenceBases != null ? strandCorrectedReferenceBases.equals(that.strandCorrectedReferenceBases) : that.strandCorrectedReferenceBases == null;
    }

    @Override
    public int hashCode() {
        int result = strand != null ? strand.hashCode() : 0;
        result = 31 * result + (strandCorrectedReferenceBases != null ? strandCorrectedReferenceBases.hashCode() : 0);
        return result;
    }
}
