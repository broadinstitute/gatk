package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.param.ParamUtils;


/**
 * This class represents integer copy number states. It is intended to be used as states of a
 * hidden Markov model.
 */
public final class IntegerCopyNumberState  {

    /**
     * Integer value of the represented copy number state
     */
    private final int copyNumber;

    /**
     * A string representation of this copy number state (used for creating human-readable copy number call lists)
     */
    private final String callString;

    /**
     * An allele representation of this copy number state (used for VCF creation)
     */
    private final Allele allele;

    /**
     * The state with this copy number represents a reference allele
     */
    private static final int REFERENCE_COPY_NUMBER_VALUE = 2;

    public IntegerCopyNumberState(final int copyNumber) {
        this.copyNumber = ParamUtils.isPositiveOrZero(copyNumber, "The integer copy number state" +
                " must be non-negative");
        this.callString = toCallString(copyNumber);
        allele = createAlleleGivenCopyNumber(copyNumber);
    }

    public Allele toAllele() {
        return allele;
    }

    public String getCallString() {
        return callString;
    }

    public double getScalar() {
        return copyNumber;
    }

    public int getCopyNumber() { return copyNumber; }

    private static String toCallString(final int copyNumber) {
        return String.valueOf(copyNumber);
    }

    private static String toAlleleString(final int copyNumber) {
        if (copyNumber == REFERENCE_COPY_NUMBER_VALUE) {
            return "N";
        } else {
            return "<CN_" + String.valueOf(copyNumber) + ">";
        }
    }

    private static Allele createAlleleGivenCopyNumber(final int copyNumber) {
        if (copyNumber == REFERENCE_COPY_NUMBER_VALUE) {
            return Allele.create(toAlleleString(copyNumber), true);
        } else {
            return Allele.create(toAlleleString(copyNumber), false);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        return copyNumber == ((IntegerCopyNumberState) o).copyNumber;
    }

    @Override
    public int hashCode() {
        return copyNumber;
    }

    @Override
    public String toString() {
        return GermlineCNVNamingConstants.COPY_NUMBER_STATE_STRING_START + Integer.toString(copyNumber);
    }
}