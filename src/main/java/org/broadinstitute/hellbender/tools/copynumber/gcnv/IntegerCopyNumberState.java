package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * This class represents integer copy number states.
 */
public final class IntegerCopyNumberState  {

    /**
     * Integer value of the represented copy number state
     */
    private final int copyNumber;

    /**
     * An allele representation of this copy number state (used for VCF creation)
     */



    public IntegerCopyNumberState(final int copyNumber) {
        this.copyNumber = ParamUtils.isPositiveOrZero(copyNumber, "The integer copy number state" +
                " must be non-negative");
    }

    public Allele toAllele(final int refCopyNumber) {
        return createAlleleGivenCopyNumber(copyNumber, refCopyNumber);
    }

    public int getCopyNumber() { return copyNumber; }

    private static String toAlleleString(final int copyNumber, final int refCopyNumber) {
        if (copyNumber == refCopyNumber) {
            return "N";
        } else {
            return "<CN_" + String.valueOf(copyNumber) + ">";
        }
    }

    private static Allele createAlleleGivenCopyNumber(final int copyNumber,
                                                      final int refCopyNumber) {
        if (copyNumber == refCopyNumber) {
            return Allele.create(toAlleleString(copyNumber, refCopyNumber), true);
        } else {
            return Allele.create(toAlleleString(copyNumber, refCopyNumber), false);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        return copyNumber == ((IntegerCopyNumberState) o).copyNumber;
    }

    @Override
    public int hashCode() {
        return copyNumber;
    }

    @Override
    public String toString() {
        return GermlineCNVNamingConstants.COPY_NUMBER_TABLE_COLUMN_PREFIX + Integer.toString(copyNumber);
    }
}