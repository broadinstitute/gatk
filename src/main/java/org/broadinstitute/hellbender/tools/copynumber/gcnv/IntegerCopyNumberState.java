package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;

/**
 * This class represents integer copy number states.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberState implements Comparable<IntegerCopyNumberState> {
    /**
     * Integer value of the represented copy number state
     */
    private final int copyNumber;

    public IntegerCopyNumberState(final int copyNumber) {
        this.copyNumber = ParamUtils.isPositiveOrZero(copyNumber, "The integer copy number state" +
                " must be non-negative");
    }

    public int getCopyNumber() { return copyNumber; }

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
    public int compareTo(@Nonnull final IntegerCopyNumberState other) {
        return copyNumber - other.copyNumber;
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