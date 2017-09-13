package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.utils.Utils;

public final class PairedStrandedIntervals {
    final StrandedInterval left;
    final StrandedInterval right;

    public PairedStrandedIntervals(final StrandedInterval left, final StrandedInterval right) {
        Utils.validate(left != null, "Can't construct PairedStrandedInterval with a null left interval");
        Utils.validate(right != null, "Can't construct PairedStrandedInterval with a null right interval");
        this.left = left;
        this.right = right;
    }

    public StrandedInterval getLeft() {
        return left;
    }

    public StrandedInterval getRight() {
        return right;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final PairedStrandedIntervals that = (PairedStrandedIntervals) o;

        if (!left.equals(that.left)) return false;
        return right.equals(that.right);
    }

    @Override
    public int hashCode() {
        int result = left.hashCode();
        result = 31 * result + right.hashCode();
        return result;
    }
}
