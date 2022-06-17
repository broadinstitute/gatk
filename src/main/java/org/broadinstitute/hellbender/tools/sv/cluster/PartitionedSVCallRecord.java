package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

public final class PartitionedSVCallRecord extends SVCallRecord {
    private final boolean primaryVariant;

    public PartitionedSVCallRecord(final SVCallRecord record, final boolean primaryVariant) {
        super(record);
        this.primaryVariant = primaryVariant;
    }

    public static PartitionedSVCallRecord create(final SVCallRecord record, final boolean primaryVariant) {
        return new PartitionedSVCallRecord(record, primaryVariant);
    }

    public boolean isPrimaryVariant() {
        return primaryVariant;
    }
}
