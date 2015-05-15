package org.broadinstitute.hellbender.utils.reference;

import org.broadinstitute.hellbender.utils.SimpleInterval;

public class ReferenceBases {

    private final byte[] bases;
    private final SimpleInterval interval;

    public ReferenceBases( final byte[] bases, final SimpleInterval interval ) {
        this.bases = bases;
        this.interval = interval;
    }

    public byte[] getBases() {
        return bases;
    }

    public SimpleInterval getInterval() {
        return interval;
    }
}
