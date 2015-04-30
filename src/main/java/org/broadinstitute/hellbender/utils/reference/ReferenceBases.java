package org.broadinstitute.hellbender.utils.reference;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.Arrays;

public class ReferenceBases implements Serializable {

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

    public ReferenceBases getSubset(SimpleInterval interval) {
        // I don't think this is quite right...

        // First check bounds, then if it's a proper subset, then return that subset.

        // Correct check?
        if (!this.interval.contains(interval)) {
            throw new GATKException("Reference doesn't match read...");
        }
        int start = interval.getStart() - this.interval.getStart();
        int end = interval.getEnd() - this.interval.getStart();
        return new ReferenceBases(Arrays.copyOfRange(this.bases, start, end), interval);
    }
}
