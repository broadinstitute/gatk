package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;

import java.util.Collection;
import java.util.function.Function;

public class SVCallRecordDeduplicator<T extends SVCallRecord> extends SVDeduplicator<T> {

    public SVCallRecordDeduplicator(final Function<Collection<T>,T> collapser, final SAMSequenceDictionary dictionary) {
        super(collapser, dictionary);
    }

    @Override
    public boolean itemsAreIdentical(final T a, final T b) {
        return super.itemsAreIdentical(a, b)
                && a.getType().equals(b.getType())
                && a.getStrandA() == b.getStrandA()
                && a.getStrandB() == b.getStrandB();
    }
}
