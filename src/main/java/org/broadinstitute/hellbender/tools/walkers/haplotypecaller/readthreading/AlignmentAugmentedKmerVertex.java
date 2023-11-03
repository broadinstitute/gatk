package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.apache.commons.lang.math.IntRange;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseVertex;

public class AlignmentAugmentedKmerVertex extends BaseVertex {
    private final IntRange alignmentRange;

    public AlignmentAugmentedKmerVertex(final byte[] sequence, final IntRange alignmentRange) {
        super(sequence);
        this.alignmentRange = alignmentRange;
    }

    public boolean compatibleWith(final AlignmentAugmentedKmerVertex other) {
        return seqEquals(other) && alignmentRange.overlapsRange(other.alignmentRange);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AlignmentAugmentedKmerVertex that = (AlignmentAugmentedKmerVertex) o;
        if (hashCode() != that.hashCode()){
            return false;
        }

        return seqEquals(that) && this.alignmentRange.equals(that.alignmentRange);
    }

    @Override
    public int hashCode() {
        return cachedHashCode + alignmentRange.hashCode();
    }

    @Override
    public String toString() {
        return getSequenceString();
    }
}
