package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Collection;
import java.util.Objects;

public final class IntrachromosomalBreakpointPair implements Serializable {
    public static final long serialVersionUID = 1L;
    private final int contig;
    private final int first;
    private final int second;
    private final Collection<String> firstAssembledContigs;
    private final Collection<String> secondAssembledContigs;

    public IntrachromosomalBreakpointPair(final int contig, final int first, final int second,
                                          final Collection<String> firstAssembledContigs,
                                          final Collection<String> secondAssembledContigs) {
        Utils.nonNull(firstAssembledContigs, "First breakpoint contig list cannot be null");
        Utils.nonNull(secondAssembledContigs, "Second breakpoint contig list cannot be null");
        this.contig = contig;
        this.first = first;
        this.second = second;
        this.firstAssembledContigs = firstAssembledContigs;
        this.secondAssembledContigs = secondAssembledContigs;
    }

    public Collection<String> getFirstAssembledContigs() {
        return firstAssembledContigs;
    }

    public Collection<String> getSecondAssembledContigs() {
        return secondAssembledContigs;
    }

    public int getContig() {
        return contig;
    }

    public int getFirst() {
        return first;
    }

    public int getSecond() { return second; }

    public String getString(final SAMSequenceDictionary dictionary) {
        final SAMSequenceRecord record = dictionary.getSequence(contig);
        if (record == null) {
            throw new IllegalArgumentException("Sequence with index " + contig + " could not be found in the sequence dictionary");
        }
        return record.getSequenceName() + "_" + first + "_" + second;
    }

    public SVInterval getInterval() {
        if (first < second) {
            return new SVInterval(contig, first, second);
        }
        return new SVInterval(contig, second, first);
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof IntrachromosomalBreakpointPair)) return false;
        if (this == o) return true;
        IntrachromosomalBreakpointPair that = (IntrachromosomalBreakpointPair) o;
        return contig == that.contig &&
                first == that.first &&
                second == that.second &&
                Objects.equals(firstAssembledContigs, that.firstAssembledContigs) &&
                Objects.equals(secondAssembledContigs, that.secondAssembledContigs);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contig, first, second, firstAssembledContigs, secondAssembledContigs);
    }
}
