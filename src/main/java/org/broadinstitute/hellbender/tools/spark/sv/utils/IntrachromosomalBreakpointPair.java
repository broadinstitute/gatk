package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;

public final class IntrachromosomalBreakpointPair {
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
}
