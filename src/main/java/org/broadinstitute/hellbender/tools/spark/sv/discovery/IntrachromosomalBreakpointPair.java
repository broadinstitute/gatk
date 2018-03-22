package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

import java.util.List;

public final class IntrachromosomalBreakpointPair {
    private final int contig;
    private final int first;
    private final int second;
    private final List<String> firstContigs;
    private final List<String> secondContigs;

    public IntrachromosomalBreakpointPair(final int contig, final int first, final int second, final List<String> firstContigs, final List<String> secondContigs) {
        this.contig = contig;
        this.first = first;
        this.second = second;
        this.firstContigs = firstContigs;
        this.secondContigs = secondContigs;
    }

    public List<String> getFirstContigs() {
        return firstContigs;
    }

    public List<String> getSecondContigs() {
        return secondContigs;
    }

    public int getContig() {
        return contig;
    }

    public int getFirst() {
        return first;
    }

    public int getSecond() { return second; }

    public String getString(final SAMSequenceDictionary dictionary) {
        return dictionary.getSequence(contig).getSequenceName() + "_" + first + "_" + second;
    }

    public SVInterval getInterval() {
        if (first < second) {
            return new SVInterval(contig, first, second);
        }
        return new SVInterval(contig, second, first);
    }
}
