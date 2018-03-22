package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

public class LargeSimpleSV {

    protected final SimpleSVType.TYPES type;
    protected final int start;
    protected final int end;
    protected final int contigId;
    protected final int readPairEvidence;
    protected final int splitReadEvidence;
    protected final int readPairCounterEvidence;
    protected final int splitReadCounterEvidence;
    protected final IntrachromosomalBreakpointPair breakpoints;

    public LargeSimpleSV(final SimpleSVType.TYPES type,
                         final int start,
                         final int end,
                         final int contigId,
                         final int readPairEvidence,
                         final int splitReadEvidence,
                         final int readPairCounterEvidence,
                         final int splitReadCounterEvidence,
                         final IntrachromosomalBreakpointPair breakpoints) {
        this.type = type;
        this.start = start;
        this.end = end;
        this.contigId = contigId;
        this.readPairEvidence = readPairEvidence;
        this.splitReadEvidence = splitReadEvidence;
        this.readPairCounterEvidence = readPairCounterEvidence;
        this.splitReadCounterEvidence = splitReadCounterEvidence;
        this.breakpoints = breakpoints;
    }

    public double getScore(final double counterEvidencePseudocount) {
        return computeScore(readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, counterEvidencePseudocount);
    }

    public static double computeScore(final int readPairEvidence, final int splitReadEvidence, final int readPairCounterEvidence, final int splitReadCounterEvidence, final double counterEvidencePseudocount) {
        return (readPairEvidence + splitReadEvidence) / Math.max(readPairCounterEvidence + splitReadCounterEvidence, counterEvidencePseudocount);
    }

    public SVInterval getInterval() {
        return new SVInterval(contigId, start, end);
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getContigId() {
        return contigId;
    }

    public int getReadPairEvidence() {
        return readPairEvidence;
    }

    public int getSplitReadEvidence() {
        return splitReadEvidence;
    }

    public IntrachromosomalBreakpointPair getBreakpoints() { return breakpoints; }

    public int getSize() {
        return end - start;
    }

    public SimpleSVType.TYPES getType() {
        return type;
    }

    public String getString(final SAMSequenceDictionary dictionary, final double counterEvidencePseudocount) {
        return dictionary.getSequence(contigId).getSequenceName() + "\t" + start + "\t" + end + "\t" + type + "\t" + readPairEvidence + "\t" + splitReadEvidence +
                "\t" + readPairCounterEvidence + "\t" + splitReadCounterEvidence + "\t" + getScore(counterEvidencePseudocount) + "\t" + (breakpoints == null ? "none" : breakpoints.getString(dictionary));
    }
}
