package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

import java.io.Serializable;
import java.util.Collection;
import java.util.Objects;

/**
 * Represents a simple structural variant (e.g. deletion, duplication, etc.) and associated evidence
 */
public class LargeSimpleSV implements Serializable {

    public static final long serialVersionUID = 1L;
    protected final SimpleSVType.TYPES type;
    protected final int start;
    protected final int end;
    protected final int contigId;
    protected final int readPairEvidence;

    protected final int splitReadEvidence;
    protected final int readPairCounterEvidence;
    protected final int splitReadCounterEvidence;
    protected final Collection<EvidenceTargetLink> supportingEvidence;
    protected final IntrachromosomalBreakpointPair breakpoints;

    public LargeSimpleSV(final SimpleSVType.TYPES type,
                         final int start,
                         final int end,
                         final int contigId,
                         final int readPairEvidence,
                         final int splitReadEvidence,
                         final int readPairCounterEvidence,
                         final int splitReadCounterEvidence,
                         final IntrachromosomalBreakpointPair breakpoints,
                         final Collection<EvidenceTargetLink> supportingEvidence) {
        this.type = type;
        this.start = start;
        this.end = end;
        this.contigId = contigId;
        this.readPairEvidence = readPairEvidence;
        this.splitReadEvidence = splitReadEvidence;
        this.readPairCounterEvidence = readPairCounterEvidence;
        this.splitReadCounterEvidence = splitReadCounterEvidence;
        this.breakpoints = breakpoints;
        this.supportingEvidence = supportingEvidence;
    }

    public double getScore(final double counterEvidencePseudocount) {
        return computeScore(readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, counterEvidencePseudocount, breakpoints != null);
    }

    public static double computeScore(final int readPairEvidence, final int splitReadEvidence, final int readPairCounterEvidence, final int splitReadCounterEvidence, final double counterEvidencePseudocount, final boolean hasBreakpoints) {
        return (readPairEvidence + splitReadEvidence) / Math.max(readPairCounterEvidence + splitReadCounterEvidence, counterEvidencePseudocount) + (hasBreakpoints ? 0.001 : 0);
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

    public Collection<EvidenceTargetLink> getSupportingEvidence() { return supportingEvidence; }

    public int getSize() {
        return end - start;
    }

    public SimpleSVType.TYPES getType() {
        return type;
    }

    public LargeSimpleSV copy() {
        final LargeSimpleSV copy = new LargeSimpleSV(type, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints, supportingEvidence);
        return copy;
    }

    public static String getBedHeader() {
        return "CONTIG\tSTART\tEND\tLEN\tTYPE\tE_RP\tE_SR\tCE_RP\tCE_SR\tSCORE\tBRKPTS";
    }

    public String toBedString(final SAMSequenceDictionary dictionary, final double counterEvidencePseudocount) {
        return dictionary.getSequence(contigId).getSequenceName() + "\t" + start + "\t" + end + "\t" + (end-start) + "\t" + type + "\t" + readPairEvidence + "\t" + splitReadEvidence +
                "\t" + readPairCounterEvidence + "\t" + splitReadCounterEvidence + "\t" + getScore(counterEvidencePseudocount) + "\t" + (breakpoints == null ? "NONE" : breakpoints.getString(dictionary));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof LargeSimpleSV)) return false;
        LargeSimpleSV that = (LargeSimpleSV) o;
        return start == that.start &&
                end == that.end &&
                contigId == that.contigId &&
                readPairEvidence == that.readPairEvidence &&
                splitReadEvidence == that.splitReadEvidence &&
                readPairCounterEvidence == that.readPairCounterEvidence &&
                splitReadCounterEvidence == that.splitReadCounterEvidence &&
                type == that.type &&
                Objects.equals(breakpoints, that.breakpoints);
    }

    @Override
    public int hashCode() {
        return Objects.hash(type, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints);
    }
}
