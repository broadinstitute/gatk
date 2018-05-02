package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.testng.collections.Lists;

import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * Represents a simple structural variant (e.g. deletion, duplication, etc.) and associated evidence
 */
public final class LargeSimpleSV implements Serializable {

    public enum EvidenceType {
        CALL_IMPRECISE,
        CALL_PRECISE,
        BREAKPOINT_PAIR,
        READ_PAIR
    }

    public static final long serialVersionUID = 1L;

    protected final SimpleSVType.TYPES eventType;
    protected final EvidenceType evidenceType;

    protected final int start;
    protected final int end;
    protected final int contigId;

    protected final int readPairEvidence;
    protected final int splitReadEvidence;
    protected final int readPairCounterEvidence;
    protected final int splitReadCounterEvidence;

    protected final Collection<EvidenceTargetLink> supportingEvidence;
    protected final IntrachromosomalBreakpointPair breakpoints;

    public LargeSimpleSV(final SimpleSVType.TYPES eventType,
                         final EvidenceType evidenceType,
                         final int start,
                         final int end,
                         final int contigId,
                         final int readPairEvidence,
                         final int splitReadEvidence,
                         final int readPairCounterEvidence,
                         final int splitReadCounterEvidence,
                         final IntrachromosomalBreakpointPair breakpoints,
                         final Collection<EvidenceTargetLink> supportingEvidence) {
        this.eventType = eventType;
        this.evidenceType = evidenceType;
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

    public SVInterval getInterval() {
        return new SVInterval(contigId, start, end);
    }

    public EvidenceType getEvidenceType() {
        return evidenceType;
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

    public SimpleSVType.TYPES getEventType() {
        return eventType;
    }

    public LargeSimpleSV copy() {
        final LargeSimpleSV copy = new LargeSimpleSV(eventType, evidenceType, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints, supportingEvidence);
        return copy;
    }

    public static String getBedHeader() {
        return "CONTIG\tSTART\tEND\tLEN\tEVENT_TYPE\tEVIDENCE_TYPE\tE_RP\tE_SR\tCE_RP\tCE_SR\tBRKPTS";
    }

    public String toBedString(final SAMSequenceDictionary dictionary) {
        return dictionary.getSequence(contigId).getSequenceName() + "\t" + start + "\t" + end + "\t" + (end-start) + "\t" + eventType + "\t" + evidenceType + "\t" + readPairEvidence + "\t" + splitReadEvidence +
                "\t" + readPairCounterEvidence + "\t" + splitReadCounterEvidence + "\t" + (breakpoints == null ? "NONE" : breakpoints.getString(dictionary));
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
                eventType == that.eventType &&
                Objects.equals(breakpoints, that.breakpoints);
    }

    @Override
    public int hashCode() {
        return Objects.hash(eventType, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpoints);
    }

    //TODO encode and decode supporting evidence and breakpoint assembly contigs
    public String encode() {
        final List<Integer> breakpointData = breakpoints == null ? Lists.newArrayList(0, 0, 0) : Lists.newArrayList(breakpoints.getContig(), breakpoints.getFirst(), breakpoints.getSecond());
        final List<String> breakpointDataStrings = breakpointData.stream().map(String::valueOf).collect(Collectors.toList());
        return eventType.name() + "\t" + evidenceType.name() + "\t" + contigId + "\t" + start + "\t" + end + "\t" + readPairEvidence + "\t" + splitReadEvidence + "\t" + readPairCounterEvidence + "\t" + splitReadCounterEvidence + "\t" + String.join(",", breakpointDataStrings);

    }

    public static LargeSimpleSV decode(final String line) {
        final String[] tokens = line.trim().split("\t");
        final SimpleSVType.TYPES eventType = SimpleSVType.TYPES.valueOf(tokens[0]);
        final EvidenceType evidenceType = EvidenceType.valueOf(tokens[1]);
        final int contigId = Integer.valueOf(tokens[2]);
        final int start = Integer.valueOf(tokens[3]);
        final int end = Integer.valueOf(tokens[4]);
        final int readPairEvidence = Integer.valueOf(tokens[5]);
        final int splitReadEvidence = Integer.valueOf(tokens[6]);
        final int readPairCounterEvidence = Integer.valueOf(tokens[7]);
        final int splitReadCounterEvidence = Integer.valueOf(tokens[8]);
        final String[] breakpointData = tokens[9].split(",");
        final IntrachromosomalBreakpointPair breakpointPair = new IntrachromosomalBreakpointPair(Integer.valueOf(breakpointData[0]), Integer.valueOf(breakpointData[1]), Integer.valueOf(breakpointData[2]), Collections.emptyList(), Collections.emptyList());
        return new LargeSimpleSV(eventType, evidenceType, start, end, contigId, readPairEvidence, splitReadEvidence, readPairCounterEvidence, splitReadCounterEvidence, breakpointPair, Collections.emptyList());
    }
}
