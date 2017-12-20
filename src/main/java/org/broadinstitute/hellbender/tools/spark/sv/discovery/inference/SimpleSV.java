package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.IntrachromosomalBreakpointPair;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class SimpleSV implements Locatable {

    protected final SimpleSVType.TYPES type;
    protected final int start;
    protected final int end;
    protected final int contigId;
    protected final String contig;
    protected final int readPairEvidence;
    protected final int splitReadEvidence;
    protected final int readPairCounterEvidence;
    protected final int splitReadCounterEvidence;
    protected final List<CopyRatio> coverage;
    protected final List<Integer> copyNumberStates;
    protected final IntrachromosomalBreakpointPair breakpoints;

    public SimpleSV(final SimpleSVType.TYPES type,
                    final int start,
                    final int end,
                    final int contigId,
                    final String contig,
                    final int readPairEvidence,
                    final int splitReadEvidence,
                    final int readPairCounterEvidence,
                    final int splitReadCounterEvidence,
                    final List<CopyRatio> coverage,
                    final List<Integer> copyNumberStates,
                    final IntrachromosomalBreakpointPair breakpoints) {
        this.type = type;
        this.start = start;
        this.end = end;
        this.contigId = contigId;
        this.contig = contig;
        this.readPairEvidence = readPairEvidence;
        this.splitReadEvidence = splitReadEvidence;
        this.readPairCounterEvidence = readPairCounterEvidence;
        this.splitReadCounterEvidence = splitReadCounterEvidence;
        this.coverage = coverage;
        this.copyNumberStates = copyNumberStates;
        this.breakpoints = breakpoints;
    }

    private double getMedianCoverage(final List<CopyRatio> coverage) {
        if (coverage.isEmpty()) return 0;
        final List<Double> sortedCoverage = new ArrayList<>(coverage.stream().map(CopyRatio::getLog2CopyRatioValue).collect(Collectors.toList()));
        Collections.sort(sortedCoverage);
        final int numCounts = coverage.size();
        final int midIndex = numCounts / 2;
        if (numCounts % 2 == 0) return (sortedCoverage.get(midIndex - 1) + sortedCoverage.get(midIndex)) / 2.0;
        return sortedCoverage.get(midIndex);
    }

    private double getInterquartileRange(final List<CopyRatio> coverage) {
        if (coverage.isEmpty()) return 0;
        final List<Double> sortedCoverage = new ArrayList<>(coverage.stream().map(CopyRatio::getLog2CopyRatioValue).collect(Collectors.toList()));
        Collections.sort(sortedCoverage);
        final int size = sortedCoverage.size();
        final int leftIndex = size / 4;
        final int rightIndex = (size * 3) / 4;
        return Math.abs(sortedCoverage.get(rightIndex) - sortedCoverage.get(leftIndex));
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

    public String getContig() {
        return contig;
    }

    public int getReadPairEvidence() {
        return readPairEvidence;
    }

    public int getSplitReadEvidence() {
        return splitReadEvidence;
    }

    public List<CopyRatio> getCoverage() {
        return coverage;
    }

    public double medianCoverage() {
        return getMedianCoverage(coverage);
    }

    public double interquartileRange() {
        return getInterquartileRange(coverage);
    }

    public boolean isSupportedByBreakpoints() {
        return breakpoints != null;
    }

    public IntrachromosomalBreakpointPair getBreakpoints() {
        return breakpoints;
    }

    public int getSize() {
        return end - start;
    }

    public SimpleSVType.TYPES getType() {
        return type;
    }

    public String getString(final ReadMetadata readMetadata) {
        final String contigName = readMetadata.getContigName(contigId);
        final String bins = String.join(",", coverage.stream().map(bin -> String.valueOf(bin.getStart())).collect(Collectors.toList()));
        final String binCounts = String.join(",", coverage.stream().map(bin -> String.valueOf(bin.getLog2CopyRatioValue())).collect(Collectors.toList()));
        final String copyNumberStatesString = String.join(",", copyNumberStates.stream().map(String::valueOf).collect(Collectors.toList()));
        return contigName + "\t" + start + "\t" + end + "\t" + type + "\t" +
                readPairEvidence + "\t" + splitReadEvidence + "\t" + readPairCounterEvidence + "\t" + splitReadCounterEvidence + "\t" + medianCoverage() + "\t" +
                (isSupportedByBreakpoints() ? 1 : 0) + "\t" + bins + "\t" + binCounts + "\t" + copyNumberStatesString;
    }
}
