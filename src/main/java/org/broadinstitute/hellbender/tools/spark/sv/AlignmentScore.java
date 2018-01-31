package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class to represent and calculate the aligned contig score.
 */
class AlignmentScore {

    public static final int STRAND_SWITCH_COST = 60;
    public static final double MATCH_COST = 0.01;
    public static final int MISMATCH_COST = 40;
    public static final double GAP_OPEN_COST = 60;
    public static final double GAP_EXTEND_COST = 10;

    final int totalReversals;
    final int totalIndels;
    final int totalMatches;
    final int totalMismatches;
    final int totalIndelLength;

    public AlignmentScore(final int matches, final int mismatches, final int indels, final int totalIndelLength, final int reversals) {
        this.totalReversals = reversals;
        this.totalIndels = indels;
        this.totalMatches = matches;
        this.totalMismatches = mismatches;
        this.totalIndelLength = totalIndelLength;
    }

    public static AlignmentScore valueOf(final String str) {
        final String[] parts = Utils.nonNull(str).split("[:,]");
        Utils.validateArg(parts.length == 6, "the input string has the wrong number of components");
        int nextIdx = 0;
        // we only check that the first value is a valid double (score), we don't keep this value as is derivable from
        // the other based on penalties.
        ParamUtils.isDouble(parts[nextIdx++], "score is not a valid double value");
        final int matches = ParamUtils.isPositiveOrZeroInteger(parts[nextIdx++], "matches is not a valid positive integer");
        final int misMatches = ParamUtils.isPositiveOrZeroInteger(parts[nextIdx++], "misMatches is not a valid positive integer");
        final int indels = ParamUtils.isPositiveOrZeroInteger(parts[nextIdx++], "indels is not a valid positive integer");
        final int indelLenghts = ParamUtils.isPositiveOrZeroInteger(parts[nextIdx++], "indel-length is not a valid positive integer");
        final int reversals = ParamUtils.isPositiveOrZeroInteger(parts[nextIdx], "reversals is not a valid positive integer");
        return new AlignmentScore(matches, misMatches, indels, indelLenghts, reversals);
    }

    public static AlignmentScore calculate(final int sequenceLength, final List<AlignmentInterval> alignmentIntervals) {
        final List<AlignmentInterval> intervals = alignmentIntervals.stream()
                .sorted(Comparator.comparing(ai -> ai.startInAssembledContig))
                .collect(Collectors.toList());
        int totalReversals = 0;
        int totalIndels = 0;
        int totalMatches = 0;
        int totalMismatches = 0;
        int totalIndelLength = 0;
        final int forwardAlignedBases = intervals.stream().filter(ai -> ai.forwardStrand && ai.mapQual > 15).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
        final int reverseAlignedBases = intervals.stream().filter(ai -> !ai.forwardStrand && ai.mapQual > 15).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
        final int direction = -Integer.compare(reverseAlignedBases, forwardAlignedBases);
        for (int i = 0; i < intervals.size(); i++) {
            final AlignmentInterval ai = intervals.get(i);
            if (i > 0) {
                final AlignmentInterval prev = intervals.get(i - 1);
                if (prev.forwardStrand != ai.forwardStrand) {
                    totalReversals++;
                    if (direction != (ai.forwardStrand ? 1 : -1)) {
                       totalIndelLength += CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig);
                       totalIndels++;
                    }
                } else {
                    final AlignmentInterval left = ai.forwardStrand ? prev : ai;
                    final AlignmentInterval right = ai.forwardStrand ? ai : prev;
                    final int refIndelLength = right.referenceSpan.getStart() - left.referenceSpan.getEnd();
                    final int ctgIndelLength = prev.endInAssembledContig - ai.startInAssembledContig;
                    if (refIndelLength != 1) {
                        totalIndels++;
                        totalIndelLength += Math.abs(refIndelLength - 1);
                    }
                    if (ctgIndelLength != 1) {
                        totalIndels++;
                        totalIndelLength += Math.abs(ctgIndelLength - 1);
                    }
                }
            }
            final int matches = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isAlignment())
                    .mapToInt(CigarElement::getLength).sum();
            final int misMatches = ai.getBaseMismatches();
            final int indelCount = (int) ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isIndel())
                    .count();
            final int indelLengthSum = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isIndel())
                    .mapToInt(CigarElement::getLength).sum();
            totalIndels += indelCount;
            totalMatches += matches;
            totalMismatches += misMatches;
            totalIndelLength += indelLengthSum;
        }
        if (intervals.isEmpty()) {
            totalIndelLength += sequenceLength;
            totalIndels++;
        } else {
            if (intervals.get(0).startInAssembledContig > 1) {
                totalIndelLength += Math.min(intervals.get(0).startInAssembledContig - 1, intervals.get(0).referenceSpan.getStart() - 1);
                totalIndels++;
            }
            if (intervals.get(intervals.size() - 1).endInAssembledContig < sequenceLength) {
                totalIndelLength += sequenceLength - intervals.get(intervals.size() - 1).endInAssembledContig;
                totalIndels++;
            }
        }
        return new AlignmentScore(totalMatches, totalMismatches, totalIndels, totalIndelLength, totalReversals);

    }

    public static AlignmentScore calculate(final byte[] ref, final byte[] seq, final List<AlignmentInterval> alignmentIntervals) {
        final List<AlignmentInterval> intervals = alignmentIntervals.stream()
                .sorted(Comparator.comparing(ai -> ai.startInAssembledContig))
                .collect(Collectors.toList());
        int totalReversals = 0;
        int totalIndels = 0;
        int totalMatches = 0;
        int totalMismatches = 0;
        int totalIndelLength = 0;
        final int forwardAlignedBases = intervals.stream().filter(ai -> ai.forwardStrand && ai.mapQual > 15).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
        final int reverseAlignedBases = intervals.stream().filter(ai -> !ai.forwardStrand && ai.mapQual > 15).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
        final int direction = -Integer.compare(reverseAlignedBases, forwardAlignedBases);
        for (int i = 0; i < intervals.size(); i++) {
            final AlignmentInterval ai = intervals.get(i);
            if (i > 0) {
                final AlignmentInterval prev = intervals.get(i - 1);
                if (prev.forwardStrand != ai.forwardStrand) {
                    totalReversals++;
                   if (direction != (ai.forwardStrand ? 1 : -1)) {
                        totalIndelLength += CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig);
                   //     totalIndels++;
                   }
                } else {
                    final AlignmentInterval left = ai.forwardStrand ? prev : ai;
                    final AlignmentInterval right = ai.forwardStrand ? ai : prev;
                    final int refIndelLength = right.referenceSpan.getStart() - left.referenceSpan.getEnd();
                    final int ctgIndelLength = prev.endInAssembledContig - ai.startInAssembledContig;
                    if (refIndelLength != 1) {
                        totalIndels++;
                        totalIndelLength += Math.abs(refIndelLength - 1);
                    }
                    if (ctgIndelLength != 1) {
                        totalIndels++;
                        totalIndelLength += Math.abs(ctgIndelLength - 1);
                    }
                }
            }
            final int matches = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isAlignment())
                    .mapToInt(CigarElement::getLength).sum();
            final int misMatches = calculateMismatches(ref, seq, ai);
            final int indelCount = (int) ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isIndel())
                    .count();
            final int indelLengthSum = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isIndel())
                    .mapToInt(CigarElement::getLength).sum();
            totalIndels += indelCount;
            totalMatches += matches;
            totalMismatches += misMatches;
            totalIndelLength += indelLengthSum;
        }
        if (intervals.isEmpty()) {
            totalIndelLength += seq.length;
            totalIndels++;
        } else {
            if (intervals.get(0).startInAssembledContig > 1) {
                totalIndelLength += Math.min(intervals.get(0).startInAssembledContig - 1, intervals.get(0).referenceSpan.getStart() - 1);
                totalIndels++;
            }
            if (intervals.get(intervals.size() - 1).endInAssembledContig < seq.length) {
                totalIndelLength += seq.length - intervals.get(intervals.size() - 1).endInAssembledContig;
                totalIndels++;
            }
        }
        return new AlignmentScore(totalMatches, totalMismatches, totalIndels, totalIndelLength, totalReversals);
    }

    private static int calculateMismatches(final byte[] ref, final byte[] seq, final AlignmentInterval ai) {
        int refOffset = ai.referenceSpan.getStart() - 1;
        int direction = ai.forwardStrand ? 1 : -1;
//        int seqOffset = direction == 1 ? (ai.startInAssembledContig - 1 - CigarUtils.countLeftHardClippedBases(ai.cigarAlong5to3DirectionOfContig))
//                : (ai.endInAssembledContig - 1 - CigarUtils.countRightHardClippedBases(ai.cigarAlong5to3DirectionOfContig));
        int seqOffset = direction == 1 ? (ai.startInAssembledContig - 1) : (ai.endInAssembledContig - 1);
        final List<CigarElement> elements = ai.cigarAlongReference().getCigarElements();
        int index = 0;
        while (index < elements.size() && elements.get(index).getOperator().isClipping()) {
            index++;
        }
        int result = 0;
        while (index < elements.size() && !elements.get(index).getOperator().isClipping()) {
            final CigarElement element = elements.get(index++);
            if (element.getOperator().isAlignment()) {
                for (int i = 0, j = 0; i < element.getLength(); i++, j += direction) {
                    if (ai.forwardStrand && Nucleotide.decode(ref[refOffset + i]) != Nucleotide.decode(seq[seqOffset + j])) {
                        result++;
                    } else if (!ai.forwardStrand && !Nucleotide.decode(ref[refOffset + i]).isComplementOf(Nucleotide.decode(seq[seqOffset + j]))) {
                        result++;
                    }
                }
            }
            if (element.getOperator().consumesReferenceBases()) refOffset += element.getLength();
            if (element.getOperator().consumesReadBases()) seqOffset += element.getLength() * direction;
        }
        return result;
    }


    public double getValue() {
        return -0.1 * Math.round(totalMatches * MATCH_COST
                + totalMismatches * MISMATCH_COST
                + totalIndels * GAP_OPEN_COST
                + (totalIndelLength - totalIndels) * GAP_EXTEND_COST
                + totalReversals * STRAND_SWITCH_COST);
    }

    public String toString() {
        return  getValue() + ":" + Utils.join(",", totalMatches, totalMismatches,
                totalIndels, totalIndelLength, totalReversals);
    }
}
