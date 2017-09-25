package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Class to represent and calculate the aligned contig score.
 */
class AlignmentScore {

    public static final int STRAND_SWITCH_COST = 60;
    public static final double MATCH_COST = 0.01;
    public static final int MISMATCH_COST = 30;
    public static final double GAP_OPEN_COST = 45;
    public static final double GAP_EXTEND_COST = 3;

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
        for (int i = 0; i < intervals.size(); i++) {
            final AlignmentInterval ai = intervals.get(i);
            if (i > 0) {
                final AlignmentInterval prev = intervals.get(i - 1);
                if (prev.forwardStrand != ai.forwardStrand) {
                    totalReversals++;
                } else {
                    final AlignmentInterval left = ai.forwardStrand ? prev : ai;
                    final AlignmentInterval right = ai.forwardStrand ? ai : prev;
                    if (left.referenceSpan.getEnd() < right.referenceSpan.getStart()) {
                        totalIndels++;
                        totalIndelLength += right.referenceSpan.getStart() - left.referenceSpan.getEnd();
                    }
                    if (left.endInAssembledContig < right.startInAssembledContig) {
                        totalIndels++;
                        totalIndelLength += right.startInAssembledContig - left.endInAssembledContig - 1;
                    }
                }
            }
            final int matches = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().isAlignment())
                    .mapToInt(CigarElement::getLength).sum();
            final int misMatches = ai.mismatches;
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
                totalIndelLength += intervals.get(0).startInAssembledContig - 1;
                totalIndels++;
            }
            if (intervals.get(intervals.size() - 1).endInAssembledContig < sequenceLength) {
                totalIndelLength += sequenceLength - intervals.get(intervals.size() - 1).endInAssembledContig;
                totalIndels++;
            }
        }
        return new AlignmentScore(totalMatches, totalMismatches, totalIndels, totalIndelLength, totalReversals);

    }


    public double getValue() {
        return -(int) Math.round(totalMatches * MATCH_COST
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
