package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.io.Serializable;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class to represent and calculate the aligned contig score.
 */
@DefaultSerializer(RealignmentScore.Serializer.class)
public final class RealignmentScore {

   // private static final long serialVersionUID = 1L;

    /**
     * Phred scaled penalty for matched base calls:
     * <br/>
     * This is typically a small number such as 0.01
     */
    public final double matchPenalty;

    /**
     * Phred scaled penalty for mismatched base calls:
     * <br/>
     * This is typically a medium large number e.g. 30.
     */
    public final double mismatchPenalty;

    /**
     * Phred scaled penalty for a indel start.
     * <br/>
     * This is typically a large number e.g. 45
     */
    public final double gapOpenPenalty;

    /**
     * Phred scaled penalty for a indel 1bp extension.
     * <br/>
     * This is typically a small number e.g 10.
     */
    public final double gapExtendPenalty;

    /**
     * Number of indels found.
     */
    public final int numberOfIndels;

    /**
     * Number of base-call matches found.
     */
    public final int numberOfMatches;

    /**
     * Number of base-call mismatches found.
     */
    public final int numberOfMismatches;

    /**
     * Total length of indels (open + extensions).
     */
    public final int indelLengthSum;

    /**
     * Indicate the strand the record is supposed to align too.
     */
    public final Strand strand;

    /**
     * The total score value in Phred scale.
     */
    private final double value;

    public RealignmentScore(final RealignmentScoreParameters parameters, final int matches, final int mismatches, final int indels, final int totalIndelLength, final Strand strand) {
        ParamUtils.isPositiveOrZero(matches, "number of matches cannot be negative");
        ParamUtils.isPositiveOrZero(mismatches, "number of mismatches cannot be negative");
        ParamUtils.isPositiveOrZero(indels, "number of indels cannot be negative");
        ParamUtils.isPositiveOrZero(totalIndelLength - indels, "total length of indels minus the number of indels cannot be negative");
        Utils.nonNull(strand);
        this.gapExtendPenalty = parameters.gapExtendPenalty;
        this.gapOpenPenalty = parameters.gapOpenPenalty;
        this.mismatchPenalty = parameters.mismatchPenalty;
        this.matchPenalty = parameters.matchPenalty;
        this.numberOfIndels = indels;
        this.numberOfMatches = matches;
        this.numberOfMismatches = mismatches;
        this.indelLengthSum = totalIndelLength;
        this.strand = strand;
        if (this.numberOfMatches + this.numberOfMismatches == 0) {
            this.value = Double.POSITIVE_INFINITY;
        } else {
            this.value = numberOfIndels * gapOpenPenalty
                   + numberOfMatches * matchPenalty
                   + numberOfMismatches * mismatchPenalty
                   + (indelLengthSum - numberOfIndels) * gapExtendPenalty;
        }
    }

    public RealignmentScore(final Kryo kryo, final Input input) {
        matchPenalty = input.readDouble();
        mismatchPenalty = input.readDouble();
        gapOpenPenalty = input.readDouble();
        gapExtendPenalty = input.readDouble();
        numberOfMatches = input.readInt();
        numberOfMismatches = input.readInt();
        numberOfIndels = input.readInt();
        indelLengthSum = input.readInt();
        strand = Strand.decode(input.readChar());
        value = input.readDouble();
    }

    public static RealignmentScore decode(final String str, final RealignmentScoreParameters parameters) {
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
        final Strand strand     = Strand.decode(parts[nextIdx]);
        return new RealignmentScore(parameters, matches, misMatches, indels, indelLenghts, strand);
    }

    private static RealignmentScore calculate(final RealignmentScoreParameters parameters,
                                              final int direction, final byte[] ref, final byte[] seq,
                                              final List<AlignmentInterval> intervals) {
        int totalIndels = 0;
        int totalMatches = 0;
        int totalMismatches = 0;
        int totalIndelLength = 0;
        if (intervals.isEmpty()) { // no intervals = all the sequence is a big insertion (score will be +Inf, Prob -Inf see #getPhredValue)
            totalIndelLength = seq.length;
            totalIndels = 1;
        } else {
            // deal with inter-interval events.
            for (final AlignmentInterval ai : intervals) {
                final int totalAligned = ai.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                        .filter(ce -> ce.getOperator().isAlignment())
                        .mapToInt(CigarElement::getLength).sum();
                final int misMatches = calculateMismatches(ref, seq, ai);
                final int matches = totalAligned - misMatches;
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
            // deal with intra-adjacent intervals (indel) events:
            for (int i = 1; i < intervals.size(); i++) {
                final AlignmentInterval ai = intervals.get(i);
                final AlignmentInterval prev = intervals.get(i - 1);
                final AlignmentInterval left = direction > 0 ? prev : ai; // left on the reference
                final AlignmentInterval right = direction > 0 ? ai : prev; // right on the reference
                final int refIndelLength = right.referenceSpan.getStart() - left.referenceSpan.getEnd() - 1;
                final int ctgIndelLength = prev.endInAssembledContig - ai.startInAssembledContig - 1;
                if (refIndelLength != 0 || ctgIndelLength != 0) {
                    final int indelLength = Math.max(Math.abs(refIndelLength) + Math.max(0, -ctgIndelLength), Math.abs(ctgIndelLength));
                    // The max(0, -ctgIndelLength) is to correct of short overlaps on the contig due to short "unclipped" soft-clips.
                    totalIndels++;
                    totalIndelLength += indelLength;
                }
            }
            // Treat clippings as indels:
            if (intervals.get(0).startInAssembledContig > 1) {
                final int indelLength = Math.min(intervals.get(0).startInAssembledContig - 1, intervals.get(0).referenceSpan.getStart());
                totalIndelLength += indelLength;
                totalIndels++;
            }
            if (intervals.get(intervals.size() - 1).endInAssembledContig < seq.length) {
                final int indelLength = seq.length - intervals.get(intervals.size() - 1).endInAssembledContig;
                totalIndelLength += indelLength;
                totalIndels++;
            }
        }
        return new RealignmentScore(parameters, totalMatches, totalMismatches, totalIndels, totalIndelLength,
                direction == 1 ? Strand.POSITIVE : Strand.NEGATIVE);
    }

    public static RealignmentScore calculate(final RealignmentScoreParameters parameters, final byte[] ref, final byte[] seq, final List<AlignmentInterval> intervals) {
        final List<AlignmentInterval> sortedIntervals = intervals == null ? Collections.emptyList()
                : intervals.stream()
                .sorted(Comparator.comparing(ai -> ai.startInAssembledContig))
                .collect(Collectors.toList());

        final int forwardAlignedBases = intervals.stream().filter(ai -> ai.forwardStrand && ai.mapQual >= parameters.minimumMappingQuality).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
        final int reverseAlignedBases = intervals.stream().filter(ai -> !ai.forwardStrand && ai.mapQual >= parameters.minimumMappingQuality).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
        if (forwardAlignedBases > reverseAlignedBases * parameters.dominantStrandAlignedBaseCountRatio ) {
            return calculate(parameters, 1, ref, seq, sortedIntervals.stream().filter(ai -> ai.forwardStrand).collect(Collectors.toList()));
        } else if (reverseAlignedBases > forwardAlignedBases * parameters.dominantStrandAlignedBaseCountRatio) {
            return calculate(parameters, -1, ref, seq, sortedIntervals.stream().filter(ai -> !ai.forwardStrand).collect(Collectors.toList()));
        } else {
            final int allForwardAlignedBases = intervals.stream().filter(ai -> ai.forwardStrand).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
            final int allReverseAlignedBases = intervals.stream().filter(ai -> !ai.forwardStrand).mapToInt(ai -> CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig)).sum();
            if (allForwardAlignedBases > allReverseAlignedBases * parameters.dominantStrandAlignedBaseCountRatio) {
                return calculate(parameters, 1, ref, seq, sortedIntervals.stream().filter(ai -> ai.forwardStrand).collect(Collectors.toList()));
            } else if (allReverseAlignedBases > allForwardAlignedBases * parameters.dominantStrandAlignedBaseCountRatio) {
                return calculate(parameters, -1, ref, seq, sortedIntervals.stream().filter(ai -> !ai.forwardStrand).collect(Collectors.toList()));
            } else {
                final RealignmentScore forward = calculate(parameters, 1, ref, seq, sortedIntervals.stream().filter(ai -> ai.forwardStrand).collect(Collectors.toList()));
                final RealignmentScore reverse = calculate(parameters, -1, ref, seq, sortedIntervals.stream().filter(ai -> !ai.forwardStrand).collect(Collectors.toList()));
                return forward.getPhredValue() <= reverse.getPhredValue() ? forward : reverse;
            }
        }
    }

    /**
     * Calculate the number of mismatching bases calls.
     * @param ref reference sequence.
     * @param refOffset first base in the reference aligned.
     * @param seq contig/query sequence.
     * @param seqOffset first base index in the ctg aligned.
     * @param length length of the target region.
     * @param forward whether the ctg is alinged on the forward (true) or reverse (false) strands.
     * @return 0 or greater but never more than {@code length}.
     */
     private static int calculateMismatches(final byte[] ref, final int refOffset,
                                            final byte[] seq, final int seqOffset,
                                            final int length, final boolean forward) {
         final int stop = refOffset + length;
         int result = 0;
         if (forward) {
             for (int i = refOffset, j = seqOffset; i < stop; i++, j++) {
                 if (!Nucleotide.intersect(ref[i], seq[j])) {
                     result++;
                 }
             }
         } else {
            for (int i = refOffset, j = seqOffset; i < stop; i++, j--) {
                if (!Nucleotide.intersect(Nucleotide.complement(ref[i]), seq[j])) {
                    result++;
                }
            }
        }
        return result;
    }

    /**
     * Calculate base call mismatches across an alignment interval.
     * @param ref the reference sequence.
     * @param seq the unclipped aligned sequence.
     * @param ai the alignment interval.
     * @return 0 or greater and never more than the number of bases aligned (in practice a number much smaller than that).
     */
    private static int calculateMismatches(final byte[] ref, final byte[] seq, final AlignmentInterval ai) {
        int refOffset = ai.referenceSpan.getStart() - 1;
        int direction = ai.forwardStrand ? 1 : -1;
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
                result += calculateMismatches(ref, refOffset, seq, seqOffset, element.getLength(), ai.forwardStrand);
            }
            if (element.getOperator().consumesReferenceBases()) refOffset += element.getLength();
            if (element.getOperator().consumesReadBases()) seqOffset += element.getLength() * direction;
        }
        return result;
    }

    public String toString() {
        return getPhredValue() + ":" + Utils.join(",", numberOfMatches, numberOfMismatches,
                numberOfIndels, indelLengthSum, strand.encodeAsString());
    }

    public double getPhredValue() {
        return value;
    }

    public double getLog10Prob() {
        return value * -.1;
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<RealignmentScore> {

        @Override
        public void write(Kryo kryo, Output output, RealignmentScore object) {
            output.writeDouble(object.matchPenalty);
            output.writeDouble(object.mismatchPenalty);
            output.writeDouble(object.gapOpenPenalty);
            output.writeDouble(object.gapExtendPenalty);
            output.writeInt(object.numberOfMatches);
            output.writeInt(object.numberOfMismatches);
            output.writeInt(object.numberOfIndels);
            output.writeInt(object.indelLengthSum);
            output.writeChar(object.strand.encodeAsChar());
            output.writeDouble(object.value);
        }

        @Override
        public RealignmentScore read(Kryo kryo, Input input, Class<RealignmentScore> type) {
            return new RealignmentScore(kryo, input);
        }
    }
}