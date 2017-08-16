package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;

/**
 * Contains information as to how a template maps against a
 * haplotype or contig.
 */
public final class TemplateMapping implements Serializable {

    private static final long serialVersionUID = 1L;

    private static final Comparator<List<AlignmentInterval>> LEFT_RIGHT_ALIGNMENT_COMPARATOR =
            Comparator.comparingInt(TemplateMapping::unclippedStart)
                    .thenComparingInt(TemplateMapping::unclippedEnd)
                    .thenComparingInt(TemplateMapping::clippedStart)
                    .thenComparingInt(TemplateMapping::clippedEnd);

    public String name;
    public final ReadPairOrientation pairOrientation;
    public OptionalDouble firstAlignmentScore = OptionalDouble.empty();
    public OptionalDouble secondAlignmentScore = OptionalDouble.empty();
    public List<AlignmentInterval> firstAlignmentIntervals;
    public List<AlignmentInterval> secondAlignmentIntervals;
    public final OptionalInt insertSize;
    public final int minCoordinate;
    public final int maxCoordinate;

    public static TemplateMapping fromAlignments(final RealignmentScoreParameters realignmentScoreParameters, final SVHaplotype haplotype,
                                                 final byte[] firstBases, final List<AlignmentInterval> firstIntervals,
                                                 final byte[] secondBases, final List<AlignmentInterval> secondIntervals) {
        Utils.nonNull(firstIntervals);
        Utils.nonNull(secondIntervals);
        ParamUtils.isPositive(firstBases.length, "the first length must be 1 or greater");
        ParamUtils.isPositive(secondBases.length, "the second length must be 1 or greater");
        if (firstIntervals.isEmpty() && secondIntervals.isEmpty()) {
            return new TemplateMapping();
        } else if (secondIntervals.isEmpty()) {
            return new TemplateMapping(score(realignmentScoreParameters, haplotype, firstBases, firstIntervals), clippedStart(firstIntervals), clippedEnd(firstIntervals), firstIntervals, true);
        } else if (firstIntervals.isEmpty()) {
            return new TemplateMapping(score(realignmentScoreParameters, haplotype, secondBases, secondIntervals), clippedStart(secondIntervals), clippedEnd(secondIntervals), secondIntervals, false);
        } else {
            final Pair<List<AlignmentInterval>, List<AlignmentInterval>> sortedAlignments
                    = sortLeftRightAlignments(firstIntervals, secondIntervals);
            final Pair<Strand, Strand> strands = new ImmutablePair<>(
                    strand(sortedAlignments.getLeft()), strand(sortedAlignments.getRight()));
            final ReadPairOrientation orientation = ReadPairOrientation.of(strands.getLeft(), strands.getRight());
            if (orientation.isProper()) {
                return new TemplateMapping(score(realignmentScoreParameters, haplotype, firstBases, firstIntervals),
                                                      score(realignmentScoreParameters, haplotype, secondBases, secondIntervals), clippedStart(sortedAlignments.getLeft()), clippedEnd(sortedAlignments.getRight()),
                        firstIntervals, secondIntervals,
                        unclippedEnd(sortedAlignments.getRight()) - unclippedStart(sortedAlignments.getLeft()));
            } else {
                return new TemplateMapping(score(realignmentScoreParameters, haplotype, firstBases, firstIntervals),
                        score(realignmentScoreParameters, haplotype, secondBases, secondIntervals), clippedStart(sortedAlignments.getLeft()), clippedEnd(sortedAlignments.getRight()), firstIntervals, secondIntervals, orientation);
            }
        }
    }

    public TemplateMapping(final double firstAlignment, final double secondAlignment, final int minCoordinate, final int maxCoordinate, final List<AlignmentInterval> firstAlignmentIntervals, final List<AlignmentInterval> secondAlignmentIntervals, final ReadPairOrientation orientation) {
        Utils.nonNull(orientation);
        if (orientation.isProper()) {
            throw new IllegalArgumentException("you cannot create a mapping information object with proper orientation without indicating the insert size");
        }
        firstAlignmentScore = OptionalDouble.of(firstAlignment);
        secondAlignmentScore = OptionalDouble.of(secondAlignment);
        this.firstAlignmentIntervals = firstAlignmentIntervals;
        this.secondAlignmentIntervals = secondAlignmentIntervals;
        pairOrientation = orientation;
        insertSize = OptionalInt.empty();
        this.minCoordinate = minCoordinate;
        this.maxCoordinate = maxCoordinate;
    }

    public TemplateMapping(final double alignment, final int minCoordinate, final int maxCoordinate, final List<AlignmentInterval> intervals, final boolean isFirst) {
        if (isFirst) {
            firstAlignmentScore = OptionalDouble.of(alignment);
            firstAlignmentIntervals = intervals;
        } else {
            secondAlignmentScore = OptionalDouble.of(alignment);
            secondAlignmentIntervals = intervals;
        }
        pairOrientation = ReadPairOrientation.XX;
        insertSize = OptionalInt.empty();
        this.minCoordinate = minCoordinate;
        this.maxCoordinate = maxCoordinate;
    }

    public TemplateMapping(final double firstAlignment, final double secondAlignment, final int minCoordinate, final int maxCoordinate, final List<AlignmentInterval> firstIntervals, final List<AlignmentInterval> secondIntervals, final int insertSize) {
        Utils.nonNull(firstAlignment);
        Utils.nonNull(secondAlignment);
        if (insertSize < 1) {
            throw new IllegalArgumentException("the input insert size cannot be less than 1");
        }
        firstAlignmentScore = OptionalDouble.of(firstAlignment);
        secondAlignmentScore = OptionalDouble.of(secondAlignment);
        firstAlignmentIntervals = firstIntervals;
        secondAlignmentIntervals = secondIntervals;
        this.insertSize = OptionalInt.of(insertSize);
        pairOrientation = ReadPairOrientation.PROPER;
        this.minCoordinate = minCoordinate;
        this.maxCoordinate = maxCoordinate;
    }

    public TemplateMapping() {
        firstAlignmentScore = OptionalDouble.empty();
        secondAlignmentScore = OptionalDouble.empty();
        insertSize = OptionalInt.empty();
        pairOrientation = ReadPairOrientation.XX;
        minCoordinate =  Integer.MAX_VALUE;
        maxCoordinate = Integer.MIN_VALUE;
    }

    private static double score(final RealignmentScoreParameters parameters, final SVHaplotype haplotype, final byte[] seq, final List<AlignmentInterval> intervals) {
        return intervals.isEmpty() ? Double.NaN : RealignmentScore.calculate(parameters, haplotype.getBases(), seq, intervals).getLog10Prob();

    }

    private static Pair<List<AlignmentInterval>, List<AlignmentInterval>> sortLeftRightAlignments(
            final List<AlignmentInterval> first, final List<AlignmentInterval> second) {
        final int cmp = LEFT_RIGHT_ALIGNMENT_COMPARATOR.compare(first, second);
        return (cmp <= 0) ? new ImmutablePair<>(first, second)
                : new ImmutablePair<>(second, first);
    }

    private static Strand strand(final List<AlignmentInterval> firstIntervals) {
        final int mappedBasesOrientation = firstIntervals.stream()
                .mapToInt(ai -> (ai.forwardStrand ? 1 : -1) * CigarUtils.countAlignedBases(ai.cigarAlong5to3DirectionOfContig))
                .sum();
        if (mappedBasesOrientation != 0) {
            return mappedBasesOrientation < 0 ? Strand.NEGATIVE : Strand.POSITIVE;
        } else { // tie-break:
            final Comparator<AlignmentInterval> comparator0 = Comparator.comparingInt(a -> CigarUtils.countAlignedBases(a.cigarAlong5to3DirectionOfContig));
            final Comparator<AlignmentInterval> comparator1 = comparator0.thenComparingInt(a -> a.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                    .filter(ce -> ce.getOperator() == CigarOperator.I)
                    .mapToInt(CigarElement::getLength)
                    .sum());
            final Comparator<AlignmentInterval> comparator = comparator1.thenComparingInt(a -> a.startInAssembledContig).reversed();

            final boolean forwardStrand = firstIntervals.stream().sorted(comparator).findFirst().get().forwardStrand;
            return forwardStrand ? Strand.POSITIVE : Strand.NEGATIVE;
        }
    }

    private static int unclippedStart(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getStart() - (ai.forwardStrand ? CigarUtils.countLeftClippedBases(ai.cigarAlong5to3DirectionOfContig)
                        : CigarUtils.countRightClippedBases(ai.cigarAlong5to3DirectionOfContig)))
                .min().getAsInt();
    }

    private static int clippedStart(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getStart())
                .min().getAsInt();
    }

    private static int unclippedEnd(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getEnd() + (ai.forwardStrand ? CigarUtils.countRightClippedBases(ai.cigarAlong5to3DirectionOfContig)
                        : CigarUtils.countLeftClippedBases(ai.cigarAlong5to3DirectionOfContig)))
                .max().getAsInt();
    }

    private static int clippedEnd(final List<AlignmentInterval> intervals) {
        return intervals.stream()
                .mapToInt(ai -> ai.referenceSpan.getEnd())
                .max().getAsInt();
    }

    public boolean crossesBreakPoint(final int[] breakPoints) {
        if (minCoordinate > maxCoordinate) return false;
        for (final int breakPoint : breakPoints) {
            if (minCoordinate <= breakPoint && maxCoordinate >= breakPoint) {
                return true;
            }
        }
        return false;
    }

    @SuppressWarnings("unused") // may be useful in the future.
    public boolean crossesBreakPointCountingClippedBases(final int[] breakPoints) {
        int minCoordinate = Integer.MAX_VALUE;
        int maxCoordinate = 0;
        for (final AlignmentInterval interval : Utils.concat(firstAlignmentIntervals, secondAlignmentIntervals)) {
            final SimpleInterval referenceSpan = interval.referenceSpan;
            final int minRefPos = referenceSpan.getStart() - CigarUtils.countLeftClippedBases(interval.cigarAlongReference());
            final int maxRefPos = referenceSpan.getStart() - CigarUtils.countRightClippedBases(interval.cigarAlongReference());
            if (minRefPos < minCoordinate) minCoordinate = minRefPos;
            if (maxRefPos > maxCoordinate) maxCoordinate = maxRefPos;
        }
        for (final int breakPoint : breakPoints) {
            if (breakPoint >= minCoordinate && breakPoint <= maxCoordinate) return true;
        }
        return false;
    }
}
