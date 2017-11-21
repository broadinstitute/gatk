package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.report.GATKReportColumnFormat;
import scala.Tuple2;
import scala.Tuple3;

import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.function.ToDoubleBiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 6/9/17.
 */
public class TemplateMappingInformation implements Serializable {

    private static final long serialVersionUID = 1L;

    private static final Comparator<List<AlignmentInterval>> LEFT_RIGHT_ALIGNMENT_COMPARATOR =
            Comparator.comparingInt(TemplateMappingInformation::unclippedStart)
                    .thenComparingInt(TemplateMappingInformation::unclippedEnd)
                    .thenComparingInt(TemplateMappingInformation::clippedStart)
                    .thenComparingInt(TemplateMappingInformation::clippedEnd);

    public String name;
    public final ReadPairOrientation pairOrientation;
    public OptionalDouble firstAlignmentScore = OptionalDouble.empty();
    public OptionalDouble secondAlignmentScore = OptionalDouble.empty();
    public List<AlignmentInterval> firstAlignmentIntervals;
    public List<AlignmentInterval> secondAlignmentIntervals;
    public final OptionalInt insertSize;
    public final int minCoordinate;
    public final int maxCoordinate;

    public static TemplateMappingInformation fromAlignments(final List<AlignmentInterval> firstIntervals, final int firstLength,
                                                      final List<AlignmentInterval> secondIntervals, final int secondLength) {
        Utils.nonNull(firstIntervals);
        Utils.nonNull(secondIntervals);
        ParamUtils.isPositive(firstLength, "the first length must be 1 or greater");
        ParamUtils.isPositive(secondLength, "the second length must be 1 or greater");
        if (firstIntervals.isEmpty() && secondIntervals.isEmpty()) {
            return new TemplateMappingInformation();
        } else if (secondIntervals.isEmpty()) {
            return new TemplateMappingInformation(score(firstIntervals, firstLength), unclippedStart(firstIntervals), unclippedEnd(firstIntervals), firstIntervals, true);
        } else if (firstIntervals.isEmpty()) {
            return new TemplateMappingInformation(score(secondIntervals, secondLength), unclippedStart(secondIntervals), unclippedEnd(secondIntervals), secondIntervals, false);
        } else {
            final Pair<List<AlignmentInterval>, List<AlignmentInterval>> sortedAlignments
                    = sortLeftRightAlignments(firstIntervals, secondIntervals);
            final Pair<Strand, Strand> strands = new ImmutablePair<>(
                    strand(sortedAlignments.getLeft()), strand(sortedAlignments.getRight()));
            final ReadPairOrientation orientation = ReadPairOrientation.fromStrands(strands.getLeft(), strands.getRight());
            if (orientation.isProper()) {
                return new TemplateMappingInformation(score(firstIntervals, firstLength),
                                                      score(secondIntervals, secondLength), unclippedStart(sortedAlignments.getLeft()), unclippedEnd(sortedAlignments.getRight()),
                        firstIntervals, secondIntervals,
                        unclippedEnd(sortedAlignments.getRight()) - unclippedStart(sortedAlignments.getLeft()));
            } else {
                return new TemplateMappingInformation(score(firstIntervals, firstLength),
                        score(secondIntervals, secondLength), unclippedStart(sortedAlignments.getLeft()), unclippedEnd(sortedAlignments.getRight()), firstIntervals, secondIntervals, orientation);
            }
        }
    }

    public TemplateMappingInformation(final double firstAlignment, final double secondAlignment, final int minCoordinate, final int maxCoordinate, final List<AlignmentInterval> firstAlignmentIntervals, final List<AlignmentInterval> secondAlignmentIntervals, final ReadPairOrientation orientation) {
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

    public TemplateMappingInformation(final double alignment, final int minCoordinate, final int maxCoordinate, final List<AlignmentInterval> intervals, final boolean isFirst) {
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

    public TemplateMappingInformation(final double firstAlignment, final double secondAlignment, final int minCoordinate, final int maxCoordinate,  final List<AlignmentInterval> firstIntervals, final List<AlignmentInterval> secondIntervals, final int insertSize) {
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

    public TemplateMappingInformation() {
        firstAlignmentScore = OptionalDouble.empty();
        secondAlignmentScore = OptionalDouble.empty();
        insertSize = OptionalInt.empty();
        pairOrientation = ReadPairOrientation.XX;
        minCoordinate =  Integer.MAX_VALUE;
        maxCoordinate = Integer.MIN_VALUE;
    }

    private static double score(final List<AlignmentInterval> intervals, final int length) {
        return intervals.isEmpty() ? Double.NaN : AlignmentScore.calculate(length, intervals).getValue();

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

    public boolean crossesBreakPoint(final int[] breakPoints, final int fragment) {
        ParamUtils.inRange(fragment, 0, 1, "fragment out of range");
        final List<AlignmentInterval> intervals = fragment == 0 ? firstAlignmentIntervals : secondAlignmentIntervals;
        if (intervals == null || intervals.isEmpty()) return false;
        int minCoordinate = Integer.MAX_VALUE;
        int maxCoordinate = Integer.MIN_VALUE;
        for (final AlignmentInterval ai : intervals) {
            final SimpleInterval refSpan = ai.referenceSpan;
            final Cigar cigar = ai.cigarAlongReference();
            final int minRef = refSpan.getStart() - CigarUtils.countLeftClippedBases(cigar);
            final int maxRef = refSpan.getEnd() + CigarUtils.countRightClippedBases(cigar);
            if (minRef < minCoordinate) minCoordinate = minRef;
            if (maxRef > maxCoordinate) maxCoordinate = maxRef;
        }
        for (final int breakPoint : breakPoints) {
            if (breakPoint >= minCoordinate && breakPoint <= maxCoordinate) {
                return true;
            }
        }
        return false;
    }

}
