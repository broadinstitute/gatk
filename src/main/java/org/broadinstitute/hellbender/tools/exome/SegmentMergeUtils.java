package org.broadinstitute.hellbender.tools.exome;

import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.exception.InsufficientDataException;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Helper/utility class for merging segments in a {@link SegmentedModel}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SegmentMergeUtils {
    private enum MergeDirection {
        LEFT, RIGHT, NONE
    }

    //random generator for breaking ties in small-segment merging
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));
    //Kolmogorov-Smirnov test for small-segment merging
    private static final int KOLMOGOROV_SMIRNOV_TEST_RANDOM_SEED = 42;
    private static final RandomGenerator rngKolmogorovSmirnovTest =
            RandomGeneratorFactory.createRandomGenerator(new Random(KOLMOGOROV_SMIRNOV_TEST_RANDOM_SEED));
    private static final KolmogorovSmirnovTest kolmogorovSmirnovTest =
            new KolmogorovSmirnovTest(rngKolmogorovSmirnovTest);

    private SegmentMergeUtils() {}

    /**
     * Returns a new segment specified by the outermost breakpoints of the two segments to be joined.
     * Segments must be on the same chromosome.
     * @param segment1  first segment to be joined
     * @param segment2  second segment to be joined
     * @return          a new segment constructed by joining the two segments
     */
    public static SimpleInterval mergeSegments(final SimpleInterval segment1, final SimpleInterval segment2) {
        if (!segment1.getContig().equals(segment2.getContig())) {
            throw new IllegalArgumentException(String.format("Cannot join segments " +
                            segment1.getContig() + ":%d-%d and " +
                            segment2.getContig() + ":%d-%d on different chromosomes.",
                    segment1.getStart(), segment1.getEnd(), segment2.getStart(), segment2.getEnd()));
        }
        final int start = Math.min(segment1.getStart(), segment2.getStart());
        final int end = Math.max(segment1.getEnd(), segment2.getEnd());

        return new SimpleInterval(segment1.getContig(), start, end);
    }

    /**
     * Returns the number of segments that contain a number of targets below a given threshold.
     * @param segments              list of segments
     * @param targets               targets to be segmented
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      number of segments containing a number of targets strictly less than the threshold
     */
    public static int countSmallSegments(final List<SimpleInterval> segments,
                                         final TargetCollection <TargetCoverage> targets,
                                         final int targetNumberThreshold) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        Utils.nonNull(targets, "The collection of targets cannot be null.");
        return (int) segments.stream()
                .filter(s -> targets.targetCount(s) < targetNumberThreshold)
                .count();
    }

    /**
     * Given a list of segments, returns a new, modifiable list of segments with the small segments (i.e., those
     * containing less than a specified number of targets) dropped.  Does not modify the input collections.
     * @param segments              list of segments (will not be modified)
     * @param targets               targets to be segmented
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      new list of segments with small segments dropped, never {@code null}
     */
    public static List<SimpleInterval> dropSmallSegments(final List<SimpleInterval> segments,
                                                         final TargetCollection <TargetCoverage> targets,
                                                         final int targetNumberThreshold) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        Utils.nonNull(targets, "The collection of targets cannot be null.");
        return segments.stream().filter(s -> targets.targetCount(s) >= targetNumberThreshold)
                .collect(Collectors.toList());
    }

    /**
     * Returns a new, modifiable list of segments with small segments (i.e., those containing less than a specified
     * number of targets) merged.  The list of segments is traversed from beginning to end.  Upon arrival at a small
     * segment, the segment is repeatedly merged with its closest neighboring segment until it is above threshold,
     * then the traversal resumes.  Does not modify the input collections.
     * @param segments              original list of segments (will not be modified)
     * @param genome                linear target-coverage and SNP-allele-count data to be segmented
     * @param targetNumberThreshold number of targets below which a segment is considered small
     * @return                      new list of segments with small segments merged, never {@code null}
     */
    public static List<SimpleInterval> mergeSmallSegments(final List<SimpleInterval> segments,
                                                          final Genome genome,
                                                          final int targetNumberThreshold) {
        Utils.nonNull(segments, "The list of segments cannot be null.");
        Utils.nonNull(genome.getTargets(), "The collection of targets contained in the genome cannot be null.");
        Utils.nonNull(genome.getSNPs(), "The collection of SNPs contained in the genome cannot be null.");
        final List<SimpleInterval> mergedSegments = new ArrayList<>(segments);
        int index = 0;
        while (index < mergedSegments.size()) {
            //if current segment is small, merge it with an adjacent segment
            if (genome.getTargets().targetCount(mergedSegments.get(index)) < targetNumberThreshold) {
                final MergeDirection direction = SmallSegments.calculateMergeDirection(mergedSegments, genome, index);
                if (direction == MergeDirection.LEFT) {
                    //current = merge(left, current), remove left, stay on current during next iteration
                    mergedSegments.set(index, mergeSegments(mergedSegments.get(index - 1), mergedSegments.get(index)));
                    mergedSegments.remove(index - 1);
                    index -= 2;
                } else if (direction == MergeDirection.RIGHT) {
                    //current = merge(current, right), remove right, stay on current during next iteration
                    mergedSegments.set(index, mergeSegments(mergedSegments.get(index), mergedSegments.get(index + 1)));
                    mergedSegments.remove(index + 1);
                    index--;
                }
            }
            index++; //if no merge performed, go to next segment during next iteration
        }
        //contigs containing only a single small segment do not get merged; drop these segments
        return dropSmallSegments(mergedSegments, genome.getTargets(), targetNumberThreshold);
    }

    /**
     * Contains private methods for small-segment merging.
     *
     * <p>
     *     The goal of small-segment merging is to merge all small segments (i.e., those with a number of targets below
     *     a given threshold) with adjacent segments.  This is done to match the output of circular binary segmentation,
     *     which does not generate segments with less than 3 targets.
     * </p>
     *
     * <p>
     *     For each small segment, the question we want to answer is: which of the two adjacent segments (i.e., those to
     *     the left or right of the small segment in the center) is "closer" to the small segment, and hence, should be
     *     merged with the small segment?  We would like to determine this based on the target coverages, SNP allele
     *     counts, and genomic locations of each of the three segments; note that each segment may be missing either
     *     coverages or allele counts (but not both).
     * </p>
     *
     * <p>
     *     However, this question is not statistically well defined, so we will use the following ad-hoc procedure:
     *     <ol>
     *         <li>
     *             We first check that all segments are on the same chromosome.  If only one of the adjacent segments is
     *             on the same chromosome, it is merged with by default.  If neither of the adjacent segments is on the
     *             same chromosome, no merge is performed.
     *         </li>
     *         <li>
     *             We then examine SNP allele counts (as opposed to target coverages---this is because the center
     *             segment has a small number of targets, by construction, and so the ability to determine which segment
     *             is closer using target coverages is inherently limited).  In particular, we examine the empirical
     *             alternate-allele fractions in each of the segments; these should have a bimodal (unimodal)
     *             distribution for unbalanced (balanced) segemnts.  The Kolmogorov-Smirnov test statistic, which
     *             measures the similarity of two data sets (it is zero (unity) for identical (dissimilar) data sets),
     *             is used to construct two distances between the alternate-allele fractions in the left and center
     *             segments and those in the right and center segments, respectively.  The distances are used to
     *             construct a pair of scores (with corrections for sizes of the data sets) for merging with the
     *             left and right segments, respectively; the sum of the scores is normalized to unity.  The adjacent
     *             segment with the higher score has alternate-allele fractions that are more similar, and hence it
     *             should be merged with the center small segment.
     *         </li>
     *         <li>
     *             However, if the Kolmogorov-Smirnov distances are not sufficiently dissimilar, if they are both close
     *             to unity (i.e., if the alternate-allele fractions in neither the left nor the right segment overlap
     *             significantly with those in the center segment), or if there are not enough SNPs to calculate the
     *             Kolmogorov-Smirnov distances (the Apache Commons implementation requires at least 2 data points
     *             in each data set), we instead use the empirical inverse minor-allele fractions (which, ideally,
     *             are proportional to total copy ratio and have a distribution that is practically unimodal) in each of
     *             the three segments.  Two distances between the two pairs of data sets are constructed using the
     *             Hodges-Lehmann estimator (which gives a measure of the difference in the location parameters of two
     *             data sets), and these distances are used to construct a pair of scores as above.
     *         </li>
     *         <li>
     *             If the scores generated from the SNP allele counts are too similar or if any of the segments is
     *             missing SNPs, then we attempt to use the target coverages instead.  Here, we simply use the
     *             Hodges-Lehmann estimator to construct two distances and a corresponding pair of scores as above.
     *         </li>
     *         <li>
     *             If the scores generated from the target coverages are too similar or if any of the segments is
     *             missing targets, we then simply use genomic distance (defined between adjacent breakpoints) to decide
     *             which adjacent segment is closer.  For consistency, we also convert the two genomic distances into a
     *             pair of scores that sum to unity.
     *         </li>
     *         <li>
     *             In the unlikely event that all of the above scores are equal, we randomly choose one of the adjacent
     *             segments and merge it with the center small segment.
     *         </li>
     *     </ol>
     * </p>
     */
    private static final class SmallSegments {
        //if larger of SNP scores is above SNP_SCORE_THRESHOLD, will use SNP scores as small-segment merge scores;
        //e.g., if SNP_SCORE_THRESHOLD = 0.55, then (0.6, 0.4) would be used, but (0.52, 0.48) would not.
        //SNP_SCORE_THRESHOLD should be in [0, 1], typically slightly above 0.5.
        private static final double SNP_SCORE_THRESHOLD = 0.55;

        //analogously, TARGET_SCORE_THRESHOLD serves the same purpose for target scores.
        private static final double TARGET_SCORE_THRESHOLD = 0.55;

        //SNP_KOLMOGOROV_SMIRNOV_DISTANCE_THRESHOLD is the Kolmogorov-Smirnov distance above which it is determined
        //there is no appreciable overlap between SNP alternate-allele--fraction data sets and the Hodges-Lehmann
        //estimator of inverse minor-allele fraction is used instead.
        //SNP_KOLMOGOROV_SMIRNOV_DISTANCE_THRESHOLD should be in [0, 1], typically slightly below 1.
        private static final double SNP_KOLMOGOROV_SMIRNOV_DISTANCE_THRESHOLD = 0.95;

        //SNP_KOLMOGOROV_SMIRNOV_DISTANCE_DIFFERENCE_THRESHOLD sets the difference in the Kolmogorov-Smirnov distance
        //below which SNP alternate-allele--fraction data sets are determined to be too similar and the Hodges-Lehmann
        //estimator of inverse minor-allele fraction is used instead.
        //SNP_KOLMOGOROV_SMIRNOV_DISTANCE_DIFFERENCE_THRESHOLD should be a small fraction in [0, 1].
        private static final double SNP_KOLMOGOROV_SMIRNOV_DISTANCE_DIFFERENCE_THRESHOLD = 0.05;

        //scores with an absolute difference less than DOUBLE_EQUALITY_EPSILON are taken to be equal
        private static final double DOUBLE_EQUALITY_EPSILON = 1E-6;

        /**
         * Returns genomic distance between two segments (specified by their indices in a list), which should be given
         * in order from left to right and are assumed to be non-overlapping.
         * This distance is infinite if the segments are on different chromosomes.
         * @param segments      list of segments
         * @param leftIndex     index of left segment
         * @param rightIndex    index of right segment
         * @return              genomic distance
         *                      (Double.POSITIVE_INFINITY if segments are on different chromosomes or
         *                      when checking distance to left (right) of first (last) segment)
         */
        private static double genomicDistance(final List<SimpleInterval> segments,
                                              final int leftIndex, final int rightIndex) {
            if (leftIndex < 0 || leftIndex >= segments.size() || rightIndex < 0 || rightIndex >= segments.size()) {
                return Double.POSITIVE_INFINITY;
            }
            final SimpleInterval leftSegment = segments.get(leftIndex);
            final SimpleInterval rightSegment = segments.get(rightIndex);
            if (leftSegment.getContig().equals(rightSegment.getContig())) {
                final double distance = rightSegment.getStart() - leftSegment.getEnd();
                if (distance < 0) {
                    throw new GATKException.ShouldNeverReachHereException("Negative genomic distance found during " +
                            "small-segment merging; segments are out of order or overlap.");
                }
                return distance;
            }
            return Double.POSITIVE_INFINITY;
        }

        /**
         * Given a segment specified by an index, returns a pair of scores for adjacent segments based on genomic
         * distance; except for edge cases, the sum of the scores will be unity.
         * Scores are Double.NEGATIVE_INFINITY for adjacent segments that are on different chromosomes.
         * @param segments  list of segments
         * @param index     index of the center segment to consider
         * @return          scores for adjacent segments based on genomic distance
         *                  (Double.NEGATIVE_INFINITY for adjacent segments that are on different chromosomes)
         */
        private static Pair<Double, Double> calculateGenomicDistanceScores(final List<SimpleInterval> segments,
                                                                           final int index) {
            final double leftDistance = genomicDistance(segments, index - 1, index);
            final double rightDistance = genomicDistance(segments, index, index + 1);

            //edge cases (divide-by-zero, segments on edges of chromosomes)
            if (leftDistance == 0. && rightDistance == 0.) {
                //should never reach this if segments are non-overlapping
                return Pair.of(0.5, 0.5);
            }
            if (leftDistance == Double.POSITIVE_INFINITY && rightDistance == Double.POSITIVE_INFINITY) {
                return Pair.of(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            }
            if (leftDistance == Double.POSITIVE_INFINITY) {
                return Pair.of(Double.NEGATIVE_INFINITY, 1.);
            }
            if (rightDistance == Double.POSITIVE_INFINITY) {
                return Pair.of(1., Double.NEGATIVE_INFINITY);
            }

            final double leftScore = 1. - leftDistance / (leftDistance + rightDistance);
            final double rightScore = 1. - rightDistance / (leftDistance + rightDistance);

            return Pair.of(leftScore, rightScore);
        }

        /**
         * Returns the effective square-root data-set size used in calculating the Kolmogorov-Smirnov distance.
         * @param data1     first data set
         * @param data2     second data set
         * @return          effective square-root data-set size used in calculating the Kolmogorov-Smirnov distance
         */
        private static double sqrtN(final double[] data1, final double[] data2) {
            final double n1 = data1.length;
            final double n2 = data2.length;
            return Math.sqrt((n1 * n2) / (n1 + n2));
        }

        /**
         * Returns a distance between two data sets given by the two-sample Kolmogorov-Smirnov test statistic.
         * @param data1     first data set
         * @param data2     second data set
         * @return          distance between data sets given by the two-sample Kolmogorov-Smirnov test statistic
         */
        private static double kolmogorovSmirnovDistance(final double[] data1, final double[] data2) {
            return kolmogorovSmirnovTest.kolmogorovSmirnovStatistic(data1, data2);
        }

        /**
         * Returns a list of the alternate-allele fractions for SNPs in a given segment.
         * @param segment   segment to consider
         * @param snps      SNP-allele-count data to be segmented
         * @return          list of alternate-allele fractions for SNPs in the segment
         */
        private static List<Double> calculateAAFs(final SimpleInterval segment,
                                                  final TargetCollection<AllelicCount> snps) {
            return snps.targets(segment).stream().map(AllelicCount::toAltAlleleFraction).collect(Collectors.toList());
        }

        /**
         * Returns a list of the inverse minor-allele fractions (which, ideally, are proportional to total copy ratio)
         * for SNPs in a given segment. Double.MIN_VALUE is added to the minor-allele fraction to avoid divide-by-zero
         * cases.
         * @param segment   segment to consider
         * @param snps      SNP-allele-count data to be segmented
         * @return          list of inverse minor-allele fractions for SNPs in the segment
         */
        private static List<Double> calculateInverseMAFs(final SimpleInterval segment,
                                                         final TargetCollection<AllelicCount> snps) {
            return snps.targets(segment).stream().map(a -> 1. / (Double.MIN_VALUE + a.toMinorAlleleFraction()))
                    .collect(Collectors.toList());
        }

        /**
         * Calculates the distance between two data sets based on the Hodges-Lehmann estimator.
         * @param data1     first data set
         * @param data2     second data set
         * @return          distance between data sets based on the Hodges-Lehmann estimator
         */
        private static double hodgesLehmannDistance(final double[] data1, final double[] data2) {
            double[] differences = new double[data1.length * data2.length];
            for (int i = 0; i < data1.length; i++) {
                for (int j = 0; j < data2.length; j++) {
                    differences[i * data2.length + j] = data1[i] - data2[j];
                }
            }
            return new Median().evaluate(differences);
        }

        /**
         * Given three data sets (corresponding to left, center, and right segments), returns a pair of scores based on
         * the Hodges-Lehmann estimator between (left, center) and (center, right); the sum of the scores will be unity.
         * @param leftData      data set for left segment
         * @param centerData    data set for center segment
         * @param rightData     data set for right segment
         * @return              pair of scores based on the Hodges-Lehmann estimator
         */
        private static Pair<Double, Double> calculateHodgesLehmannScores(final double[] leftData,
                                                                         final double[] centerData,
                                                                         final double[] rightData) {
            final double leftDistance = hodgesLehmannDistance(leftData, centerData);
            final double rightDistance = hodgesLehmannDistance(centerData, rightData);

            if (leftDistance == 0. && rightDistance == 0.) {
                return Pair.of(0.5, 0.5);
            }

            //if center segment is above or below both left and right segments,
            //assign score 1 to the closer segment and 0 to the other
            if (leftDistance * rightDistance < 0) {
                return Math.abs(leftDistance) < Math.abs(rightDistance) ? Pair.of(1., 0.) : Pair.of(0., 1.);
            }
            return Pair.of(1. - Math.abs(leftDistance / (leftDistance + rightDistance)),
                           1. - Math.abs(rightDistance / (leftDistance + rightDistance)));
        }

        /**
         * Given a segment specified by an index, returns a pair of scores for adjacent segments based on 2-sample
         * Kolmogorov-Smirnov tests of the observed alternate-allele fractions; except for edge cases, the sum of the
         * scores will be unity. All segments are assumed to be on the same chromosome.
         * If any of the three segments is missing SNPs, both scores are Double.NEGATIVE_INFINITY.
         * @param segments  list of segments
         * @param snps      SNP-allele-count data to be segmented
         * @param index     index of the center segment to consider
         * @return          scores for adjacent segments based on 2-sample Kolmogorov-Smirnov tests
         */
        private static Pair<Double, Double> calculateSNPScores(final List<SimpleInterval> segments,
                                                               final TargetCollection<AllelicCount> snps,
                                                               final int index) {
            final SimpleInterval leftSegment = segments.get(index - 1);
            final SimpleInterval centerSegment = segments.get(index);
            final SimpleInterval rightSegment = segments.get(index + 1);

            //check if any segment is missing SNPs
            if (snps.targetCount(leftSegment) == 0 || snps.targetCount(centerSegment) == 0 ||
                    snps.targetCount(rightSegment) == 0) {
                return Pair.of(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            }

            //calculate Kolmogorov-Smirnov distances based on alternate-allele-fractions in each segment
            final double[] leftAAFs = Doubles.toArray(calculateAAFs(leftSegment, snps));
            final double[] centerAAFs = Doubles.toArray(calculateAAFs(centerSegment, snps));
            final double[] rightAAFs = Doubles.toArray(calculateAAFs(rightSegment, snps));
            try {
                //if not enough alternate-allele fractions in any segment to use Apache Commons implementation of
                //kolmogorovSmirnovStatistic, kolmogorovSmirnovDistance will throw an unchecked
                //InsufficientDataException and the code after the catch block will be executed
                final double leftKSDistance = kolmogorovSmirnovDistance(leftAAFs, centerAAFs);
                final double rightKSDistance = kolmogorovSmirnovDistance(centerAAFs, rightAAFs);

                //edge case (divide-by-zero)
                if (leftKSDistance == 0. && rightKSDistance == 0.) {
                    //this is unlikely to occur using the Apache Commons implementation of kolmogorovSmirnovStatistic
                    return Pair.of(0.5, 0.5);
                }

                //if alternate-allele fractions in all segments are too similar or do not overlap appreciably,
                //will not return anything here and the code after the catch block will be executed
                if (Math.abs(leftKSDistance - rightKSDistance) > SNP_KOLMOGOROV_SMIRNOV_DISTANCE_DIFFERENCE_THRESHOLD &&
                        (leftKSDistance < SNP_KOLMOGOROV_SMIRNOV_DISTANCE_THRESHOLD ||
                                rightKSDistance < SNP_KOLMOGOROV_SMIRNOV_DISTANCE_THRESHOLD)) {
                    final double leftSqrtN = sqrtN(leftAAFs, centerAAFs);
                    final double rightSqrtN = sqrtN(centerAAFs, rightAAFs);

                    return Pair.of(
                            1. - leftKSDistance * leftSqrtN /
                                    (leftKSDistance * leftSqrtN + rightKSDistance * rightSqrtN),
                            1. - rightKSDistance * rightSqrtN /
                                    (leftKSDistance * leftSqrtN + rightKSDistance * rightSqrtN));
                }
            } catch (final InsufficientDataException e) {
                //do nothing here, continue below to use Hodges-Lehmann scores instead
            }
            //use Hodges-Lehmann scores computed using inverse minor-allele fractions
            //(which, ideally, are proportional to total copy ratio)
            //will be executed after an InsufficientDataException or alternate-allele fractions in all segments
            //are too similar or do not overlap appreciably
            final double[] leftInverseMAFs = Doubles.toArray(calculateInverseMAFs(leftSegment, snps));
            final double[] centerInverseMAFs = Doubles.toArray(calculateInverseMAFs(centerSegment, snps));
            final double[] rightInverseMAFs = Doubles.toArray(calculateInverseMAFs(rightSegment, snps));

            return calculateHodgesLehmannScores(leftInverseMAFs, centerInverseMAFs, rightInverseMAFs);
        }

        /**
         * Returns a list of the coverages for targets in a given segment.
         * @param segment   segment to consider
         * @param targets   linear target-coverage data to be segmented
         * @return          list of the coverages for targets in the segment
         */
        private static List<Double> makeCoverageList(final SimpleInterval segment,
                                                     final TargetCollection<TargetCoverage> targets) {
            return targets.targets(segment).stream().map(t -> t.getCoverage()).collect(Collectors.toList());
        }

        /**
         * Given a segment specified by an index, returns a pair of scores for adjacent segments based on the
         * Hodges-Lehmann estimators between the observed linear target coverages; except for edge cases,
         * the sum of the scores will be unity. All segments are assumed to be on the same chromosome.
         * If any of the three segments is missing targets, both scores are Double.NEGATIVE_INFINITY.
         * @param segments  list of segments
         * @param targets   linear target-coverage data to be segmented
         * @param index     index of the center segment to consider
         * @return          scores for adjacent segments based on Hodges-Lehmann estimators
         */
        private static Pair<Double, Double> calculateTargetScores(final List<SimpleInterval> segments,
                                                                  final TargetCollection<TargetCoverage> targets,
                                                                  final int index) {
            final SimpleInterval leftSegment = segments.get(index - 1);
            final SimpleInterval centerSegment = segments.get(index);
            final SimpleInterval rightSegment = segments.get(index + 1);

            //check if any segment is missing targets
            if (targets.targetCount(leftSegment) == 0 || targets.targetCount(centerSegment) == 0 ||
                    targets.targetCount(rightSegment) == 0) {
                return Pair.of(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
            }

            final double[] leftCoverages = Doubles.toArray(makeCoverageList(leftSegment, targets));
            final double[] centerCoverages = Doubles.toArray(makeCoverageList(centerSegment, targets));
            final double[] rightCoverages = Doubles.toArray(makeCoverageList(rightSegment, targets));

            return calculateHodgesLehmannScores(leftCoverages, centerCoverages, rightCoverages);
        }

        /**
         * Given a segment specified by an index, returns the small-segment merge scores for merging with adjacent
         * segments.
         * @param segments  list of segments
         * @param genome    linear target-coverage and SNP-allele-count data to be segmented
         * @param index     index of the center segment to consider
         * @return          pair of small-segment merge scores for merging with adjacent segments
         *                  (Double.NEGATIVE_INFINITY for adjacent segments that are on different chromosomes)
         */
        private static Pair<Double, Double> calculateScores(final List<SimpleInterval> segments,
                                                            final Genome genome,
                                                            final int index) {
            final Pair<Double, Double> genomicDistanceScores = calculateGenomicDistanceScores(segments, index);

            //if either adjacent segment on different chromosome or out of range, no need to do non-parametric tests
            if (genomicDistanceScores.getLeft() == Double.NEGATIVE_INFINITY ||
                    genomicDistanceScores.getRight() == Double.NEGATIVE_INFINITY) {
                return genomicDistanceScores;
            }

            //first, try to use SNP scores based on 2-sample Kolmogorov-Smirnov test (or Hodges-Lehmann estimator, if no
            //appreciable overlap between samples, Kolmogorov-Smirnov distances are too similar, or not enough SNPs
            //to use Apache Commons implementation of kolmogorovSmirnovStatistic)
            final TargetCollection<AllelicCount> snps = genome.getSNPs();
            final Pair<Double, Double> snpScores = calculateSNPScores(segments, snps, index);
            if (Math.max(snpScores.getLeft(), snpScores.getRight()) > SNP_SCORE_THRESHOLD) {
                return snpScores;
            }

            //if any of the three segments is missing SNPs or the SNP scores are too similar,
            //try to use target scores based on Hodges-Lehmann estimator
            final TargetCollection<TargetCoverage> targets = genome.getTargets();
            final Pair<Double, Double> targetScores = calculateTargetScores(segments, targets, index);
            if (Math.max(targetScores.getLeft(), targetScores.getRight()) > TARGET_SCORE_THRESHOLD) {
                return targetScores;
            }

            //if any of the three segments is missing targets or the target scores are too similar,
            //use genomic-distance scores
            return genomicDistanceScores;
        }

        /**
         * Given a segment specified by an index, returns the direction of the adjacent segment with which it should be
         * merged (i.e., the adjacent segment with higher small-segment merge score given by
         * {@link SegmentMergeUtils.SmallSegments#calculateScores}).
         * If both adjacent segments are on different chromosomes than the specified segment, returns
         * MergeDirection.NONE.
         * @param segments  list of segments
         * @param genome    linear target-coverage and SNP-allele-count data to be segmented
         * @param index     index of the center segment to consider
         * @return          direction of the adjacent segment with higher small-segment merge score
         *                  (MergeDirection.NONE if both adjacent segments are on different chromosomes)
         */
        private static MergeDirection calculateMergeDirection(final List<SimpleInterval> segments,
                                                              final Genome genome,
                                                              final int index) {
            Utils.validIndex(index, segments.size());
            final Pair<Double, Double> scores = calculateScores(segments, genome, index);
            if (scores.getLeft() > scores.getRight()) {
                return MergeDirection.LEFT;
            }
            if (scores.getRight() > scores.getLeft()) {
                return MergeDirection.RIGHT;
            }
            if (scores.getLeft() == Double.NEGATIVE_INFINITY && scores.getRight() == Double.NEGATIVE_INFINITY) {
                return MergeDirection.NONE;
            }
            if (Math.abs(scores.getLeft() - scores.getRight()) <= DOUBLE_EQUALITY_EPSILON) {
                //in the unlikely event that all scores are equal, merge randomly
                final boolean isLeftScoreHigher = rng.nextBoolean();
                return isLeftScoreHigher ? MergeDirection.LEFT : MergeDirection.RIGHT;
            }
            throw new GATKException.ShouldNeverReachHereException("Something went wrong during small-segment merging.");
        }
    }
}