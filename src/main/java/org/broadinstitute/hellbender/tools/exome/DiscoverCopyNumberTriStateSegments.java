package org.broadinstitute.hellbender.tools.exome;

import com.google.common.collect.Iterators;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModel;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModelArgumentCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Discover potential copy-number segments in input read-counts as the most likely segment sequences
 * based on {@link CopyNumberTriStateHiddenMarkovModel} HMM model.
 *
 * <p>
 *     You normally want to run this tool on values normalized using {@link SubtractCoverageComponents}
 *     which should get rid of systematic biases due to sequencing and capture technology
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        programGroup = CopyNumberProgramGroup.class,
        summary = "Finds possible locations for rare copy number variation events in germline samples using a HMM",
        oneLineSummary = "Discovers possible location of copy number variation"
)
public final class DiscoverCopyNumberTriStateSegments extends CommandLineProgram {

    public static final String ZSCORE_DIMENSION_FULL_NAME = "standardizeBy";
    public static final String ZSCORE_DIMENSION_SHORT_NAME = "standardizeBy";

    private static final double INV_LN_10 = 1.0 / Math.log(10);

    /**
     * Maximum reportable output quality score; higher quality scores will
     * capped to this value.
     */
    public static final double MAX_QUAL_SCORE = 199.99999999999;

    @ArgumentCollection
    protected CopyNumberTriStateHiddenMarkovModelArgumentCollection modelArguments =
            new CopyNumberTriStateHiddenMarkovModelArgumentCollection();

    /**
     * Z-scores transformation dimensions.
     */
    public enum ZScoreDimension {

        /**
         * Input coverage must be transformed in z-scores sample by sample.
         */
        SAMPLE,

        /**
         * Input coverage must be transformed in z-scores target by target.
         */
        TARGET,

        /**
         * No transformation must be performed.
         */
        NONE
    }

    @Argument(
            doc = "In what dimension to apply the zscore transformation",
            fullName = ZSCORE_DIMENSION_FULL_NAME,
            shortName = ZSCORE_DIMENSION_SHORT_NAME,
            optional = true
    )
    protected ZScoreDimension zscoreDimension = ZScoreDimension.SAMPLE;

    @Argument(
            doc = "Input counts",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputFile;

    @Argument(
            doc = "Output counts",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection(() -> inputFile);

    @Override
    protected Object doWork() {
        final CopyNumberTriStateHiddenMarkovModel model = modelArguments.createModel();
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);
        final ReadCountCollection inputCounts = readAndSortCountsByTargetCoordinates(targets);
        applyZScoreTransformation(inputCounts);
        checkForMissingTargetsInInputCounts(targets, inputCounts);

        final Map<String, List<CopyNumberTriStateSegment>> allSegmentsBySampleName =
                calculateBestPathSegments(model, inputCounts);

        final Iterator<CopyNumberTriStateSegmentRecord> allSegmentRecordsSortedByCoordinates =
                composeATargetSortedSegmentRecordIterator(targets, allSegmentsBySampleName);

        writeSegmentsInOutputFile(allSegmentRecordsSortedByCoordinates);

        return "SUCCESS";
    }

    /**
     * Transform read counts to z-scores.
     * @param inputCounts the input read-counts, modified in-situ.
     */
    private void applyZScoreTransformation(final ReadCountCollection inputCounts) {
        final RealMatrix counts = inputCounts.counts();
        switch (zscoreDimension) {
            case SAMPLE:
                standardizeBySample(counts);
                break;
            case TARGET:
                standardizeByTarget(counts);
                break;
        }
    }

    private void standardizeByTarget(final RealMatrix counts) {
        final double[] rowMeans = GATKProtectedMathUtils.rowMeans(counts);
        final double[] rowStdDevs = GATKProtectedMathUtils.rowStdDevs(counts);

        counts.walkInColumnOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return (value - rowMeans[row]) / rowStdDevs[row];
            }
        });
    }

    private void standardizeBySample(final RealMatrix counts) {
        final double[] columnMeans = GATKProtectedMathUtils.columnMeans(counts);
        final double[] columnStdDev = GATKProtectedMathUtils.columnStdDevs(counts);

        counts.walkInColumnOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return (value - columnMeans[column]) / columnStdDev[column];
            }
        });
    }

    /**
     * Writes sorted segments into the output file provided by the user.
     * @param allSegmentsSorted iterator that must produce segment records in order.
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CouldNotCreateOutputFile if there is any
     * issue creating or writing into the output file provided by the user.
     */
    private void writeSegmentsInOutputFile(final Iterator<CopyNumberTriStateSegmentRecord> allSegmentsSorted) {
        try (final CopyNumberTriStateSegmentRecordWriter writer = new CopyNumberTriStateSegmentRecordWriter(outputFile)) {
            long count = 0;
            final Set<String> samples = new HashSet<>();
            while (allSegmentsSorted.hasNext()) {
                final CopyNumberTriStateSegmentRecord record = allSegmentsSorted.next();
                samples.add(record.getSampleName());
                writer.writeRecord(record);
                count++;
            }
            logger.info(String.format("Found a total of %d segments across %d samples", count, samples.size()));

        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problems writing into the output", ex);
        }
    }

    /**
     * Create a merging iterator that sorts all the segments across samples by
     * their coordinates.
     * <p>
     *     The order between contigs is determined by the input target collection.
     * </p>
     * <p>
     *     The input map of sample name to segment list is assumed to contain all
     *     segments within each value list already in the correct order.
     * </p>
     *
     * @param targets the target collection.
     * @param allSegments map with all segment lists coming from different samples.
     * @return never null.
     */
    private Iterator<CopyNumberTriStateSegmentRecord> composeATargetSortedSegmentRecordIterator(final TargetCollection<Target> targets, final Map<String, List<CopyNumberTriStateSegment>> allSegments) {
        final List<Iterator<CopyNumberTriStateSegmentRecord>> allSegmentRecordIterators =
                allSegments.entrySet().stream()
                    .map(e -> e.getValue().stream()
                            .map(s -> new CopyNumberTriStateSegmentRecord(e.getKey(), s))
                            .iterator())
                    .collect(Collectors.toList());

        final Comparator<CopyNumberTriStateSegmentRecord> recordComparator = (o1, o2) -> {
            final CopyNumberTriStateSegment s1 = o1.getSegment();
            final CopyNumberTriStateSegment s2 = o2.getSegment();
            // Using the index-range.from we make sure we sort first by the contig over
            // the start position in the same order as contigs are present in
            // the input data or target collection.
            final IndexRange ir1 = targets.indexRange(o1.getSegment());
            final IndexRange ir2 = targets.indexRange(o2.getSegment());

            final int fromCmp = Integer.compare(ir1.from, ir2.from);
            if (fromCmp != 0) {
                return fromCmp;
            }
            // if fromCmp == 0, they must be in the same contig,
            // then we can sort based on segment start and then end.
            final int startCmp = Integer.compare(s1.getStart(), s2.getStart());
            if (startCmp != 0) {
                return startCmp;
            }
            return Integer.compare(s1.getEnd(), s2.getEnd());
        };

        return Iterators.mergeSorted(allSegmentRecordIterators, recordComparator);
    }

    private Map<String, List<CopyNumberTriStateSegment>> calculateBestPathSegments(final CopyNumberTriStateHiddenMarkovModel model, final ReadCountCollection inputCounts) {
        final List<String> sampleNames = inputCounts.columnNames();
        final Map<String, List<CopyNumberTriStateSegment>> allSegments = new LinkedHashMap<>(sampleNames.size());
        final List<Target> targets = inputCounts.targets();
        for (int i = 0; i < sampleNames.size(); i++) {
            final String sampleName = sampleNames.get(i);
            final List<Double> inputValues = DoubleStream.of(inputCounts.counts().getColumn(i)).boxed().collect(Collectors.toList());
            final List<CopyNumberTriState> bestPath = ViterbiAlgorithm.apply(inputValues, targets, model);
            final List<Pair<IndexRange, CopyNumberTriState>> bestPathTargetIndexRanges = condenseBestPathIntoTargetIndexAndStatePairs(bestPath, targets);
            final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult = ForwardBackwardAlgorithm.apply(inputValues, targets, model);
            final List<CopyNumberTriStateSegment> bestPathSegmentList = composeSegments(fbResult, bestPathTargetIndexRanges);
            allSegments.put(sampleName, bestPathSegmentList);
        }
        return allSegments;
    }

    /**
     * Given a plausible sequence of hidden copy-number states, it condenses it to the corresponding segments represented
     * as a pair of target index-range (which encloses the targets in the segment) and the hidden copy-number state
     * for that segment.
     *
     * @param bestPath a plausible copy number state sequence across the targets included in {@code targets}.
     * @param targets the collection of targets to take in consideration. Needed to make sure that segments don't expand
     *                across contigs.
     * @return never {@code null}.
     */
    private List<Pair<IndexRange, CopyNumberTriState>> condenseBestPathIntoTargetIndexAndStatePairs(
            final List<CopyNumberTriState> bestPath, final List<Target> targets) {

        if (bestPath.isEmpty()) {
            return Collections.emptyList();
        }

        // We approximate the expected number of segments as the square root of the path length.
        // very arbitrary but probably better than assuming there are as many segments as targets.
        final List<Pair<IndexRange, CopyNumberTriState>> result = new ArrayList<>((int) Math.ceil(Math.sqrt(bestPath.size())));

        int currentStartIndex = 0; // contains the start index of the segment being traversed.
        final ListIterator<CopyNumberTriState> pathIterator = bestPath.listIterator();
        final ListIterator<Target> targetIterator = targets.listIterator();
        String currentContig = targetIterator.next().getContig();
        CopyNumberTriState currentState = pathIterator.next();
        while (pathIterator.hasNext()) {
            final CopyNumberTriState nextState = pathIterator.next();
            final String nextContig = targetIterator.next().getContig();
            final boolean contigChanged = ! currentContig.equals(nextContig);
            final boolean stateChanged = nextState != currentState;
            if (contigChanged || stateChanged) {
                final int newCurrentStartIndex = pathIterator.previousIndex();
                result.add(new ImmutablePair<>(new IndexRange(currentStartIndex, newCurrentStartIndex), currentState));
                currentStartIndex = newCurrentStartIndex;
                currentState = nextState;
                currentContig = nextContig;
            }
        }
        result.add(new ImmutablePair<>(new IndexRange(currentStartIndex, bestPath.size()), currentState));
        return result;
    }

    /**
     * Compose the list of segments based on a inferred best hidden state sequence and
     *   the result of running the forward-backward algorithm.
     * @param fbResult the result of the Forward-backward algorithm on the same input data.
     * @param bestPathSegments the best hidden-state path along the int targets and count sequence.
     * @return never {@code null}.
     */
    private List<CopyNumberTriStateSegment> composeSegments(final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult, final List<Pair<IndexRange, CopyNumberTriState>> bestPathSegments) {
        return bestPathSegments.stream()
                .map(ir -> composeSegment(fbResult, ir.getLeft(), ir.getRight()))
                .collect(Collectors.toList());
    }

    /**
     * Composes the segment calculating all the corresponding quality scores,
     *
     * @param fbResult the {@link ForwardBackwardAlgorithm} execution result object.
     * @param targetIndexRange the target position index range of the segment.
     * @param call the call for the segment.
     * @return never {@code null}
     */
    private CopyNumberTriStateSegment composeSegment(final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult,
                                                     final IndexRange targetIndexRange, final CopyNumberTriState call) {
        final List<Double> values = fbResult.data().subList(targetIndexRange.from, targetIndexRange.to);
        final double mean = values.stream().mapToDouble(d -> d).average().orElse(Double.NaN);
        final double stdDev = GATKProtectedMathUtils.stdDev(values);
        final List<Target> targets = fbResult.positions().subList(targetIndexRange.from, targetIndexRange.to);
        final int length = targetIndexRange.size();
        final double logExactProbability = fbResult.logProbability(targetIndexRange.from, targetIndexRange.to, call);
        final double logSomeProbability = logSomeProbability(targetIndexRange.from, length, call, fbResult);
        final double logStartProbability = logStartProbability(targetIndexRange.from, call, fbResult);
        final double logEndProbability = logEndProbability(targetIndexRange.to - 1, call, fbResult);
        final double logNeutralProbability = fbResult.logProbability(targetIndexRange.from, Collections.nCopies(targets.size(), CopyNumberTriState.NEUTRAL));
        return new CopyNumberTriStateSegment(targets, mean, stdDev, call,
                logToPhredScore(logExactProbability, true),
                logToPhredScore(logSomeProbability, true),
                logToPhredScore(logStartProbability, true),
                logToPhredScore(logEndProbability, true),
                logToPhredScore(logNeutralProbability, false));
    }

    /**
     * Threshold used to determine best way to calculate log(1- exp(a))
     * based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
     */
    private static final double LN_1_M_EXP_THRESHOLD = - Math.log(2);

    /**
     * Maximum tolerated log probability value.
     * <p>Values between this
     * constant and 0.0 as considered as a 0.0.
     * </p>
     * Values above this threshold
     * are considered to indicate a log probability calculation problem that is
     * worth to be reported as an error or warning.
     */
    private static final double MAX_LOG_PROB = 1e-3;

    /**
     * Transform a log scaled probability (x) into the Phred scaled
     * equivalent or its complement (1-x) Phred scaled equivalent.
     * <p>
     *     This method tolerates probabilities slightly larger than 1.0
     *     (> 0.0 in log scale) which may occur occasionally due to
     *     float point calculation rounding.
     * </p>
     * <p>
     *     The value returned is a phred score capped by {@link #MAX_QUAL_SCORE}.
     * </p>
     *
     * @param rawLogProb the probability.
     * @param complement whether to return the direct Phred transformation ({@code false})
     *                    or its complenent ({@code true)}.
     * @return a values between 0 and {@link #MAX_QUAL_SCORE}.
     * @throws GATKException if {@code rawLogProb} is larger than {@link #MAX_LOG_PROB}.
     */
    private double logToPhredScore(final double rawLogProb, final boolean complement) {
        if (rawLogProb > MAX_LOG_PROB) {
            throw new GATKException(
                    String.format("numerical instability problem: the log-probability is too large: %g > 0.0 (with maximum tolerance %g)", rawLogProb, MAX_LOG_PROB));
        }
        // make sure that log probs are less than 1 in linear scale.
        // there are cases in that they are just over 0.0 due to float point precision.
        final double logProbEqOrLessThan0 = Math.min(0.0, rawLogProb);

        // Accurate way to calculate log(1-exp(a))
        // based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
        final double finalLogProb = complement
                ? logComplement(logProbEqOrLessThan0)
                : logProbEqOrLessThan0;

        final double absoluteQualScore = QualityUtils.phredScaleLog10ErrorRate(finalLogProb * INV_LN_10);
        return Math.min(MAX_QUAL_SCORE, absoluteQualScore);
    }

    /**
     * Calculates the complement of a log probability.
     *
     * <p>
     *     With complement of {@code x} we mean: {@code log(1-log(x))}.
     * </p>
     * @param x the input log probability.
     * @return {@code log(1-log(x))}
     */
    private double logComplement(final double x) {
        return x >= LN_1_M_EXP_THRESHOLD
            ? Math.log(-Math.expm1(x))
            : Math.log1p(-Math.exp(x));
    }

    /**
     * Calculates the probability of a hidden state switch at the beginning of the segment
     * from a different state to the call state.
     */
    private double logStartProbability(final int firstTarget,
                                       final CopyNumberTriState call,
                                       final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        if (firstTarget == 0) {
            return fbResult.logProbability(firstTarget, call);
        } else {
            final EnumSet<CopyNumberTriState> onlyCallState = EnumSet.of(call);
            final EnumSet<CopyNumberTriState> allStatesExceptCall = EnumSet.complementOf(onlyCallState);
            return fbResult.logConstrainedProbability(firstTarget - 1, Arrays.asList(allStatesExceptCall, onlyCallState));
        }
    }

    /**
     * Calculates the probability of a hidden state switch at the end of the segment
     * from the call state to a different state.
     */
    private double logEndProbability(final int lastTarget,
                                     final CopyNumberTriState call,
                                     final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        if (lastTarget + 1 == fbResult.positions().size()) {
            return fbResult.logProbability(lastTarget, call);
        } else {
            final EnumSet<CopyNumberTriState> onlyCallState = EnumSet.of(call);
            final EnumSet<CopyNumberTriState> allStatesExceptCall = EnumSet.complementOf(onlyCallState);
            return fbResult.logConstrainedProbability(lastTarget, Arrays.asList(onlyCallState, allStatesExceptCall));
        }
    }

    private double logSomeProbability(final int firstTarget, final int length, final CopyNumberTriState call,
                                          final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        // trivial case when the segment has length 0.
        if (length == 0) {
            return 0;
        } else {
            final Set<CopyNumberTriState> otherStates = EnumSet.complementOf(EnumSet.of(call));
            final List<Set<CopyNumberTriState>> otherStatesConstraints = Collections.nCopies(length, otherStates);
            final double logOtherStates = fbResult.logConstrainedProbability(firstTarget, otherStatesConstraints);
            return logComplement(logOtherStates);
        }
    }

    private void checkForMissingTargetsInInputCounts(final TargetCollection<Target> targets, final ReadCountCollection inputCounts) {
        if (inputCounts.targets().size() != targets.targetCount()) {
            logger.warn(String.format("The input counts are missing some of the targets required (%d): e.g. %s",
                    inputCounts.targets().size() - targets.targetCount(), targets.targets().stream()
                            .filter(t -> !inputCounts.targets().contains(t))
                            .limit(5)));
        }
    }

    private ReadCountCollection readAndSortCountsByTargetCoordinates(final TargetCollection<Target> targets) {

        // read the input read-counts.
        final ReadCountCollection originalReadCounts = readInputCounts(targets);

        // Get the input read count target positions in the input targets collection.
        final int[] sortedPositions = originalReadCounts.targets().stream()
                .mapToInt(t -> targets.index(t.getName())).toArray();

        // Check whether the input is in the same order as the targets in the input target collection.
        // If not we rearrange them accordingly.
        for (int i = 1, previousIndex = -1; i < sortedPositions.length; i++) {
            if (sortedPositions[i] == -1) {
                return originalReadCounts.arrangeTargets(targets.targets());
            } else {
                if (previousIndex >= sortedPositions[i]) {
                    return originalReadCounts.arrangeTargets(targets.targets());
                }
                previousIndex = sortedPositions[i];
            }
        }
        return originalReadCounts;
    }

    private ReadCountCollection readInputCounts(final TargetCollection<Target> targets) {
        try {
             return ReadCountCollectionUtils.parse(inputFile, targets, true);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, ex);
        }
    }
}
