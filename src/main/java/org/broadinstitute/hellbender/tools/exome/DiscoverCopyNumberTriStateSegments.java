package org.broadinstitute.hellbender.tools.exome;

import com.google.common.collect.Iterators;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModel;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModelArgumentCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
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

    public static final String INPUT_ARE_ZSCORES_FULL_NAME = "inputAreZscores";

    public static final String INPUT_ARE_ZSCORES_SHORT_NAME = "zscores";

    private static final double INV_LN_10 = 1.0 / Math.log(10);

    /**
     * Maximum reportable output quality score; higher quality scores will
     * capped to this value.
     */
    public static final double MAX_QUAL_SCORE = 199.99999999999;

    @ArgumentCollection
    protected CopyNumberTriStateHiddenMarkovModelArgumentCollection modelArguments =
            new CopyNumberTriStateHiddenMarkovModelArgumentCollection();

    @Argument(
            doc = "Indicate whether the input is expressed in z-scores",
            fullName = INPUT_ARE_ZSCORES_FULL_NAME,
            shortName = INPUT_ARE_ZSCORES_SHORT_NAME,
            optional = false
    )
    protected boolean inputAreZscores;

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
        checkForMissingTargetsInInputCounts(targets, inputCounts);

        final Map<String, List<CopyNumberTriStateSegment>> allSegmentsBySampleName =
                calculateBestPathSegments(model, inputCounts);

        final Iterator<CopyNumberTriStateSegmentRecord> allSegmentRecordsSortedByCoordinates =
                composeATargetSortedSegmentRecordIterator(targets, allSegmentsBySampleName);

        writeSegmentsInOutputFile(allSegmentRecordsSortedByCoordinates);

        return "SUCCESS";
    }

    /**
     * Writes sorted segments into the output file provided by the user.
     * @param allSegmentsSorted iterator that must produce segment records in order.
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CouldNotCreateOutputFile if there is any
     * issue creating or writing into the output file provided by the user.
     */
    private void writeSegmentsInOutputFile(final Iterator<CopyNumberTriStateSegmentRecord> allSegmentsSorted) {
        try (final CopyNumberTriStateSegmentRecordWriter writer = new CopyNumberTriStateSegmentRecordWriter(outputFile)) {
            while (allSegmentsSorted.hasNext()) {
                writer.writeRecord(allSegmentsSorted.next());
            }
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

    private Map<String, List<CopyNumberTriStateSegment>> calculateBestPathSegments(CopyNumberTriStateHiddenMarkovModel model, ReadCountCollection inputCounts) {
        final List<String> sampleNames = inputCounts.columnNames();
        final Map<String, List<CopyNumberTriStateSegment>> allSegments = new LinkedHashMap<>(sampleNames.size());
        long recordCount = 0;
        final List<Target> targets = inputCounts.targets();
        for (int i = 0; i < sampleNames.size(); i++) {
            final String sampleName = sampleNames.get(i);
            final double[] inputValues = inputCounts.counts().getColumn(i);
            final List<CopyNumberTriStateSegment> bestPathSegmentList = getCopyNumberTriStateSegments(model, targets, inputValues);
            recordCount += bestPathSegmentList.size();
            allSegments.put(sampleName, bestPathSegmentList);
        }
        logger.info(String.format("Found a total of %d segments across %d samples", recordCount, inputCounts.columnNames().size()));
        return allSegments;
    }

    /**
     * Calculates the list of copy-number change segments as the one present in the
     * most likely segment sequence given a {@link CopyNumberTriStateHiddenMarkovModel model}.
     * @param model the HMM model to be used to discover copy-number change segments.
     * @param targets the input targets sequence.
     * @param inputValues the input coverage values.
     * @return never {@code null}.
     */
    private List<CopyNumberTriStateSegment> getCopyNumberTriStateSegments(final CopyNumberTriStateHiddenMarkovModel model, final List<Target> targets, double[] inputValues) {
        final double mean = MathUtils.mean(inputValues, 0, inputValues.length);
        final double stdev =  MathUtils.stddev(inputValues, 0, inputValues.length);
        final List<Double> inputValuesList = DoubleStream.of(inputValues).boxed().collect(Collectors.toList());
        final List<Double> standardizedValuesList = inputAreZscores ? inputValuesList :
                inputValuesList.stream().map(v -> (v - mean) / stdev).collect(Collectors.toList());
        final List<CopyNumberTriState> bestPath = ViterbiAlgorithm.apply(standardizedValuesList, targets, model);
        final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult = ForwardBackwardAlgorithm.apply(standardizedValuesList, targets, model);
        return composeSegments(bestPath, fbResult, inputValuesList);
    }

    /**
     * Compose the list of segments based on a inferred best hidden state sequence and
     *   the result of running the forward-backward algorithm.
     * @param bestPath the best hidden-state path along the int targets and count sequence.
     * @param fbResult the result of the Forward-backward algorithm on the same input data.
     * @return never {@code null}.
     */
    private List<CopyNumberTriStateSegment> composeSegments(final List<CopyNumberTriState> bestPath,
                                                            final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult,
                                                            final List<Double> inputValues) {
        final List<Target> targets = fbResult.positions();
        // trivial case when the number of targets analyzed is 0.
        if (targets.size() == 0) {
            return new ArrayList<>(0);
        }

        final List<Double> coverage = fbResult.data();
        final List<CopyNumberTriStateSegment> result = new ArrayList<>(bestPath.size());
        final List<Target> currentSegmentTargets = new ArrayList<>(bestPath.size());
        final List<Double> currentSegmentData = new ArrayList<>(bestPath.size());
        final List<Double> currentSegmentInputValues = new ArrayList<>(bestPath.size());
        CopyNumberTriState currentSegmentState = bestPath.get(0);
        Target lastTargetInPreviousSegment = null;
        currentSegmentTargets.add(targets.get(0));
        currentSegmentData.add(coverage.get(0));
        currentSegmentInputValues.add(inputValues.get(0));
        for (int i = 1; i < targets.size(); i++) {
            final CopyNumberTriState thisState = bestPath.get(i);
            final Double thisDatum = coverage.get(i);
            final Target target = targets.get(i);
            if (currentSegmentState == thisState
                    && target.getContig().equals(currentSegmentTargets.get(0).getContig())) {
                currentSegmentTargets.add(target);
                currentSegmentData.add(thisDatum);
                currentSegmentInputValues.add(inputValues.get(i));
            } else {
                result.add(composeSegment(currentSegmentTargets, currentSegmentState,
                                          lastTargetInPreviousSegment, target, fbResult, currentSegmentInputValues));
                lastTargetInPreviousSegment = currentSegmentTargets.get(currentSegmentTargets.size() - 1);
                currentSegmentTargets.clear();
                currentSegmentTargets.add(target);
                currentSegmentData.clear();
                currentSegmentData.add(thisDatum);
                currentSegmentInputValues.clear();
                currentSegmentInputValues.add(inputValues.get(i));
                currentSegmentState = thisState;
            }
        }
        result.add(composeSegment(currentSegmentTargets, currentSegmentState,
                                  lastTargetInPreviousSegment, null, fbResult, currentSegmentInputValues));
        return result;
    }

    /**
     * Composes the segment calculating all the corresponding quality scores,
     *  coverage mean and stddev.
     *
     * @param targets the targets enclosed in the segment.
     * @param call the segment called state.
     * @param precedingTarget the target preceding the segment's start target.
     * @param followingTarget the target following the segment's last target.
     * @param fbResult the result of running the Forward-Backward algorithm.
     * @return never {@code null}.
     */
    private CopyNumberTriStateSegment composeSegment(final List<Target> targets,
                                                     final CopyNumberTriState call,
                                                     final Target precedingTarget, final Target followingTarget,
                                                     final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult,
                                                     final List<Double> inputValues) {
        final int length = targets.size();
        final Target firstTarget = targets.get(0);
        final double mean = inputValues.stream().mapToDouble(d -> d).average().getAsDouble();
        final double stdev = inputValues.size() <= 1 ? 0 : Math.sqrt(inputValues.stream().mapToDouble(d -> (d - mean) * (d - mean)).sum() / (length - 1));
        final double logExactProbability = logExactProbability(firstTarget, length, call, fbResult);
        final double logSomeProbability = logSomeProbability(firstTarget, length, call, fbResult);
        final double logStartProbability = logStartProbability(precedingTarget, firstTarget, call, fbResult);
        final double logEndProbability = logEndProbability(targets.get(length - 1), followingTarget, call, fbResult);
        final double logNeutralProbability = fbResult.logProbability(firstTarget, Collections.nCopies(targets.size(), CopyNumberTriState.NEUTRAL));
        return new CopyNumberTriStateSegment(targets, mean, stdev, call,
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
            throw new GATKException(String.format("numerical instability problem: the log-probability is too large: %g > 0.0 (with maximum tolerance %g)", rawLogProb, MAX_LOG_PROB));
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
    private double logStartProbability(final Target precedingTarget, final Target firstTarget,
                                       final CopyNumberTriState call,
                                       final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        if (precedingTarget == null || !precedingTarget.getContig().equals(firstTarget.getContig())) {
            // at the very beginning? then we just get the posterior of the first target to
            // have the call state.
            return fbResult.logProbability(firstTarget, call);
        } else {
            final EnumSet<CopyNumberTriState> onlyCallState = EnumSet.of(call);
            final EnumSet<CopyNumberTriState> allStatesExceptCall = EnumSet.complementOf(onlyCallState);
            return fbResult.logConstrainedProbability(precedingTarget, Arrays.asList(allStatesExceptCall, onlyCallState));
        }
    }

    /**
     * Calculates the probability of a hidden state switch at the end of the segment
     * from the call state to a different state.
     */
    private double logEndProbability(final Target lastTarget, final Target followingTarget,
                                     final CopyNumberTriState call,
                                     final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        if (followingTarget == null || !followingTarget.getContig().equals(lastTarget.getContig())) {
            return fbResult.logProbability(lastTarget, call);
        } else {
            final EnumSet<CopyNumberTriState> onlyCallState = EnumSet.of(call);
            final EnumSet<CopyNumberTriState> allStatesExceptCall = EnumSet.complementOf(onlyCallState);
            return fbResult.logConstrainedProbability(lastTarget, Arrays.asList(onlyCallState, allStatesExceptCall));
        }
    }

    /**
     * This method calculates the probability that in a given interval targets/data pairs, we
     * pass at least once thru the "call state" the way XHMM does.
     *
     * <p>
     *     In XHMM the "Some" quality only considers the copy-neutral and the call state in the calculations.
     *     So if we are looking into a deletion event, paths that go thru the duplication state are ignored.
     * </p>
     * <p>
     *     XHMM does not give definition for this score when the call state is the neutral state.
     *     Here we just return 0.
     * </p>
     *
     * @param firstTarget the first target in the segment of interest.
     * @param call the call for the segment.
     * @param fbResult the result from running the {@link ForwardBackwardAlgorithm#apply forward-backward algorithm}
     *                 on the original data.
     * @return a valid log scaled probability from 0 to -Inf.
     */
    private double logXHMMSomeProbability(final Target firstTarget, final int length, final CopyNumberTriState call,
                                          final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        // trivial case when the segment has length 0.
        if (length == 0 || call == CopyNumberTriState.NEUTRAL) {
            return 0;
        } else {
            final Set<CopyNumberTriState> callOrNeutralStates =  EnumSet.of(call, CopyNumberTriState.NEUTRAL);
            final List<Set<CopyNumberTriState>> callOrNeutralConstraints = Collections.nCopies(length, callOrNeutralStates);

            final double logCallOrNeutralProbability = fbResult.logConstrainedProbability(firstTarget, callOrNeutralConstraints);
            final double logAllNeutralProbability = fbResult.logProbability(firstTarget, Collections.nCopies(length, CopyNumberTriState.NEUTRAL));

            // log-sub-exp trick: log(exp(a) - exp(b)) = a + log(1 - exp(b - a))
            return logCallOrNeutralProbability + Math.log1p(-Math.exp(logAllNeutralProbability - logCallOrNeutralProbability));
        }
    }

    private double logSomeProbability(final Target firstTarget, final int length, final CopyNumberTriState call,
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

    /**
     * Probability of the segment to be in the call state all the way through.
     *
     * @param firstTarget the first target in the segment.
     * @param length the length of the segment in targets.
     * @param call the call of the segment.
     * @param fbResult the forward-backward algorithm result object.
     * @return a log scaled probability between 0 to -Inf.
     */
    private double logExactProbability(final Target firstTarget, final int length, final CopyNumberTriState call, final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        return fbResult.logProbability(firstTarget, Collections.nCopies(length, call));
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
