package org.broadinstitute.hellbender.utils.hmm.segmentation;

import com.google.common.collect.Iterators;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;

import javax.annotation.Nonnull;
import java.io.IOException;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO github/gatk-protected issue #855 -- this is a last-minute fork of {@link HMMPostProcessor}.
 * It just does segmentation (no VCF creation)
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class HMMSegmentProcessor<DATA, STATE extends CallStringProducer & ScalarProducer, TARGET extends Target> {

    private final Logger logger = LogManager.getLogger(HMMSegmentProcessor.class);

    /**
     * Maximum tolerated log probability value. Values between this constant and 0.0 as considered as a 0.0.
     * Values above this threshold are considered to indicate a log probability calculation problem that is
     * worth to be reported as an error or warning.
     */
    public static final double MAX_LOG_PROB = 1e-3;

    /**
     * Maximum reportable output quality score; higher quality scores will
     * capped to this value.
     */
    public static final double MAX_QUAL_SCORE = 9999.99999999999;

    /**
     * Precision of Phred-scale scores
     */
    public static final double PHRED_SCORE_PRECISION = 0.1;

    /**
     * Calculating variance posterior on long segments is time consuming; we assume hidden states become
     * independent if their separation is larger than the following value
     */
    public static final double INDEPENDENT_TARGETS_SEPARATION_THRESHOLD = 0;

    private final List<String> sampleNames;
    private final List<SexGenotypeData> sampleSexGenotypes;
    private final BiFunction<SexGenotypeData, TARGET, STATE> referenceStateFactory;
    private final int numSamples;
    private final List<TargetCollection<TARGET>> sampleTargets;
    private final List<ForwardBackwardAlgorithm.Result<DATA, TARGET, STATE>> sampleForwardBackwardResults;
    private final List<List<STATE>> sampleBestPaths;
    private final Map<String, List<HiddenStateSegment<STATE, TARGET>>> allSegmentsBySampleName;

    /**
     * This "master" target collection is the lexicographically sorted union of all targets from all samples
     */
    private final TargetCollection<TARGET> masterTargetCollection;

    /**
     * Public constructor.
     *
     * @param sampleNames list of sample names
     * @param sampleSexGenotypes list of sample sex genotype data
     * @param sampleTargets list of target collection for each sample
     * @param sampleForwardBackwardResults list of forward-backward result for each sample
     * @param sampleBestPaths list of best hidden state sequence call for each sample
     * @param referenceStateFactory a bi-function from (sex genotype data, target) to the "reference" state.
     */
    public HMMSegmentProcessor(
            @Nonnull final List<String> sampleNames,
            @Nonnull final List<SexGenotypeData> sampleSexGenotypes,
            @Nonnull final BiFunction<SexGenotypeData, TARGET, STATE> referenceStateFactory,
            @Nonnull final List<TargetCollection<TARGET>> sampleTargets,
            @Nonnull final List<ForwardBackwardAlgorithm.Result<DATA, TARGET, STATE>> sampleForwardBackwardResults,
            @Nonnull final List<List<STATE>> sampleBestPaths) {

        /* basic consistency checks on the input data */
        numSamples = sampleNames.size();
        Utils.nonNull(sampleNames, "The list of sample names must be non-null");
        Utils.nonNull(sampleSexGenotypes, "The list of sample sex genotypes must be non-null");
        Utils.nonNull(referenceStateFactory, "The reference state factory must be non-null");
        Utils.nonNull(sampleTargets, "The list of sample target list must be non-null");
        Utils.nonNull(sampleForwardBackwardResults, "The list of forward-backward results must be non-null");
        Utils.nonNull(sampleBestPaths, "The list of sample best paths must be non-null");
        Utils.validateArg(numSamples > 0, "Number of samples must be > 0");
        Utils.validateArg(new HashSet<>(sampleNames).size() == numSamples, "Sample names must be unique");
        Utils.validateArg(sampleTargets.size() == numSamples, "List of targets per sample must have the" +
                " same length as the list of sample names");
        Utils.validateArg(sampleSexGenotypes.size() == numSamples, "List of sample sex genotypes must have the" +
                " same length as the list of sample names");
        Utils.validateArg(sampleForwardBackwardResults.size() == numSamples, "List of forward-backward" +
                " result per sample must have the same length as the list of sample names");
        Utils.validateArg(sampleBestPaths.size() == numSamples, "List of best paths per sample must have" +
                " the same length as the list of sample names");
        Utils.validateArg(IntStream.range(0, numSamples)
                        .allMatch(si -> sampleTargets.get(si).targetCount() == sampleBestPaths.get(si).size()),
                "Some of the provided sample best paths have different lengths than their corresponding target list");
        Utils.validateArg(IntStream.range(0, numSamples)
                .allMatch(si -> sampleNames.get(si).equals(sampleSexGenotypes.get(si).getSampleName())),
                "Some of the sample names in the sex genotype data do not correspond to the list of sample names");

        /* get the hidden states from the first sample and ensure all samples have the same hidden states */
        final Set<STATE> pooledStates = sampleForwardBackwardResults.stream()
                .map(fb -> fb.model().hiddenStates())
                .flatMap(List::stream)
                .collect(Collectors.toSet());
        Utils.validateArg(sampleForwardBackwardResults.stream().allMatch(fb ->
                        new HashSet<>(fb.model().hiddenStates()).equals(pooledStates)), "All samples must have the" +
                " same set of hidden states");

        /* TODO in the future, we may want to replace lexicographical order with natural order */
        /* assert that the target collection of each sample is lexicographically sorted */
        IntStream.range(0, numSamples)
                .forEach(si -> {
                    final List<TARGET> sampleTargetList = sampleTargets.get(si).targets();
                    Utils.validateArg(IntStream.range(0, sampleTargetList.size() - 1)
                                    .allMatch(ti -> IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR
                                            .compare(sampleTargetList.get(ti + 1), sampleTargetList.get(ti)) > 0),
                            "The target collection for sample number " + si + " is not lexicographically sorted");
                });

        /* store sample data verbatim */
        this.sampleNames = sampleNames;
        this.sampleTargets = sampleTargets;
        this.sampleForwardBackwardResults = sampleForwardBackwardResults;
        this.sampleBestPaths = sampleBestPaths;
        this.referenceStateFactory = referenceStateFactory;
        this.sampleSexGenotypes = sampleSexGenotypes;

        /* pool targets from all samples and sort them lexicographically to make a master target collection */
        masterTargetCollection = new HashedListTargetCollection<>(sampleTargets.stream()
                .map(TargetCollection::targets)
                .flatMap(List::stream)
                .distinct()
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList()));

        allSegmentsBySampleName = calculateBestPathSegments();
    }

    /**
     * Write segments to a table writer
     *
     * @param writer an instance of {@link HiddenStateSegmentRecordWriter}
     */
    public void writeSegmentsToTableWriter(@Nonnull final HiddenStateSegmentRecordWriter<STATE, TARGET> writer) {
        final Iterator<HiddenStateSegmentRecord<STATE, TARGET>> allSegmentRecordsSortedByCoordinates =
                composeTargetSortedSegmentRecordIterator(masterTargetCollection, allSegmentsBySampleName);
        logger.info("Writing segments to file...");
        writeSegments(allSegmentRecordsSortedByCoordinates, writer);
    }

    /**
     * Return segments as a list of {@link HiddenStateSegmentRecord}
     *
     * @return a list of
     */
    public List<HiddenStateSegmentRecord<STATE, TARGET>> getSegmentsAsList() {
        final List<HiddenStateSegmentRecord<STATE, TARGET>> segmentList = new ArrayList<>();
        composeTargetSortedSegmentRecordIterator(masterTargetCollection, allSegmentsBySampleName)
                .forEachRemaining(segmentList::add);
        return segmentList;
    }

    /* segmentation-related methods */

    /**
     * Takes the provided forward-backward result and the best hidden state chain for each sample and
     * generates segments
     *
     * @return a map from sample names to a list of segments
     */
    private Map<String, List<HiddenStateSegment<STATE, TARGET>>> calculateBestPathSegments() {
        final Map<String, List<HiddenStateSegment<STATE, TARGET>>> allSegments = new LinkedHashMap<>(sampleNames.size());
        IntStream.range(0, numSamples).parallel()
                .forEach(sampleIndex -> {
                    logger.info("Current sample: " + sampleNames.get(sampleIndex) + "...");
                    final String sampleName = sampleNames.get(sampleIndex);
                    final TargetCollection<TARGET> targets = sampleTargets.get(sampleIndex);
                    final List<STATE> bestPath = sampleBestPaths.get(sampleIndex);
                    final ForwardBackwardAlgorithm.Result<DATA, TARGET, STATE> fbResult = sampleForwardBackwardResults.get(sampleIndex);
                    final List<Pair<IndexRange, STATE>> bestPathTargetIndexRanges =
                            condenseBestPathIntoTargetIndexAndStatePairs(bestPath, targets);
                    final List<HiddenStateSegment<STATE, TARGET>> bestPathSegmentList =
                            composeSegments(sampleIndex, fbResult, bestPathTargetIndexRanges);
                    allSegments.put(sampleName, bestPathSegmentList);
                });
        return allSegments;
    }

    /**
     * Given a plausible sequence of hidden states, it condenses the sequence to the corresponding segments
     * represented as a pair of target index-range (which encloses the targets in the segment) and the hidden
     * states for that segment.
     *
     * @param bestPath a plausible copy number state sequence across the targets included in {@code targets}.
     * @param targets the collection of targets to take in consideration. Needed to make sure that segments don't
     *                expand across contigs.
     * @return never {@code null}.
     */
    private List<Pair<IndexRange, STATE>> condenseBestPathIntoTargetIndexAndStatePairs(final List<STATE> bestPath,
                                                                                       final TargetCollection<TARGET> targets) {
        if (bestPath.isEmpty()) {
            return Collections.emptyList();
        }

        /* estimate the number of segments */
        final long numChangePoints = IntStream.range(0, bestPath.size() - 1)
                .filter(i -> !bestPath.get(i).equals(bestPath.get(i + 1)))
                .count();
        final List<Pair<IndexRange, STATE>> result = new ArrayList<>((int) numChangePoints);

        int currentStartIndex = 0; // contains the start index of the segment being traversed.
        final ListIterator<STATE> pathIterator = bestPath.listIterator();
        final ListIterator<TARGET> targetIterator = targets.targets().listIterator();
        String currentContig = targetIterator.next().getContig();
        STATE currentState = pathIterator.next();
        while (pathIterator.hasNext()) {
            final STATE nextState = pathIterator.next();
            final String nextContig = targetIterator.next().getContig();
            final boolean contigChanged = !currentContig.equals(nextContig);
            final boolean stateChanged = !nextState.equals(currentState);
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
     * the result of running the forward-backward algorithm.
     *
     * @param fbResult the result of the forward-backward algorithm on the same input data.
     * @param bestPathSegments the best hidden-state path along the int targets and count sequence.
     * @return never {@code null}.
     */
    private List<HiddenStateSegment<STATE, TARGET>> composeSegments(
            final int sampleIndex,
            final ForwardBackwardAlgorithm.Result<DATA, TARGET, STATE> fbResult,
            final List<Pair<IndexRange, STATE>> bestPathSegments) {
        return bestPathSegments.stream()
                .map(ir -> composeHiddenStateSegmentFromTargetIndexRange(sampleIndex, fbResult, ir.getLeft(),
                        ir.getRight()))
                .collect(Collectors.toList());
    }

    /**
     * Composes the segment calculating all the corresponding quality scores,
     *
     * @param sampleIndex sample index (used for inferring the reference state)
     * @param fbResult the {@link ForwardBackwardAlgorithm} execution result object.
     * @param targetIndexRange the target position index range of the segment.
     * @param call the call for the segment.
     * @return never {@code null}
     */
    private HiddenStateSegment<STATE, TARGET> composeHiddenStateSegmentFromTargetIndexRange(
            final int sampleIndex,
            final ForwardBackwardAlgorithm.Result<DATA, TARGET, STATE> fbResult,
            final IndexRange targetIndexRange,
            final STATE call) {

        final int segmentLength = targetIndexRange.size();
        final double mean = calculateSegmentMean(targetIndexRange.from, segmentLength, fbResult);
        final double stdDev = FastMath.sqrt(calculateSegmentVariance(targetIndexRange.from, segmentLength, fbResult));

        final List<TARGET> targets = fbResult.positions().subList(targetIndexRange.from, targetIndexRange.to);
        final double logExactProbability = fbResult.logProbability(targetIndexRange.from, targetIndexRange.to, call);
        final double logSomeProbability = logSomeProbability(targetIndexRange.from, segmentLength, call, fbResult);
        final double logStartProbability = logStartProbability(targetIndexRange.from, call, fbResult);
        final double logEndProbability = logEndProbability(targetIndexRange.to - 1, call, fbResult);
        final Set<String> contigsSet = targets.stream().map(Target::getContig).collect(Collectors.toSet());
        if (contigsSet.size() > 1) {
            throw new IllegalArgumentException("The target index range for the segment encompasses multiple contigs");
        }
        final STATE referenceState = referenceStateFactory.apply(sampleSexGenotypes.get(sampleIndex), targets.get(0));
        final double logReferenceProbability = fbResult.logProbability(targetIndexRange.from,
                Collections.nCopies(segmentLength, referenceState));

        return new HiddenStateSegment<>(targets, mean, stdDev, call,
                logProbToPhredScore(logExactProbability, true),
                logProbToPhredScore(logSomeProbability, true),
                logProbToPhredScore(logStartProbability, true),
                logProbToPhredScore(logEndProbability, true),
                logProbToPhredScore(logReferenceProbability, false));
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
     * @param masterTargetCollection the target collection that encompasses targets from all samples
     * @param allSegments map with all segment lists coming from different samples.
     * @return never null.
     */
    private Iterator<HiddenStateSegmentRecord<STATE, TARGET>> composeTargetSortedSegmentRecordIterator(
            final TargetCollection<TARGET> masterTargetCollection,
            final Map<String, List<HiddenStateSegment<STATE, TARGET>>> allSegments) {
        final List<Iterator<HiddenStateSegmentRecord<STATE, TARGET>>> allSegmentRecordIterators =
                allSegments.entrySet().stream()
                        .map(e -> e.getValue().stream()
                                .map(s -> new HiddenStateSegmentRecord<>(e.getKey(), s))
                                .iterator())
                        .collect(Collectors.toList());

        final Comparator<HiddenStateSegmentRecord<STATE, TARGET>> recordComparator = (o1, o2) -> {
            final HiddenStateSegment<STATE, TARGET> s1 = o1.getSegment();
            final HiddenStateSegment<STATE, TARGET> s2 = o2.getSegment();
            // Using the index-range.from we make sure we sort first by the contig over
            // the start position in the same order as contigs are present in
            // the input data or target collection.
            final IndexRange ir1 = masterTargetCollection.indexRange(o1.getSegment());
            final IndexRange ir2 = masterTargetCollection.indexRange(o2.getSegment());

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

    /**
     * Writes sorted segments into the output file provided by the user.
     *
     * @param allSegmentsSorted iterator that must produce segment records in order.
     * @param writer a table writer
     * @throws org.broadinstitute.hellbender.exceptions.UserException.CouldNotCreateOutputFile if there is any
     * issue creating or writing into the output file provided by the user.
     */
    private void writeSegments(final Iterator<HiddenStateSegmentRecord<STATE, TARGET>> allSegmentsSorted,
                               final HiddenStateSegmentRecordWriter<STATE, TARGET> writer) {
        try {
            long count = 0;
            final Set<String> samples = new HashSet<>();
            while (allSegmentsSorted.hasNext()) {
                final HiddenStateSegmentRecord<STATE, TARGET> record = allSegmentsSorted.next();
                samples.add(record.getSampleName());
                writer.writeRecord(record);
                count++;
            }
            logger.info(String.format("Found a total of %d segments across %d samples", count, samples.size()));

        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("problems writing segments", ex);
        }
    }

    /* methods related to calculating various quality scores, probabilities, and statistics */

    /**
     * Calculates the probability that some of the targets in the range [firstTargetIndex, firstTargetIndex + segmentLength)
     * are called as a hidden state {@param call}
     *
     * @param firstTargetIndex segment start target index
     * @param segmentLength number of targets in the segment
     * @param call the called state
     * @param fbResult forward-backward result
     * @param <STATE> type of hidden state
     * @param <TARGET> type of target
     * @return a double
     */
    public static <STATE, TARGET extends Locatable> double logSomeProbability(
            final int firstTargetIndex,
            final int segmentLength,
            final STATE call,
            final ForwardBackwardAlgorithm.Result<?, TARGET, STATE> fbResult) {
        /* trivial case when the segment has length 0 */
        if (segmentLength == 0) {
            return 0;
        } else {
            final Set<STATE> otherStates = new HashSet<>(fbResult.model().hiddenStates());
            otherStates.remove(call);
            final List<Set<STATE>> otherStatesConstraints = Collections.nCopies(segmentLength, otherStates);
            final double logOtherStates = fbResult.logConstrainedProbability(firstTargetIndex, otherStatesConstraints);
            return GATKProtectedMathUtils.logProbComplement(logOtherStates);
        }
    }

    /**
     * Calculates the probability of a hidden state switch at the beginning of the segment
     * from a different state to the call state.
     *
     * @param firstTargetIndex segment start target index
     * @param call the called state
     * @param fbResult forward-backward result
     * @param <STATE> type of hidden state
     * @param <TARGET> type of target
     * @return a double
     */
    public static <STATE, TARGET extends Locatable> double logStartProbability(
            final int firstTargetIndex,
            final STATE call,
            final ForwardBackwardAlgorithm.Result<?, TARGET, STATE> fbResult) {
        if (firstTargetIndex == 0) {
            return fbResult.logProbability(firstTargetIndex, call);
        } else {
            final Set<STATE> onlyCallState = Collections.singleton(call);
            final Set<STATE> allStatesExceptCall = new HashSet<>(fbResult.model().hiddenStates());
            allStatesExceptCall.remove(call);
            return fbResult.logConstrainedProbability(firstTargetIndex - 1, Arrays.asList(allStatesExceptCall, onlyCallState));
        }
    }

    /**
     * Calculates the probability of a hidden state switch at the end of the segment
     * from the call state to a different state.
     *
     * @param lastTargetIndex segment end target index
     * @param call the called state
     * @param fbResult forward-backward result
     * @param <STATE> type of hidden state
     * @param <TARGET> type of target
     * @return a double
     */
    public static <STATE, TARGET extends Locatable> double logEndProbability(
            final int lastTargetIndex,
            final STATE call,
            final ForwardBackwardAlgorithm.Result<?, TARGET, STATE> fbResult) {
        if (lastTargetIndex + 1 == fbResult.positions().size()) {
            return fbResult.logProbability(lastTargetIndex, call);
        } else {
            final Set<STATE> onlyCallState = Collections.singleton(call);
            final Set<STATE> allStatesExceptCall = new HashSet<>(fbResult.model().hiddenStates());
            return fbResult.logConstrainedProbability(lastTargetIndex, Arrays.asList(onlyCallState, allStatesExceptCall));
        }
    }

    /**
     * Calculates the posterior expectation of the mean of hidden state values on a segment
     *
     * @param firstTargetIndex first target index of the segment
     * @param segmentLength length of the segment
     * @param fbResult result of forward-backward algorithm
     * @param <STATE> type of hidden state; must implement {@link ScalarProducer}
     * @param <TARGET> type of target; must implement {@link Locatable}
     * @throws IllegalStateException if segment length is negative
     * @return segment mean
     */
    public static <STATE extends ScalarProducer, TARGET extends Locatable> double calculateSegmentMean(
            final int firstTargetIndex,
            final int segmentLength,
            final ForwardBackwardAlgorithm.Result<?, TARGET, STATE> fbResult) {
        Utils.validateArg(segmentLength >= 0, "Segment length must be non-negative");
        /* trivial case when the segment has length 0 */
        if (segmentLength == 0) {
            return Double.NaN;
        } else {
            final List<STATE> hiddenStates = fbResult.model().hiddenStates();
            final double[] hiddenStateValues = hiddenStates.stream()
                    .mapToDouble(STATE::getScalar)
                    .toArray();
            return IntStream.range(firstTargetIndex, firstTargetIndex + segmentLength).boxed()
                    .flatMapToDouble(targetIndex -> IntStream.range(0, hiddenStates.size())
                            .mapToDouble(si -> hiddenStateValues[si] *
                                    FastMath.exp(fbResult.logProbability(targetIndex, hiddenStates.get(si)))))
                    .sum() / segmentLength;
        }
    }

    /**
     * Calculates the posterior expectation of the variance of hidden state values on a segment:
     *
     *      E[\sigma^2] = (1/TARGET) \sum_{t} E[S_t^2] - (1/TARGET^2) \sum_{t1,t2} E[S_t1 S_t2]
     *
     * Here, TARGET = {@code segmentLength}
     *
     * @param firstTargetIndex first target index of the segment
     * @param segmentLength length of the segment
     * @param fbResult result of forward-backward algorithm
     * @param <STATE> type of hidden state; must implement {@link ScalarProducer}
     * @param <TARGET> type of target; must extent {@link Target}
     * @throws IllegalStateException if segment length is negative
     * @return segment variance
     */
    public static <STATE extends ScalarProducer, TARGET extends Target> double calculateSegmentVariance(
            final int firstTargetIndex,
            final int segmentLength,
            final ForwardBackwardAlgorithm.Result<?, TARGET, STATE> fbResult) {
        Utils.validateArg(segmentLength >= 0, "Segment length must be non-negative");
        /* trivial case when the segment has length 0 */
        if (segmentLength == 0) {
            return Double.NaN;
        } else {
            final List<STATE> hiddenStates = fbResult.model().hiddenStates();
            final int numHiddenStates = hiddenStates.size();
            final double[] hiddenStateValues = hiddenStates.stream().mapToDouble(STATE::getScalar).toArray();
            final double[] means = new double[segmentLength];
            double diagonalCorrSum = 0, crossCorrSum;

            /* diagonal contribution to cross correlation sum, and individual target means */
            for (int ti = firstTargetIndex; ti < firstTargetIndex + segmentLength; ti++) {
                final int tif = ti; /* finalized to use in lambda */
                final double[] posteriorProbabilities = IntStream.range(0, numHiddenStates)
                        .mapToDouble(si -> FastMath.exp(fbResult.logProbability(tif, hiddenStates.get(si))))
                        .toArray();
                diagonalCorrSum += IntStream.range(0, numHiddenStates)
                        .mapToDouble(si -> hiddenStateValues[si] * hiddenStateValues[si] * posteriorProbabilities[si])
                        .sum();
                means[ti - firstTargetIndex] = IntStream.range(0, numHiddenStates)
                        .mapToDouble(si -> hiddenStateValues[si] * posteriorProbabilities[si])
                        .sum();
            }

            /* cross terms */
            crossCorrSum = diagonalCorrSum;
            for (int ti0 = firstTargetIndex; ti0 < firstTargetIndex + segmentLength; ti0++) {
                for (int ti1 = ti0 + 1; ti1 < firstTargetIndex + segmentLength; ti1++) {
                    final double distance = Target.calculateDistance(fbResult.positions().get(ti0),
                            fbResult.positions().get(ti1));
                    if (distance > INDEPENDENT_TARGETS_SEPARATION_THRESHOLD) {
                        /* independent state approximation */
                        crossCorrSum += 2 * means[ti0 - firstTargetIndex] * means[ti1 - firstTargetIndex];
                    } else {
                        for (int si0 = 0; si0 < numHiddenStates; si0++) {
                            for (int si1 = 0; si1 < numHiddenStates; si1++) {
                                crossCorrSum += 2 * hiddenStateValues[si0] * hiddenStateValues[si1] *
                                        FastMath.exp(fbResult.logJointProbability(ti0, ti1, hiddenStates.get(si0),
                                                hiddenStates.get(si1)));
                            }
                        }
                    }
                }
            }
            if (segmentLength == 1) { /* return the variance on the only state */
                return diagonalCorrSum - means[0] * means[0];
            } else {
                /* correct for bias by multiplying with segmentLength / (segmentLength - 1) */
                return FastMath.max(diagonalCorrSum / (segmentLength - 1)
                        - crossCorrSum / (segmentLength * (segmentLength - 1)), 0);
            }
        }
    }

    /**
     * See {@link org.broadinstitute.hellbender.utils.GATKProtectedMathUtils#logProbToPhredScore(double, boolean, double, double, double)}
     */
    private static double logProbToPhredScore(final double rawLogProb, final boolean complement) {
        return GATKProtectedMathUtils.logProbToPhredScore(rawLogProb, complement, MAX_QUAL_SCORE, MAX_LOG_PROB,
                PHRED_SCORE_PRECISION);
    }
}
