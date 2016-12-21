package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import com.google.common.collect.Iterators;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.SubtractCoverageComponents;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Discovers potential copy-number segments in input read-counts as the most likely segment sequences
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
        summary = "Find possible locations for rare copy number variation events in germline samples using a HMM",
        oneLineSummary = "Discover possible locations of copy number variation"
)
public final class DiscoverCopyNumberTriStateSegments extends CopyNumberTriStateSegmentCaller {

    private CopyNumberTriStateSegmentRecordWriter outputWriter;

    @Override
    protected void openOutput(final File outputFile, final CopyNumberTriStateHiddenMarkovModel model,
                              final TargetCollection<Target> targets, final ReadCountCollection inputCounts) {
        Utils.nonNull(outputFile);
        Utils.nonNull(model);
        Utils.nonNull(targets);
        Utils.nonNull(inputCounts);
        try {
            outputWriter = new CopyNumberTriStateSegmentRecordWriter(outputFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    @Override
    protected void closeOutput(final File outputFile) {
        try {
            outputWriter.close();
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    @Override
    protected void makeCalls(CopyNumberTriStateHiddenMarkovModel model, TargetCollection<Target> targets, ReadCountCollection inputCounts) {
        final Map<String, List<CopyNumberTriStateSegment>> allSegmentsBySampleName =
                calculateBestPathSegments(model, inputCounts);

        final Iterator<CopyNumberTriStateSegmentRecord> allSegmentRecordsSortedByCoordinates =
                composeATargetSortedSegmentRecordIterator(targets, allSegmentsBySampleName);

        writeSegmentsInOutputFile(allSegmentRecordsSortedByCoordinates);
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
                logProbToPhredScore(logExactProbability, true),
                logProbToPhredScore(logSomeProbability, true),
                logProbToPhredScore(logStartProbability, true),
                logProbToPhredScore(logEndProbability, true),
                logProbToPhredScore(logNeutralProbability, false));
    }
}
