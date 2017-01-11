package org.broadinstitute.hellbender.utils.hmm.segmentation;

import com.google.common.collect.Iterators;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.hmm.interfaces.AlleleMetadataProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.CallStringProducer;
import org.broadinstitute.hellbender.utils.hmm.interfaces.ScalarProducer;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class processes hidden state posteriors, i.e. an instance of {@link ForwardBackwardAlgorithm.Result}),
 * along with best hidden state sequence calls, e.g. the output of
 * {@link org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm}) in order to:
 * <ul>
 *     <li>
 *         Create segments from best paths for each sample along with various quality scores and write the
 *         results to a SEG file
 *     </li>
 *     <li>
 *         Create variant contexts from obtained or externally provided segments and writer the results
 *         to a VCF file
 *     </li>
 * </ul>
 *
 * <p>
 * The user must identify the reference state in order to create variant contexts. All other hidden states
 * are construed as alternative alleles.
 * </p>
 *
 * </p>
 * Each sample may have its own set of targets. The targets must be lexicographically
 * sorted (according to {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}).
 * A {@link UserException.BadInput} exception will be thrown.
 * <p>
 *
 * @param <T> target type; must extend {@link Target}
 *
 * @param <S> hidden state type; must implement {@link AlleleMetadataProducer} and {@link CallStringProducer}.
 *            For example, see the enum type {@link org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState}
 *            as a typical use case.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class HiddenMarkovModelPostProcessor<D, S extends AlleleMetadataProducer & CallStringProducer & ScalarProducer,
        T extends Target> {

    private final Logger logger = LogManager.getLogger(HiddenMarkovModelPostProcessor.class);

    /**
     * VCF header keys
     */
    public static final String DISCOVERY_KEY = "DSCVR";
    public static final String NUMBER_OF_POOLED_TARGETS_KEY = "NTARGETS";
    public static final String SOME_QUALITY_KEY = "SQ";
    public static final String START_QUALITY_KEY = "LQ";
    public static final String END_QUALITY_KEY = "RQ";
    public static final String NUMBER_OF_SAMPLE_SPECIFIC_TARGETS_KEY = "NT";
    public static final String DISCOVERY_TRUE = "Y";
    public static final String DISCOVERY_FALSE = "N";

    /**
     * Genotype qualities will be capped at this value
     */
    public static final int MAX_GQ = 99;

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
    public static final double MAX_QUAL_SCORE = 199.99999999999;

    /**
     * Precision of Phred-scale scores
     */
    public static final double PHRED_SCORE_PRECISION = 0.1;

    /**
     * Calculating variance posterior on long segments is time consuming; we assume hidden states become
     * independent if their separation is larger than the following value
     */
    public static final double INDEPENDENT_TARGETS_SEPARATION_THRESHOLD = 1000;

    /**
     * Threshold used to determine best way to calculate log(1- exp(a))
     * based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
     */
    private static final double LN_1_M_EXP_THRESHOLD = - Math.log(2);

    /**
     * Cached value of ln(10)^-1
     */
    private static final double INV_LN_10 = 1.0 / Math.log(10);

    private final List<S> allStates;
    private final S referenceState;
    private final List<S> alternativeStates;

    private final List<Allele> allAlleles;
    private final List<Allele> alternativeAlleles;

    private final List<String> sampleNames;
    private final int numSamples;
    private final List<TargetCollection<T>> sampleTargets;
    private final List<ForwardBackwardAlgorithm.Result<D, T, S>> sampleForwardBackwardResults;
    private final List<List<S>> sampleBestPaths;

    /**
     * An intermediate result of segment calling -- will be used for VCF creation if an external
     * segments file is not given
     */
    private Map<String, List<HiddenStateSegment<S, T>>> allSegmentsBySampleName = null;

    /**
     * This "master" target collection is the lexicographically sorted union of all targets from all samples
     */
    private final TargetCollection<T> masterTargetCollection;

    /**
     * Public constructor.
     *
     * @param sampleNames list of sample names
     * @param sampleTargets list of target collection for each sample
     * @param sampleForwardBackwardResults list of forward-backward result for each sample
     * @param sampleBestPaths list of best hidden state sequence call for each sample
     * @param referenceState the reference state
     */
    public HiddenMarkovModelPostProcessor(
            @Nonnull final List<String> sampleNames,
            @Nonnull final List<TargetCollection<T>> sampleTargets,
            @Nonnull final List<ForwardBackwardAlgorithm.Result<D, T, S>> sampleForwardBackwardResults,
            @Nonnull final List<List<S>> sampleBestPaths,
            @Nonnull final S referenceState) {

        /* basic consistency checks on the input data */
        numSamples = sampleNames.size();
        Utils.validateArg(numSamples > 0, "Number of samples must be > 0");
        Utils.validateArg(new HashSet<>(sampleNames).size() == numSamples, "Sample names must be unique");
        Utils.validateArg(sampleTargets.size() == numSamples, "List of targets per sample must have the" +
                " same length as the list of sample names");
        Utils.validateArg(sampleForwardBackwardResults.size() == numSamples, "List of forward-backward" +
                " result per sample must have the same length as the list of sample names");
        Utils.validateArg(sampleBestPaths.size() == numSamples, "List of best paths per sample must have" +
                " the same length as the list of sample names");
        Utils.validateArg(IntStream.range(0, numSamples)
                .allMatch(si -> sampleTargets.get(si).targetCount() == sampleBestPaths.get(si).size()),
                "Some of the provided sample best paths have different lengths than their corresponding target list");

        /* get the hidden states from the first sample and ensure all samples have the same hidden states */
        final List<S> allStatesInferred = sampleForwardBackwardResults.get(0).model().hiddenStates();
        final Set<S> allStatesInferredSet = new HashSet<>(allStatesInferred);
        Utils.validateArg(IntStream.range(1, numSamples)
                .allMatch(si -> new HashSet<>(sampleForwardBackwardResults.get(si).model().hiddenStates())
                        .equals(allStatesInferredSet)), "All samples must have the same set of hidden states");
        this.referenceState = referenceState;
        Utils.validateArg(allStatesInferredSet.contains(referenceState), "The provided reference state is not in" +
                " the HMM hidden states");
        alternativeStates = new ArrayList<>(allStatesInferred);
        alternativeStates.remove(referenceState);

        /**
         * The reference state comes first; this is required by methods used in
         * {@link #composeVariantContext(GenotypingSegment, String)}
         */
        allStates = new ArrayList<>(allStatesInferred.size());
        allStates.add(referenceState);
        allStates.addAll(alternativeStates);

        /* TODO in the future, we may want to replace lexicographical order with natural order */
        /* assert that the target collection of each sample is lexicographically sorted */
        IntStream.range(0, numSamples)
                .forEach(si -> {
                    final List<T> sampleTargetList = sampleTargets.get(si).targets();
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

        /* the list of alleles must be in the same order as the corresponding lists of hidden states */
        allAlleles = allStates.stream()
                .map(AlleleMetadataProducer::toAllele)
                .collect(Collectors.toList());
        alternativeAlleles = alternativeStates.stream()
                .map(AlleleMetadataProducer::toAllele)
                .collect(Collectors.toList());

        /* pool targets from all samples and sort them lexicographically to make a master target collection */
        masterTargetCollection = new HashedListTargetCollection<>(sampleTargets.stream()
                .map(TargetCollection::targets)
                .flatMap(List::stream)
                .distinct()
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList()));
    }

    /***********************
     * main public methods *
     ***********************/

    /**
     * Write segments to a table writer
     *
     * @param writer an instance of {@link HiddenStateSegmentRecordWriter}
     */
    public void writeSegmentsToTableWriter(@Nonnull final HiddenStateSegmentRecordWriter<S, T> writer) {
        if (allSegmentsBySampleName == null) {
            performSegmentation();
        }
        final Iterator<HiddenStateSegmentRecord<S, T>> allSegmentRecordsSortedByCoordinates =
                composeTargetSortedSegmentRecordIterator(masterTargetCollection, allSegmentsBySampleName);
        logger.info("Writing segments to file...");
        writeSegments(allSegmentRecordsSortedByCoordinates, writer);
    }

    /**
     * Create the VCF file based on segmentation of the provided input data
     *
     * @param outputWriter variant context writer
     * @param variantPrefix a prefix to attach to variant IDs
     * @param commandLine (optional) command line used to generate the data handed to this class
     */
    public void writeVariantsToVCFWriter(@Nonnull final VariantContextWriter outputWriter,
                                         @Nonnull final String variantPrefix,
                                         @Nullable final String commandLine) {
        if (allSegmentsBySampleName == null) {
            performSegmentation();
        }
        logger.info("Generating genotyping segments...");
        final List<GenotypingSegment> genotypingSegments = composeGenotypingSegments(allSegmentsBySampleName);
        composeVariantContextAndWrite(genotypingSegments, outputWriter, variantPrefix, commandLine);
    }

    /**
     * Create the VCF file based on pre-determined segments. The forward-backward result provided to the
     * constructor is used for variant context generation.
     *
     * @param segmentsFile a segment file
     * @param outputWriter variant context writer
     * @param variantPrefix a prefix to attach to variant IDs
     * @param commandLine (optional) command line used to generate the data handled to this class
     */
    public void writeVariantsToVCFWriterOnGivenSegments(@Nonnull final File segmentsFile,
                                                        @Nonnull final Function<String, S> inverseCallStringFactory,
                                                        @Nonnull final VariantContextWriter outputWriter,
                                                        @Nonnull final String variantPrefix,
                                                        @Nullable final String commandLine) {
        logger.info("Composing genotyping segments from the provided segments file...");
        final List<GenotypingSegment> genotypingSegments =
                composeGenotypingSegmentsFromSegmentsFile(segmentsFile, inverseCallStringFactory);
        composeVariantContextAndWrite(genotypingSegments, outputWriter, variantPrefix, commandLine);
    }


    /********************************
     * segmentation-related methods *
     ********************************/

    /**
     * Compose segments for all samples
     */
    private void performSegmentation() {
        logger.info("Composing segments...");
        allSegmentsBySampleName = calculateBestPathSegments();
    }

    /**
     * Takes the provided forward-backward result and the best hidden state chain for each sample and
     * generates segments
     *
     * @return a map from sample names to a list of segments
     */
    private Map<String, List<HiddenStateSegment<S, T>>> calculateBestPathSegments() {
        final Map<String, List<HiddenStateSegment<S, T>>> allSegments = new LinkedHashMap<>(sampleNames.size());
        IntStream.range(0, numSamples)
                .forEach(sampleIndex -> {
                    final String sampleName = sampleNames.get(sampleIndex);
                    final TargetCollection<T> targets = sampleTargets.get(sampleIndex);
                    final List<S> bestPath = sampleBestPaths.get(sampleIndex);
                    final ForwardBackwardAlgorithm.Result<D, T, S> fbResult = sampleForwardBackwardResults.get(sampleIndex);
                    final List<Pair<IndexRange, S>> bestPathTargetIndexRanges =
                            condenseBestPathIntoTargetIndexAndStatePairs(bestPath, targets);
                    final List<HiddenStateSegment<S, T>> bestPathSegmentList =
                            composeSegments(fbResult, bestPathTargetIndexRanges);
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
    private List<Pair<IndexRange, S>> condenseBestPathIntoTargetIndexAndStatePairs(final List<S> bestPath,
                                                                                   final TargetCollection<T> targets) {
        if (bestPath.isEmpty()) {
            return Collections.emptyList();
        }

        /* estimate the number of segments */
        final long numChangePoints = IntStream.range(0, bestPath.size() - 1)
                .filter(i -> !bestPath.get(i).equals(bestPath.get(i + 1)))
                .count();
        final List<Pair<IndexRange, S>> result = new ArrayList<>((int) numChangePoints);

        int currentStartIndex = 0; // contains the start index of the segment being traversed.
        final ListIterator<S> pathIterator = bestPath.listIterator();
        final ListIterator<T> targetIterator = targets.targets().listIterator();
        String currentContig = targetIterator.next().getContig();
        S currentState = pathIterator.next();
        while (pathIterator.hasNext()) {
            final S nextState = pathIterator.next();
            final String nextContig = targetIterator.next().getContig();
            final boolean contigChanged = !currentContig.equals(nextContig);
            final boolean stateChanged = (nextState != currentState);
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
    private List<HiddenStateSegment<S, T>> composeSegments(
            final ForwardBackwardAlgorithm.Result<D, T, S> fbResult,
            final List<Pair<IndexRange, S>> bestPathSegments) {
        return bestPathSegments.stream()
                .map(ir -> composeHiddenStateSegmentFromTargetIndexRange(fbResult, ir.getLeft(), ir.getRight()))
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
    private HiddenStateSegment<S, T> composeHiddenStateSegmentFromTargetIndexRange(
            final ForwardBackwardAlgorithm.Result<D, T, S> fbResult,
            final IndexRange targetIndexRange,
            final S call) {

        final int segmentLength = targetIndexRange.size();
        final double mean = calculateSegmentMean(targetIndexRange.from, segmentLength, fbResult);
        final double stdDev = FastMath.sqrt(calculateSegmentVariance(targetIndexRange.from, segmentLength, fbResult));

        final List<T> targets = fbResult.positions().subList(targetIndexRange.from, targetIndexRange.to);
        final double logExactProbability = fbResult.logProbability(targetIndexRange.from, targetIndexRange.to, call);
        final double logSomeProbability = logSomeProbability(targetIndexRange.from, segmentLength, call, fbResult);
        final double logStartProbability = logStartProbability(targetIndexRange.from, call, fbResult);
        final double logEndProbability = logEndProbability(targetIndexRange.to - 1, call, fbResult);
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
    private Iterator<HiddenStateSegmentRecord<S, T>> composeTargetSortedSegmentRecordIterator(
            final TargetCollection<T> masterTargetCollection,
            final Map<String, List<HiddenStateSegment<S, T>>> allSegments) {
        final List<Iterator<HiddenStateSegmentRecord<S, T>>> allSegmentRecordIterators =
                allSegments.entrySet().stream()
                        .map(e -> e.getValue().stream()
                                .map(s -> new HiddenStateSegmentRecord<>(e.getKey(), s))
                                .iterator())
                        .collect(Collectors.toList());

        final Comparator<HiddenStateSegmentRecord<S, T>> recordComparator = (o1, o2) -> {
            final HiddenStateSegment<S, T> s1 = o1.getSegment();
            final HiddenStateSegment<S, T> s2 = o2.getSegment();
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
    private void writeSegments(final Iterator<HiddenStateSegmentRecord<S, T>> allSegmentsSorted,
                               final HiddenStateSegmentRecordWriter<S, T> writer) {
        try {
            long count = 0;
            final Set<String> samples = new HashSet<>();
            while (allSegmentsSorted.hasNext()) {
                final HiddenStateSegmentRecord<S, T> record = allSegmentsSorted.next();
                samples.add(record.getSampleName());
                writer.writeRecord(record);
                count++;
            }
            logger.info(String.format("Found a total of %d segments across %d samples", count, samples.size()));

        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile("problems writing segments", ex);
        }
    }

    /***********************
     * VCF-related methods *
     ***********************/

    /**
     * Compose genotyping segments from the internally performed segmentation of the provided HMM results
     *
     * @param allSegmentsBySampleName a map from each sample name to its list of segments
     * @return a list of genotyping segments
     */
    private List<GenotypingSegment> composeGenotypingSegments(
            @Nonnull final Map<String, List<HiddenStateSegment<S, T>>> allSegmentsBySampleName) {

        /* create lexicographically sorted segment records with non-reference calls */
        final List<HiddenStateSegmentRecord<S, T>> records = allSegmentsBySampleName.entrySet().stream()
                .flatMap(entry -> entry.getValue().stream()
                        .map(seg -> new HiddenStateSegmentRecord<>(entry.getKey(), seg)))
                .filter(rec -> !rec.getSegment().getCall().equals(referenceState))
                .sorted(Comparator.comparing(HiddenStateSegmentRecord::getSegment,
                        IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR))
                .collect(Collectors.toList());

        final List<GenotypingSegment> result = new ArrayList<>();
        for (final HiddenStateSegmentRecord<S, T> record : records) {
            if (result.isEmpty()) {
                result.add(new GenotypingSegment(record, masterTargetCollection));
            } else {
                final GenotypingSegment lastSegment = result.get(result.size() - 1);
                if (lastSegment.getInterval().equals(record.getSegment().getInterval())) {
                    lastSegment.addSample(record.getSampleName());
                } else {
                    result.add(new GenotypingSegment(record, masterTargetCollection));
                }
            }
        }
        return result;
    }

    /**
     * Compose genotyping segments from a previous segmentation
     *
     * @param segmentsFile a segments file
     * @param inverseCallStringFactory a map from call string to hidden state
     * @throws UserException.BadInput if the sample names in the provided segments file does not match those on which
     * HMM output is given
     * @return a list of genotyping segments
     */
    private List<GenotypingSegment> composeGenotypingSegmentsFromSegmentsFile(
            @Nonnull final File segmentsFile,
            @Nonnull final Function<String, S> inverseCallStringFactory) {
        try (final HiddenStateSegmentRecordReader<S, T> reader =
                     new HiddenStateSegmentRecordReader<>(segmentsFile, inverseCallStringFactory)) {
            final List<HiddenStateSegmentRecord<S, T>> allRecords = reader.toList();
            final List<HiddenStateSegmentRecord<S, T>> nonRefAndSortedRecords = allRecords.stream()
                .filter(s -> !s.getSegment().getCall().equals(referenceState))
                        .sorted(Comparator.comparing(
                                HiddenStateSegmentRecord::getSegment,
                                IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR))
                        .collect(Collectors.toList());
            /* ensure that the segments have the same sample names as the one provided to the constructor */
            final Set<String> sampleNamesFromSegmentsFile = allRecords.stream()
                    .map(HiddenStateSegmentRecord::getSampleName)
                    .collect(Collectors.toSet());
            if (!new HashSet<>(sampleNames).equals(sampleNamesFromSegmentsFile)) {
                throw new UserException.BadInput("The sample names in the provided segments file does not match" +
                        " those on which HMM output is given");
            }
            final List<GenotypingSegment> result = new ArrayList<>();
            for (final HiddenStateSegmentRecord<S, T> record : nonRefAndSortedRecords) {
                if (result.isEmpty()) {
                    result.add(new GenotypingSegment(record, masterTargetCollection));
                } else {
                    final GenotypingSegment lastSegment = result.get(result.size() - 1);
                    if (lastSegment.getInterval().equals(record.getSegment().getInterval())) {
                        lastSegment.addSample(record.getSampleName());
                    } else {
                        result.add(new GenotypingSegment(record, masterTargetCollection));
                    }
                }
            }
            return result;
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(segmentsFile, ex);
        }
    }

    /**
     * Composes a variant context on a given genotyping segment
     *
     * @param segment a genotyping segment
     * @param variantPrefix a prefix to attach to the variant ID
     * @return a variant context
     */
    private VariantContext composeVariantContext(@Nonnull final GenotypingSegment segment,
                                                 @Nonnull final String variantPrefix) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        builder.alleles(allAlleles);
        builder.chr(segment.getContig());
        builder.start(segment.getStart());
        builder.stop(segment.getEnd());
        builder.id(String.format(variantPrefix + "_%s_%d_%d", segment.getContig(),
                segment.getStart(), segment.getEnd()));

        final List<Genotype> genotypes = IntStream.range(0, sampleNames.size())
                .mapToObj(sampleIndex -> {
                    final String sample = sampleNames.get(sampleIndex);
                    final ForwardBackwardAlgorithm.Result<D, T, S> fbResult =
                            sampleForwardBackwardResults.get(sampleIndex);

                    // Different samples could have different lists of targets; we must find the index
                    // range of the segment for each sample separately
                    final IndexRange sampleIndexRange = sampleTargets.get(sampleIndex).indexRange(segment);

                    final double[] log10GP = calculateLog10GP(sampleIndexRange, fbResult);
                    final double[] SQ = calculateSQ(sampleIndexRange, fbResult);
                    final double[] LQ = calculateLQ(sampleIndexRange, fbResult);
                    final double[] RQ = calculateRQ(sampleIndexRange, fbResult);
                    final GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(log10GP);
                    final int[] PL = likelihoods.getAsPLs();
                    final int GQ = calculateGQ(PL);
                    final int genotypeCall = MathUtils.maxElementIndex(log10GP);

                    final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                    genotypeBuilder.PL(PL);
                    genotypeBuilder.alleles(Collections.singletonList(allAlleles.get(genotypeCall)));
                    genotypeBuilder.attribute(DISCOVERY_KEY, segment.containingSamples.contains(sample) ?
                            DISCOVERY_TRUE : DISCOVERY_FALSE);
                    genotypeBuilder.attribute(SOME_QUALITY_KEY, SQ);
                    genotypeBuilder.attribute(START_QUALITY_KEY, LQ);
                    genotypeBuilder.attribute(END_QUALITY_KEY, RQ);
                    genotypeBuilder.attribute(NUMBER_OF_SAMPLE_SPECIFIC_TARGETS_KEY, sampleIndexRange.size());
                    genotypeBuilder.GQ(GQ);

                    return genotypeBuilder.make();
                }).collect(Collectors.toList());

        final int alleleNumber = genotypes.size();
        final int[] altAlleleCounts = alternativeAlleles.stream()
                .mapToInt(allele -> (int) genotypes.stream()
                        .filter(g -> g.getAllele(0).equals(allele))
                        .count())
                .toArray();
        final double[] altAlleleFrequencies = Arrays.stream(altAlleleCounts)
                .mapToDouble(count -> count / (double) alleleNumber)
                .toArray();

        builder.attribute(VCFConstants.END_KEY, segment.getEnd());
        builder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, altAlleleFrequencies);
        builder.attribute(VCFConstants.ALLELE_COUNT_KEY, altAlleleCounts);
        builder.attribute(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
        builder.attribute(NUMBER_OF_POOLED_TARGETS_KEY, segment.getTargetIndexes().size());
        builder.genotypes(genotypes);

        return builder.make();
    }

    /**
     * Composes the VCF header
     *
     * @param commandLine optional command line used to generate the underlying data
     * @return
     */
    private VCFHeader composeHeader(@Nullable final String commandLine) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), sampleNames);

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add description of alternative alleles */
        alternativeStates.stream().forEach(s -> s.addHeaderLineTo(result));

        /* add command line */
        if (commandLine != null) {
            result.addMetaDataLine(new VCFHeaderLine("command", commandLine));
        }

        /* header lines related to formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.Integer, "Genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_PL_KEY, VCFHeaderLineCount.G,
                VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"));
        result.addMetaDataLine(new VCFFormatHeaderLine(DISCOVERY_KEY, 1,
                VCFHeaderLineType.Character, "Y if this segment was discovered in that sample, N otherwise"));
        result.addMetaDataLine(new VCFFormatHeaderLine(SOME_QUALITY_KEY, VCFHeaderLineCount.A,
                VCFHeaderLineType.Float, "Quality of the region to contain at least one target with that call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(START_QUALITY_KEY, VCFHeaderLineCount.A,
                VCFHeaderLineType.Float, "Quality of the segment to start exactly at that location for alternative allele calls only"));
        result.addMetaDataLine(new VCFFormatHeaderLine(END_QUALITY_KEY, VCFHeaderLineCount.A,
                VCFHeaderLineType.Float, "Quality of the segment to end exactly at that location for alternative allele calls only"));
        result.addMetaDataLine(new VCFFormatHeaderLine(NUMBER_OF_SAMPLE_SPECIFIC_TARGETS_KEY, VCFHeaderLineCount.A,
                VCFHeaderLineType.Integer, "Number of targets enclosed in this variant for the sample"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1,
                VCFHeaderLineType.Integer, "Genotype call quality as the difference between the best and second best PL"));

        /* informational header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY, VCFHeaderLineCount.A,
                VCFHeaderLineType.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A,
                VCFHeaderLineType.Float, "Allele Frequency, for each ALT allele, in the same order as listed"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY, 1,
                VCFHeaderLineType.Integer, "Total number of alleles in called genotypes"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of this variant"));
        result.addMetaDataLine(new VCFInfoHeaderLine(NUMBER_OF_POOLED_TARGETS_KEY, 1,
                VCFHeaderLineType.Integer, "Number of targets enclosed in this variant coordinates pooled from" +
                " all samples (some samples may have fewer)"));

        return result;
    }

    /**
     * For each segment in genotypingSegments, compose the variant context and write to outputWriter
     *
     * @param genotypingSegments a list of genotyping segments
     * @param outputWriter a VCF writer
     * @param variantPrefix a prefix for composing variant IDs
     * @param commandLine (optional) command line used to generated the data
     */
    private void composeVariantContextAndWrite(@Nonnull final List<GenotypingSegment> genotypingSegments,
                                               @Nonnull final VariantContextWriter outputWriter,
                                               @Nonnull final String variantPrefix,
                                               @Nullable final String commandLine) {
        outputWriter.writeHeader(composeHeader(commandLine));
        for (final GenotypingSegment segment : genotypingSegments) {
            final VariantContext variant = composeVariantContext(segment, variantPrefix);
            outputWriter.add(variant);
        }
    }


    /****************************************************************************************
     * methods related to calculating various quality scores, probabilities, and statistics *
     ****************************************************************************************/

    /**
     * Calculates genotype quality from Phred scale likelihoods
     *
     * @param PL Phred scale likelihoods array
     * @return
     */
    private int calculateGQ(final int[] PL) {
        int best = PL[0];
        int secondBest = Integer.MAX_VALUE;
        for (int i = 1; i < PL.length; i++) {
            final int value = PL[i];
            if (value <= best) {
                secondBest = best;
                best = value;
            } else if (value < secondBest) {
                secondBest = value;
            }
        }
        return Math.min(secondBest - best, MAX_GQ);
    }

    /**
     * Calculates variant end quality for a genotyping segment
     *
     * @param targetIndexes the target index range for the segment for the -specific- sample
     * @param fbResult forward-backward result for the -specific- sample
     * @return
     */
    private double[] calculateRQ(final IndexRange targetIndexes,
                                 final ForwardBackwardAlgorithm.Result<D, T, S> fbResult) {
        return alternativeStates.stream()
                .mapToDouble(state ->
                        logProbToPhredScore(logEndProbability(targetIndexes.to - 1, state, fbResult), true))
                .toArray();
    }

    /**
     * Calculates variant start quality for a genotyping segment
     *
     * @param targetIndexes the target index range for the segment for the -specific- sample
     * @param fbResult forward-backward result for the -specific- sample
     * @return
     */
    private double[] calculateLQ(final IndexRange targetIndexes,
                                 final ForwardBackwardAlgorithm.Result<D, T, S> fbResult) {
        return alternativeStates.stream()
                .mapToDouble(state ->
                        logProbToPhredScore(logStartProbability(targetIndexes.from, state, fbResult), true))
                .toArray();
    }

    /**
     * Calculates variant some quality for a segment
     *
     * @param targetIndexes the target index range for the segment for the -specific- sample
     * @param fbResult forward-backward result for the -specific- sample
     * @return
     */
    private double[] calculateSQ(final IndexRange targetIndexes,
                                 final ForwardBackwardAlgorithm.Result<D, T, S> fbResult) {
        return alternativeStates.stream()
                .mapToDouble(state ->
                        logProbToPhredScore(logSomeProbability(targetIndexes.from, targetIndexes.size(), state,
                                fbResult), true))
                .toArray();
    }

    /**
     * Calculates genotyping probability in log_10 scale for a segment
     *
     * @param targetIndexes the target index range for the segment for the -specific- sample
     * @param fbResult forward-backward result for the -specific- sample
     * @return
     */
    private double[] calculateLog10GP(final IndexRange targetIndexes,
                                      final ForwardBackwardAlgorithm.Result<D, T, S> fbResult) {
        return allStates.stream()
                .mapToDouble(state -> fbResult.logProbability(targetIndexes.from, targetIndexes.to, state) * INV_LN_10)
                .map(HiddenMarkovModelPostProcessor::roundPhred)
                .toArray();
    }

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
            return logProbComplement(logOtherStates);
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
     *      E[\sigma^2] = (1/T) \sum_{t} E[S_t^2] - (1/T^2) \sum_{t1,t2} E[S_t1 S_t2]
     *
     * Here, T = {@code segmentLength}
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
                            fbResult.positions().get(ti1), Double.POSITIVE_INFINITY);
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
     * Calculates the complement of a log probability.
     *
     * <p>
     *     With complement of {@code x} we mean: {@code log(1-log(x))}.
     * </p>
     * @param x the input log probability.
     * @return {@code log(1-log(x))}
     */
    public static double logProbComplement(final double x) {
        return x >= LN_1_M_EXP_THRESHOLD
                ? Math.log(-Math.expm1(x))
                : Math.log1p(-Math.exp(x));
    }

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
     *                    or its complement ({@code true)}.
     * @return a values between 0 and {@link #MAX_QUAL_SCORE}.
     * @throws GATKException if {@code rawLogProb} is larger than {@link #MAX_LOG_PROB}.
     */
    public static double logProbToPhredScore(final double rawLogProb, final boolean complement) {
        if (rawLogProb > MAX_LOG_PROB) {
            throw new GATKException(String.format("numerical instability problem: the log-probability is too" +
                    " large: %g > 0.0 (with maximum tolerance %g)", rawLogProb, MAX_LOG_PROB));
        }
        // make sure that the probability is less than 1 in linear scale. there are cases where
        // log probability exceeds 0.0 due to floating point errors.
        final double logProbEqOrLessThan0 = Math.min(0.0, rawLogProb);

        // Accurate way to calculate log(1-exp(a))
        // based on https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
        final double finalLogProb = complement
                ? logProbComplement(logProbEqOrLessThan0)
                : logProbEqOrLessThan0;

        final double absoluteQualScore = QualityUtils.phredScaleLog10ErrorRate(finalLogProb * INV_LN_10);
        final double exactValue = Math.min(MAX_QUAL_SCORE, absoluteQualScore);
        // We round the value to the required precession.
        return roundPhred(exactValue);
    }

    /**
     * Round a Phred scaled score to precision {@link #PHRED_SCORE_PRECISION}
     *
     * @param value Phred score
     * @return rounded Phred score
     */
    public static double roundPhred(final double value) {
        return Math.round(value / PHRED_SCORE_PRECISION) * PHRED_SCORE_PRECISION;
    }

    /**
     * Use to collect information about segment to be genotyped.
     */
    private final class GenotypingSegment implements Locatable {
        /**
         * Genomic interval
         */
        private final SimpleInterval interval;

        /**
         * Index range of the segment in the master target collection
         */
        private final IndexRange targetIndexes;

        /**
         * Samples that share this segment
         */
        private final Set<String> containingSamples;

        public GenotypingSegment(final SimpleInterval interval, final IndexRange targetIndexes) {
            this.interval = Utils.nonNull(interval);
            this.targetIndexes = Utils.nonNull(targetIndexes);
            this.containingSamples = new HashSet<>();
        }

        public GenotypingSegment(final HiddenStateSegmentRecord<S, T> record,
                                 final TargetCollection<T> masterTargetCollection) {
            this(record.getSegment().getInterval(), masterTargetCollection.indexRange(record.getSegment()));
            addSample(record.getSampleName());
        }

        public void addSample(final String name) {
            this.containingSamples.add(name);
        }

        public IndexRange getTargetIndexes() {
            return targetIndexes;
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        public int getTargetCount() {
            return targetIndexes.size();
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }
    }
}
