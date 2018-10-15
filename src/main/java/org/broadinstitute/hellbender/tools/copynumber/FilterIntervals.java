package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given annotated intervals output by {@link AnnotateIntervals} and/or counts collected on those intervals output
 * by {@link CollectReadCounts}, outputs a filtered Picard interval list.  Parameters for filtering based on the
 * annotations and counts can be adjusted.  Annotation-based filters will be applied first, followed by count-based
 * filters.  The result may be passed via -L to other tools (e.g., {@link DetermineGermlineContigPloidy} and
 * {@link GermlineCNVCaller}) to mask intervals from analysis.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
           Intervals to be filtered.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *         A blacklist of regions in which intervals should always be filtered (regardless of other annotation-based
 *         or count-based filters) may also be provided via -XL; this can be used to filter pseudoautosomal regions
 *         (PARs), for example.
 *     </li>
 *     <li>
 *         (Optional) Annotated-intervals file from {@link AnnotateIntervals}.
 *         Must contain the intervals to be filtered as a subset.  Must be provided if no counts files are provided.
 *     </li>
 *     <li>
 *         (Optional) Counts files (TSV or HDF5 output of {@link CollectReadCounts}).
 *         Must contain the intervals to be filtered as a subset.  Must be provided if no annotated-intervals file
 *         is provided.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Filtered Picard interval-list file.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk FilterIntervals \
 *          -L intervals.interval_list \
 *          -XL blacklist_intervals.interval_list \
 *          -I sample_1.counts.hdf5 \
 *          -I sample_2.counts.hdf5 \
 *          ... \
 *          --annotated-intervals annotated_intervals.tsv \
 *          -O filtered_intervals.interval_list
 * </pre>
 *
 * <pre>
 *     gatk FilterIntervals \
 *          -L intervals.interval_list \
 *          --annotated-intervals annotated_intervals.tsv \
 *          -O filtered_intervals.interval_list
 * </pre>
 *
 * <pre>
 *     gatk FilterIntervals \
 *          -L intervals.interval_list \
 *          -I sample_1.counts.hdf5 \
 *          -I sample_2.counts.hdf5 \
 *          ... \
 *          -O filtered_intervals.interval_list
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Filters intervals based on annotations and/or count statistics",
        oneLineSummary = "Filters intervals based on annotations and/or count statistics",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class FilterIntervals extends CommandLineProgram {
    public static final String MINIMUM_GC_CONTENT_LONG_NAME = "minimum-gc-content";
    public static final String MAXIMUM_GC_CONTENT_LONG_NAME = "maximum-gc-content";
    public static final String MINIMUM_MAPPABILITY_LONG_NAME = "minimum-mappability";
    public static final String MAXIMUM_MAPPABILITY_LONG_NAME = "maximum-mappability";
    public static final String MINIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME = "minimum-segmental-duplication-content";
    public static final String MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME = "maximum-segmental-duplication-content";
    public static final String LOW_COUNT_FILTER_COUNT_THRESHOLD_LONG_NAME = "low-count-filter-count-threshold";
    public static final String LOW_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME = "low-count-filter-percentage-of-samples";
    public static final String EXTREME_COUNT_FILTER_MINIMUM_PERCENTILE_LONG_NAME = "extreme-count-filter-minimum-percentile";
    public static final String EXTREME_COUNT_FILTER_MAXIMUM_PERCENTILE_LONG_NAME = "extreme-count-filter-maximum-percentile";
    public static final String EXTREME_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME = "extreme-count-filter-percentage-of-samples";

    @Argument(
            doc = "Input file containing annotations for genomic intervals (output of AnnotateIntervals).  " +
                    "All intervals specified via -L must be contained.  " +
                    "Must be provided if no counts files are provided.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAnnotatedIntervalsFile = null;

    @Argument(
            doc = "Input TSV or HDF5 files containing integer read counts in genomic intervals (output of CollectReadCounts).  " +
                    "All intervals specified via -L must be contained.  " +
                    "Must be provided if no annotated-intervals file is provided.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = true
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Output Picard interval-list file containing the filtered intervals.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFilteredIntervalsFile;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection
            = new OptionalIntervalArgumentCollection();

    @Argument(
            doc = "Minimum allowed value for GC-content annotation (inclusive).",
            fullName = MINIMUM_GC_CONTENT_LONG_NAME,
            optional = true
    )
    private double minimumGCContent = 0.1;

    @Argument(
            doc = "Maximum allowed value for GC-content annotation (inclusive).",
            fullName = MAXIMUM_GC_CONTENT_LONG_NAME,
            optional = true
    )
    private double maximumGCContent = 0.9;

    @Argument(
            doc = "Minimum allowed value for mappability annotation (inclusive).",
            fullName = MINIMUM_MAPPABILITY_LONG_NAME,
            optional = true
    )
    private double minimumMappability = 0.9;

    @Argument(
            doc = "Maximum allowed value for mappability annotation (inclusive).",
            fullName = MAXIMUM_MAPPABILITY_LONG_NAME,
            optional = true
    )
    private double maximumMappability = 1.;

    @Argument(
            doc = "Minimum allowed value for segmental-duplication-content annotation (inclusive).",
            fullName = MINIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME,
            optional = true
    )
    private double minimumSegmentalDuplicationContent = 0.;

    @Argument(
            doc = "Maximum allowed value for segmental-duplication-content annotation (inclusive).",
            fullName = MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT_LONG_NAME,
            optional = true
    )
    private double maximumSegmentalDuplicationContent = 0.5;

    @Argument(
            doc = "Count-threshold parameter for the low-count filter.  Intervals with a count " +
                    "strictly less than this threshold in a percentage of samples strictly greater than " +
                    LOW_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME + " will be filtered out.  " +
                    "(This is the first count-based filter applied.)",
            fullName = LOW_COUNT_FILTER_COUNT_THRESHOLD_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int lowCountFilterCountThreshold = 5;

    @Argument(
            doc = "Percentage-of-samples parameter for the low-count filter.  Intervals with a count " +
                    "strictly less than " + LOW_COUNT_FILTER_COUNT_THRESHOLD_LONG_NAME +
                    " in a percentage of samples strictly greater than this will be filtered out.  " +
                    "(This is the first count-based filter applied.)",
            fullName = LOW_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double lowCountFilterPercentageOfSamples = 90.;

    @Argument(
            doc = "Minimum-percentile parameter for the extreme-count filter.  Intervals with a count " +
                    "that has a percentile strictly less than this in a percentage of samples strictly greater than " +
                    EXTREME_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME + " will be filtered out.  " +
                    "(This is the second count-based filter applied.)",
            fullName = EXTREME_COUNT_FILTER_MINIMUM_PERCENTILE_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double extremeCountFilterMinimumPercentile = 1.;

    @Argument(
            doc = "Maximum-percentile parameter for the extreme-count filter.  Intervals with a count " +
                    "that has a percentile strictly greater than this in a percentage of samples strictly greater than " +
                    EXTREME_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME + " will be filtered out.  " +
                    "(This is the second count-based filter applied.)",
            fullName = EXTREME_COUNT_FILTER_MAXIMUM_PERCENTILE_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double extremeCountFilterMaximumPercentile = 99.;

    @Argument(
            doc = "Percentage-of-samples parameter for the extreme-count filter.  Intervals with a count " +
                    "that has a percentile outside of [" + EXTREME_COUNT_FILTER_MINIMUM_PERCENTILE_LONG_NAME + ", " +
                    EXTREME_COUNT_FILTER_MAXIMUM_PERCENTILE_LONG_NAME + "] in a percentage of samples strictly greater than " +
                    "this will be filtered out.  (This is the second count-based filter applied.)",
            fullName = EXTREME_COUNT_FILTER_PERCENTAGE_OF_SAMPLES_LONG_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double extremeCountFilterPercentageOfSamples = 90.;

    private SimpleIntervalCollection specifiedIntervals;
    private AnnotatedIntervalCollection annotatedIntervals;

    @Override
    public Object doWork() {
        validateFilesAndResolveIntervals();
        final SimpleIntervalCollection filteredIntervals = filterIntervals();
        logger.info(String.format("Writing filtered intervals to %s...", outputFilteredIntervalsFile));
        final IntervalList filteredIntervalList = new IntervalList(filteredIntervals.getMetadata().getSequenceDictionary());
        filteredIntervals.getIntervals().forEach(i -> filteredIntervalList.add(new Interval(i)));
        filteredIntervalList.write(outputFilteredIntervalsFile);

        return "SUCCESS";
    }

    private void validateFilesAndResolveIntervals() {
        //validate input intervals and files
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);
        if (inputAnnotatedIntervalsFile == null && inputReadCountFiles.isEmpty()) {
            throw new UserException("Must provide annotated intervals or counts.");
        }
        if (inputAnnotatedIntervalsFile != null) {
            IOUtils.canReadFile(inputAnnotatedIntervalsFile);
        }
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");

        //parse inputs to resolve intervals and validate annotated intervals, if provided
        if (inputReadCountFiles.isEmpty()) {
            //only annotated intervals provided (no counts)
            annotatedIntervals = new AnnotatedIntervalCollection(inputAnnotatedIntervalsFile);
            specifiedIntervals = new SimpleIntervalCollection(
                    annotatedIntervals.getMetadata(),
                    intervalArgumentCollection.getIntervals(annotatedIntervals.getMetadata().getSequenceDictionary()));
            Utils.validateArg(specifiedIntervals.size() != 0, "At least one interval must be specified.");
            Utils.validateArg(new HashSet<>(annotatedIntervals.getIntervals()).containsAll(specifiedIntervals.getIntervals()),
                    "Annotated intervals do not contain all specified intervals.");
        } else {
            //counts provided
            final File firstReadCountFile = inputReadCountFiles.get(0);
            specifiedIntervals = CopyNumberArgumentValidationUtils.resolveIntervals(
                    firstReadCountFile, intervalArgumentCollection, logger);
            if (inputAnnotatedIntervalsFile != null) {
                //both annotated intervals and counts provided
                annotatedIntervals = CopyNumberArgumentValidationUtils.validateAnnotatedIntervalsSubset(
                        inputAnnotatedIntervalsFile, specifiedIntervals, logger);
            }
        }
    }

    private SimpleIntervalCollection filterIntervals() {
        final int numOriginalIntervals = specifiedIntervals.size();
        final boolean[] mask = new boolean[numOriginalIntervals];     //if true, filter out; each filter modifies this mask

        //apply annotation-based filters
        if (annotatedIntervals != null) {
            logger.info("Applying annotation-based filters...");
            //for present annotations, apply corresponding filters
            final List<AnnotationKey<?>> annotationKeys = annotatedIntervals.getRecords().get(0).getAnnotationMap().getKeys();
            if (annotationKeys.contains(CopyNumberAnnotations.GC_CONTENT)) {    //this should always be true, but we check it anyway
                updateMaskByAnnotationFilter(logger, annotatedIntervals, mask,
                        CopyNumberAnnotations.GC_CONTENT, "GC-content",
                        minimumGCContent, maximumGCContent);
            }
            if (annotationKeys.contains(CopyNumberAnnotations.MAPPABILITY)) {
                updateMaskByAnnotationFilter(logger, annotatedIntervals, mask,
                        CopyNumberAnnotations.MAPPABILITY, "mappability",
                        minimumMappability, maximumMappability);
            }
            if (annotationKeys.contains(CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT)) {
                updateMaskByAnnotationFilter(logger, annotatedIntervals, mask,
                        CopyNumberAnnotations.SEGMENTAL_DUPLICATION_CONTENT, "segmental-duplication-content",
                        minimumSegmentalDuplicationContent, maximumSegmentalDuplicationContent);
            }
        }

        //apply count-based filters
        if (!inputReadCountFiles.isEmpty()) {
            //get the read-count matrix (samples x specified intervals) and validate intervals and sequence dictionaries
            final RealMatrix readCountMatrix = constructReadCountMatrix(logger, inputReadCountFiles, specifiedIntervals);
            final int numSamples = readCountMatrix.getRowDimension();
            logger.info("Applying count-based filters...");

            //low-count filter: filter out intervals with a count strictly less than lowCountFilterCountThreshold
            //for strictly greater than lowCountFilterPercentageOfSamples
            IntStream.range(0, numOriginalIntervals)
                    .filter(i -> !mask[i])
                    .forEach(i -> {
                        if (Arrays.stream(readCountMatrix.getColumn(i))
                                .filter(c -> c < lowCountFilterCountThreshold)
                                .count() > lowCountFilterPercentageOfSamples * numSamples / 100.) {
                            mask[i] = true;
                        }
                    });
            logger.info(String.format("After applying low-count filter " +
                            "(intervals with a count < %d in > %s%% of samples fail), " +
                            "%d / %d intervals remain...",
                    lowCountFilterCountThreshold, lowCountFilterPercentageOfSamples,
                    countNumberPassing(mask), numOriginalIntervals));

            //extreme-count filter: filter out remaining intervals with counts that fall outside of the per-sample percentiles
            //[extremeCountMinimumPercentile, extremeCountMaximumPercentile] for strictly greater than extremeCountFilterPercentageOfSamples
            final boolean[][] percentileMask = new boolean[numSamples][numOriginalIntervals];
            for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                final double[] counts = readCountMatrix.getRow(sampleIndex);
                final double[] filteredCounts = IntStream.range(0, numOriginalIntervals)
                        .filter(i -> !mask[i])
                        .mapToDouble(i -> counts[i])
                        .toArray();
                final double extremeCountMinimumPercentileThreshold = extremeCountFilterMinimumPercentile == 0.
                        ? 0.
                        : new Percentile(extremeCountFilterMinimumPercentile).evaluate(filteredCounts);
                final double extremeCountMaximumPercentileThreshold = extremeCountFilterMaximumPercentile == 0.
                        ? 0.
                        : new Percentile(extremeCountFilterMaximumPercentile).evaluate(filteredCounts);
                for (int intervalIndex = 0; intervalIndex < numOriginalIntervals; intervalIndex++) {
                    final double count = readCountMatrix.getEntry(sampleIndex, intervalIndex);
                    if (!(extremeCountMinimumPercentileThreshold <= count && count <= extremeCountMaximumPercentileThreshold)) {
                        percentileMask[sampleIndex][intervalIndex] = true;
                    }
                }
            }
            IntStream.range(0, numOriginalIntervals)
                    .filter(i -> !mask[i])
                    .forEach(i -> {
                        if (IntStream.range(0, numSamples)
                                .filter(sampleIndex -> percentileMask[sampleIndex][i])
                                .count() > extremeCountFilterPercentageOfSamples * numSamples / 100.) {
                            mask[i] = true;
                        }
                    });
            logger.info(String.format("After applying extreme-count filter " +
                            "(intervals with a count percentile outside of [%s, %s] in > %s%% of samples fail), " +
                            "%d / %d intervals remain...",
                    extremeCountFilterMinimumPercentile, extremeCountFilterMaximumPercentile, extremeCountFilterPercentageOfSamples,
                    countNumberPassing(mask), numOriginalIntervals));
        }

        logger.info(String.format("%d / %d intervals passed all filters...", countNumberPassing(mask), numOriginalIntervals));

        //return the filtered intervals as a SimpleIntervalCollection
        return new SimpleIntervalCollection(
                specifiedIntervals.getMetadata(),
                IntStream.range(0, numOriginalIntervals)
                        .filter(i -> !mask[i])
                        .mapToObj(i -> specifiedIntervals.getRecords().get(i))
                        .collect(Collectors.toList()));
    }

    private static void updateMaskByAnnotationFilter(final Logger logger,
                                                     final AnnotatedIntervalCollection annotatedIntervals,
                                                     final boolean[] mask,
                                                     final AnnotationKey<Double> annotationKey,
                                                     final String filterName,
                                                     final double minValue,
                                                     final double maxValue) {
        IntStream.range(0, annotatedIntervals.size())
                .filter(i -> !mask[i])
                .forEach(i -> {
                    final double value = annotatedIntervals.getRecords().get(i).getAnnotationMap().getValue(annotationKey);
                    if (!(minValue <= value && value <= maxValue)) {
                        mask[i] = true;
                    }});
        logger.info(String.format("After applying %s filter (intervals with values outside of [%s, %s] fail), %d / %d intervals remain...",
                filterName, minValue, maxValue, countNumberPassing(mask), annotatedIntervals.size()));
    }

    private static RealMatrix constructReadCountMatrix(final Logger logger,
                                                       final List<File> inputReadCountFiles,
                                                       final SimpleIntervalCollection specifiedIntervals) {
        logger.info("Validating and aggregating input read-counts files...");
        final int numSamples = inputReadCountFiles.size();
        final int numIntervals = specifiedIntervals.size();
        //construct the interval subset to pull out from the read-count files
        final Set<SimpleInterval> intervalSubset = new HashSet<>(specifiedIntervals.getRecords());
        final RealMatrix readCountMatrix = new Array2DRowRealMatrix(numSamples, numIntervals);
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-counts file %s (%d / %d)", inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);
            if (!CopyNumberArgumentValidationUtils.isSameDictionary(
                    readCounts.getMetadata().getSequenceDictionary(),
                    specifiedIntervals.getMetadata().getSequenceDictionary())) {
                logger.warn(String.format("Sequence dictionary for read-counts file %s is inconsistent with those for other inputs.", inputReadCountFile));
            }
            final double[] subsetReadCounts = readCounts.getRecords().stream()
                    .filter(c -> intervalSubset.contains(c.getInterval()))
                    .mapToDouble(SimpleCount::getCount)
                    .toArray();
            Utils.validateArg(subsetReadCounts.length == intervalSubset.size(),
                    String.format("Intervals for read-count file %s do not contain all specified intervals.",
                            inputReadCountFile));
            readCountMatrix.setRow(sampleIndex, subsetReadCounts);
        }
        return readCountMatrix;
    }

    private static int countNumberPassing(final boolean[] mask) {
        final int numPassing = (int) IntStream.range(0, mask.length).filter(i -> !mask[i]).count();
        if (numPassing == 0) {
            throw new UserException.BadInput("Filtering removed all intervals.  Select less strict filtering criteria.");
        }
        return numPassing;
    }
}
