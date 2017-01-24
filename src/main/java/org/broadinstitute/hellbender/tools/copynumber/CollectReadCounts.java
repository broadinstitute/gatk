package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.*;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Perform binned or simple count coverage collection on single-sample BAM file. " +
                "There are two distinct modes this tool can run in: if read count type is SIMPLE_COUNT then" +
                " the tool will output a single count per interval, i.e. number of reads intersecting each interval. " +
                "If read count type is set to BINNED, then tool requires at least one binning covariate (GC or FL)" +
                "specified to run. The output will be number of reads that fall into each covariate bin " +
                "(possibly a multidimensional bin if multiple covariates are specified).",
        oneLineSummary = "Perform binned or simple count coverage collection on single-sample BAM file",
        programGroup = CopyNumberProgramGroup.class
)
public final class CollectReadCounts extends ReadWalker {

    public static final String LINE_SEPARATOR = "\n";
    public static final double GC_MIN_BIN_VALUE = 0.;
    public static final double GC_MAX_BIN_VALUE = 1.0;

    public static final String READ_COUNT_TYPE_SHORT_NAME = "RCT";
    public static final String READ_COUNT_TYPE_LONG_NAME = "readCountType";

    public static final String INCLUDE_GC_BINS_SHORT_NAME = "GC";
    public static final String INCLUDE_GC_BINS_FULL_NAME = "includeGCBins";

    public static final String NUMBER_GC_BINS_SHORT_NAME = "NGC";
    public static final String NUMBER_GC_BINS_FULL_NAME = "numberGCBins";

    public static final String INCLUDE_FRAGMENT_LENGTH_BINS_SHORT_NAME = "FL";
    public static final String INCLUDE_FRAGMENT_LENGTH_BINS_FULL_NAME = "includeFragmentLengthBins";

    public static final String NUMBER_FRAGMENT_LENGTH_BINS_SHORT_NAME = "FLN";
    public static final String NUMBER_FRAGMENT_LENGTH_BINS_FULL_NAME = "numberFragmentLengthBins";

    public static final String FRAGMENT_LENGTH_MIN_BIN_VALUE_SHORT_NAME = "FLMIN";
    public static final String FRAGMENT_LENGTH_MIN_BIN_VALUE_FULL_NAME = "fragmentLengthMinValue";

    public static final String FRAGMENT_LENGTH_MAX_BIN_VALUE_SHORT_NAME  = "FLMAX";
    public static final String FRAGMENT_LENGTH_MAX_BIN_VALUE_FULL_NAME = "fragmentLengthMaxValue";

    @Argument(
            doc = "Output tab separated file with the counts",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected File output = null;

    @Argument(
            doc = "Type of read count collection. It can be SIMPLE_COUNT for simple" +
                    " integer count output, or BINNED for a binned read count output",
            shortName = READ_COUNT_TYPE_SHORT_NAME,
            fullName = READ_COUNT_TYPE_LONG_NAME,
            optional = false
    )
    protected ReadCountType readCountType;

    @Argument(
            doc = "Include GC content(GC) binning. Note that if read count type argument is set to BINNED, then " +
                    "at least one the the binning must be included. ",
            shortName = INCLUDE_GC_BINS_SHORT_NAME,
            fullName = INCLUDE_GC_BINS_FULL_NAME,
            optional = true
    )
    protected boolean includeGCBinning = false;

    @Argument(
            doc = "Number of GC bins. GC bins will be distributed uniformly.",
            shortName = NUMBER_GC_BINS_SHORT_NAME,
            fullName = NUMBER_GC_BINS_FULL_NAME,
            optional = true
    )
    protected int numGCBins = 10;

    @Argument(
            doc = "Include fragment length(FL) binning. Note that if read count type argument is set to BINNED, then " +
                    "at least one the the binning must be included. ",
            shortName = INCLUDE_FRAGMENT_LENGTH_BINS_SHORT_NAME,
            fullName = INCLUDE_FRAGMENT_LENGTH_BINS_FULL_NAME,
            optional = true
    )
    protected boolean includeFragmentLengthBinning = false;

    @Argument(
            doc = "Number of fragment length bins. FL bins will be distributed uniformly.",
            shortName = NUMBER_FRAGMENT_LENGTH_BINS_SHORT_NAME,
            fullName = NUMBER_FRAGMENT_LENGTH_BINS_FULL_NAME,
            optional = true
    )
    protected int numFragmentLengthBins = 5;

    @Argument(
            doc = "The start value of the fragment length binning",
            shortName = FRAGMENT_LENGTH_MIN_BIN_VALUE_SHORT_NAME,
            fullName = FRAGMENT_LENGTH_MIN_BIN_VALUE_FULL_NAME,
            optional = true
    )
    protected int fragmentLengthBinStartValue = 0;

    @Argument(
            doc = "The end value of the fragment length binning",
            shortName = FRAGMENT_LENGTH_MAX_BIN_VALUE_SHORT_NAME,
            fullName = FRAGMENT_LENGTH_MAX_BIN_VALUE_FULL_NAME,
            optional = true
    )
    protected int fragmentLengthBinEndValue = 500;

    /**
     * Writer to the main output file indicated by {@link #output}.
     */
    private PrintWriter outputWriter;

    /**
     * Reference to the logger.
     */
    private static final Logger logger = LogManager.getLogger(CalculateTargetCoverage.class);

    /**
     * Sequence dictionary
     */
    private SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

    /**
     * List of non-overlapping intervals
     */
    private CachedBinarySearchIntervalList<SimpleInterval> intervalList;

    /**
     * Map from targets to their corresponding read count data
     */
    private Map<SimpleInterval, ReadCountData> readCountDataMap;

    /**
     * Collection of read count bins, it could be {@code null} if type of read count collection is not binned
     */
    private ReadCountCovariateBinCollection covariateBinCollection = null;

    /**
     * List of binning configurations specified by the user input
     */
    private List<ReadCountCovariateBinningConfiguration> covariateBinningConfigurations;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>(super.getDefaultReadFilters());
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);

        return filters;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // validate arguments
        if (readArguments.getReadFilesNames().size() != 1) {
            throw new UserException.BadInput("This tool only accepts a single bam/sam/cram as input");
        }
        SampleCollection sampleCollection = new SampleCollection(getHeaderForReads());
        if (sampleCollection.sampleCount() > 1) {
            throw new UserException.BadInput("We do not support BAM files with multiple samples");
        }

        // initialize fields
        if (readCountType == ReadCountType.BINNED) {
            if (!(includeGCBinning || includeFragmentLengthBinning)) {
                throw new UserException(String.format("Bad argument combination: type of the read count collection" +
                        " is %s, but no binning options are declared", ReadCountType.BINNED.getReadCountTypeName()));
            }
            covariateBinCollection = buildCovariateBinCollection();
        }
        final String sampleName = sampleCollection.sampleIds().get(0);
        logger.log(Level.INFO, "Reading targets locations from intervals...");

        intervalList = new CachedBinarySearchIntervalList<>(intervalArgumentCollection.getIntervals(sequenceDictionary));

        readCountDataMap = buildReadCountDataMap();
        // Open output files and write headers:
        outputWriter = openOutputWriter(output, composeMatrixOutputHeader(getCommandLine(), sampleName));

        // Next we start the traversal:
        logger.log(Level.INFO, "Collecting read counts ...");
    }

    private ReadCountCovariateBinCollection buildCovariateBinCollection() {
        final List<ReadCountCovariateBinningConfiguration> configurations = new ArrayList<>();
        if (includeGCBinning) {
            configurations.add(ReadCountCovariateBinningConfiguration.GC_CONTENT.setParameters(GC_MIN_BIN_VALUE, GC_MAX_BIN_VALUE, numGCBins));
        }
        if (includeFragmentLengthBinning) {
            configurations.add(ReadCountCovariateBinningConfiguration.FRAGMENT_LENGTH.setParameters(fragmentLengthBinStartValue, fragmentLengthBinEndValue, numFragmentLengthBins));
        }
        this.covariateBinningConfigurations = configurations;
        return new ReadCountCovariateBinCollection(covariateBinningConfigurations);
    }

    /**
     * Populate the read count data map by instantiating a read count data container for each interval
     *
     * @return never {@code null}.
     */
    private Map<SimpleInterval, ReadCountData> buildReadCountDataMap() {
        final Map<SimpleInterval, ReadCountData> result = new HashMap<>();
        intervalList.getSortedIntervals().stream()
                .forEach(interval -> result.put(interval, ReadCountDataFactory.getReadCountDataObject(readCountType, interval, covariateBinCollection)));
        return result;
    }

    /**
     * Opens the output file for writing with a print-writer.
     *
     * @param output     the output file.
     * @param headerText to be printed immediately after opening the writer.
     * @return never {@code null}.
     * @throws UserException.CouldNotCreateOutputFile if there was some problem creating or overwriting {@code output}.
     */
    private PrintWriter openOutputWriter(final File output, final String headerText) {
        try {
            final PrintWriter result = new PrintWriter(output);
            result.println(headerText);
            result.flush();
            return result;
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, e);
        }
    }

    /**
     * Composes the main output header.
     *
     * @param commandLine      the tool command line.
     * @param sampleName        the name of the sample
     * @return never {@code null}.
     */
    private String composeMatrixOutputHeader(final String commandLine, final String sampleName) {
        final String formatString = String.join(LINE_SEPARATOR,
                "##fileFormat    = tsv",
                "##commandLine   = %s",
                "##title         = Read counts per target",
                ("##" + ReadCountFileHeaderKey.READ_COUNT_TYPE.getHeaderKeyName() + " = %s"),
                ("##" + ReadCountFileHeaderKey.SAMPLE_NAME.getHeaderKeyName() + "    = %s"));

        StringBuilder header = new StringBuilder(String.format(formatString, commandLine, readCountType.toString(), sampleName));

        //append the binning configuration descriptor if read count collection type is binned
        if (this.readCountType == ReadCountType.BINNED) {
            final String binningConfigurationDescriptor = covariateBinningConfigurations.stream().
                    map(config -> config.toString()).reduce("", (s1, s2) -> s1.concat(s2));
            final String binningConfigurationLine = ("##" + ReadCountFileHeaderKey.BINNING_CONFIGURATION.getHeaderKeyName()
                    + " = " + binningConfigurationDescriptor);
            header.append(LINE_SEPARATOR + binningConfigurationLine);
        }

        return header.toString();
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        final SimpleInterval readLocation = referenceContext.getInterval();
        intervalList.findIntersectionRange(readLocation).
                forEach(intervalIndex -> readCountDataMap.get(intervalList.getSortedIntervals().get(intervalIndex)).updateReadCount(read));
    }

    @Override
    public Object onTraversalSuccess() {
        logger.log(Level.INFO, "Collecting read counts done.");
        logger.log(Level.INFO, "Writing counts ...");

        List<String> readCountTableColumns = new ArrayList<>();

        readCountTableColumns.addAll(ReadCountTableColumn.COLUMNS.names());
        readCountTableColumns.addAll(ReadCountDataFactory.getColumnsOfReadCountType(readCountType, covariateBinCollection).names());
        final TableColumnCollection columns = new TableColumnCollection(readCountTableColumns);

        try {
            BinnedReadCountWriter writer = new BinnedReadCountWriter(outputWriter, columns);
            for(SimpleInterval interval: intervalList.getSortedIntervals()) {
                writer.writeRecord(readCountDataMap.get(interval));
            }
        } catch(IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(output, "Could not create output file");
        }

        logger.log(Level.INFO, "Writing counts done.");

        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (outputWriter != null) {
            outputWriter.close();
        }
    }

    /**
     * This is the stripped down version of {@link HashedListTargetCollection} designed to only support binary search
     * of {@link Locatable} objects
     *
     * @param <E> a locatable object
     */
    private class CachedBinarySearchIntervalList<E extends Locatable> {

        private List<E> sortedIntervals;

        private final Comparator<Locatable> intervalComparator = IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

        private int lastBinarySearchPosition = -1;

        CachedBinarySearchIntervalList(final List<E> unsortedIntervalList) {
            Utils.nonEmpty(unsortedIntervalList, "Interval list cannot be empty");
            Utils.containsNoNull(unsortedIntervalList, "intervals may not be null");
            sortedIntervals = unsortedIntervalList.stream().sorted(intervalComparator).collect(Collectors.toList());
            checkForIntervalOverlaps(sortedIntervals);
        }

        private void checkForIntervalOverlaps(final List<E> sortedIntervals) {
            final OptionalInt failureIndex = IntStream.range(1, sortedIntervals.size())
                    .filter(i -> IntervalUtils.overlaps(sortedIntervals.get(i-1),sortedIntervals.get(i)))
                    .findFirst();

            if (failureIndex.isPresent()) {
                final int index = failureIndex.getAsInt();
                throw new IllegalArgumentException(
                        String.format("input intervals contain at least two overlapping intervals: %s and %s",
                                sortedIntervals.get(index-1),sortedIntervals.get(index)));
            }
        }

        private IndexRange findIntersectionRange(final Locatable location) {
            Utils.nonNull(location, "the input location cannot be null");
            final int searchIndex = cachedBinarySearch(location);
            if (searchIndex < 0) {
                //TODO output proper empty interval with start=stop position that indicates where location would be inserted was an interval there
                return new IndexRange(0, 0);
            } else {
                final int firstOverlappingIndex = extendSearchIndexBackwards(location, searchIndex);
                final int lastOverlappingIndex = extendSearchIndexForward(location, searchIndex);
                return new IndexRange(firstOverlappingIndex, lastOverlappingIndex + 1);
            }
        }

        /**
         * Get list of sorted intervals
         *
         * @return sorted intervals
         */
        private List<E> getSortedIntervals() {
            return sortedIntervals;
        }

        private int cachedBinarySearch(final Locatable location) {
            if (lastBinarySearchPosition < 0) {
                return lastBinarySearchPosition = uncachedBinarySearch(location);
            } else {
                if (IntervalUtils.overlaps(sortedIntervals.get(lastBinarySearchPosition), location)) {
                    return lastBinarySearchPosition;
                } else {
                    final int candidate = lastBinarySearchPosition + 1;
                    if (candidate == sortedIntervals.size()) {
                        return lastBinarySearchPosition = uncachedBinarySearch(location);
                    } else if (IntervalUtils.overlaps(sortedIntervals.get(candidate), location)) {
                        return lastBinarySearchPosition = candidate;
                    } else {
                        return lastBinarySearchPosition = uncachedBinarySearch(location);
                    }
                }
            }
        }

        private int uncachedBinarySearch(final Locatable location) {
            final int searchResult = Collections.binarySearch(sortedIntervals, location, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
            if (searchResult >= 0) {
                return searchResult;
            } else {
                final int insertIndex = - (searchResult + 1);
                if (insertIndex < sortedIntervals.size() && IntervalUtils.overlaps(sortedIntervals.get(insertIndex),location)) {
                    return insertIndex;
                } if (insertIndex > 0 && IntervalUtils.overlaps(sortedIntervals.get(insertIndex - 1),location)) {
                    return insertIndex - 1;
                } else {
                    return searchResult;
                }
            }
        }

        /**
         * Looks for the last index in {@link #sortedIntervals} that has an overlap with the input {@code location}.
         * starting at {@code startIndex} and assuming that the element at that index has an
         * overlap with {@code location}.
         */
        private int extendSearchIndexForward(final Locatable location, final int startIndex) {
            final ListIterator<E> it = sortedIntervals.listIterator(startIndex + 1);
            while (it.hasNext()) {
                final E next = it.next();
                if (!IntervalUtils.overlaps(location, next)) {
                    return it.previousIndex() - 1;
                }
            }
            return it.previousIndex();
        }

        /**
         * Looks for the first index in {@link #sortedIntervals} that has an overlap with the input {@code location}
         * starting at {@code startIndex} and assuming that the element at that index has an overlap with {@code location}.
         */
        private int extendSearchIndexBackwards(final Locatable location, final int startIndex) {
            final ListIterator<E> it = sortedIntervals.listIterator(startIndex);
            while (it.hasPrevious()) {
                final E previous = it.previous();
                if (!IntervalUtils.overlaps(location, previous)) {
                    return it.nextIndex() + 1;
                }
            }
            return it.nextIndex();
        }
    }
}
