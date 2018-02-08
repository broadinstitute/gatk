package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.*;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.BinnedReadCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.BinningSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleBinningSampleLocatableMetadata;
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
@ExperimentalFeature
public final class CollectBinnedReadCounts extends ReadWalker {

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
    protected File outputFile = null;

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
     * Reference to the logger.
     */
    private static final Logger logger = LogManager.getLogger(CollectBinnedReadCounts.class);

    /**
     * Sequence dictionary TODO this will be replaced with SampleLocatableMetadata
     */
    private SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

    /**
     * Metadata contained in the BAM file.
     */
    private BinningSampleLocatableMetadata metadata;

    /**
     * List of non-overlapping intervals
     */
    private CachedBinarySearchIntervalList<SimpleInterval> intervalList;

    /**
     * Map from targets to their corresponding read count data
     */
    private final Map<SimpleInterval, ReadCountData> readCountDataMap = new LinkedHashMap<>();

    /**
     * Collection of read count bins, it could be {@code null} if type of read count collection is not binned
     */
    private ReadCountCovariateBinCollection covariateBinCollection = null;

    /**
     * List of binning configurations specified by the user input
     */
    private List<ReadCountCovariateBinningConfiguration> covariateBinningConfigurations;

    /**
     * Reference data source used to compute reference related covariate metrics such as GC content
     */
    private ReferenceDataSource referenceDataSource;
    
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
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        // validate arguments
        if (readArguments.getReadFilesNames().size() != 1) {
            throw new UserException.BadInput("This tool only accepts a single bam/sam/cram as input");
        }
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

        referenceDataSource = ReferenceDataSource.of(referenceArguments.getReferencePath());

        // initialize fields
        if (readCountType == ReadCountType.BINNED) {
            if (!(includeGCBinning || includeFragmentLengthBinning)) {
                throw new UserException(String.format("Bad argument combination: type of the read count collection" +
                        " is %s, but no binning options are declared", ReadCountType.BINNED.getReadCountTypeName()));
            }
            covariateBinCollection = buildCovariateBinCollection();
        }
        metadata = new SimpleBinningSampleLocatableMetadata(MetadataUtils.readSampleName(getHeaderForReads()),
                getHeaderForReads().getSequenceDictionary(), readCountType, covariateBinningConfigurations);

        if (!CopyNumberArgumentValidationUtils.isSameDictionary(metadata.getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in BAM does not match the master sequence dictionary.");
        }

        logger.log(Level.INFO, "Reading intervals...");
        final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);

        //TODO replace this with the OverlapDetector
        intervalList = new CachedBinarySearchIntervalList<>(intervalArgumentCollection.getIntervals(sequenceDictionary));

        intervals.stream().forEach(interval -> readCountDataMap.put(
                interval, ReadCountDataFactory.getReadCountDataObject(readCountType, interval, covariateBinCollection)));

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

        final BinnedReadCountCollection binnedReadCountCollection = new BinnedReadCountCollection(metadata,
                new ArrayList<>(readCountDataMap.values()), columns);
        binnedReadCountCollection.write(outputFile);

        logger.log(Level.INFO, "Writing counts done.");

        return "SUCCESS";
    }





    /**
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
