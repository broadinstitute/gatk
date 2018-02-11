package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
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
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleBinningSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.utils.CachedOverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.utils.ReadOrientation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.*;

/**
 * TODO
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

    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

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
     * Sequence dictionary
     */
    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Metadata contained in the BAM file.
     */
    private BinningSampleLocatableMetadata metadata;

    /**
     * Overlap detector
     */
    private CachedOverlapDetector<SimpleInterval> overlapDetector;

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

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>(super.getDefaultReadFilters());
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.FIRST_OF_PAIR); // this will make sure we don't double count
        filters.add(ReadFilterLibrary.PROPERLY_PAIRED);
        filters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        // this will only keep reads in pairs that are properly oriented and mapped on same chromosome
        // and lie within a few standard deviations from the mean of fragment size distributions
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


    //TODO clean up this code
    @Override
    public void onTraversalStart() {
        // validate arguments
        if (readArguments.getReadFilesNames().size() != 1) {
            throw new UserException.BadInput("This tool only accepts a single bam/sam/cram as input");
        }
        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

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

        sequenceDictionary = getBestAvailableSequenceDictionary();
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(metadata.getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in BAM does not match the master sequence dictionary.");
        }

        logger.log(Level.INFO, "Reading intervals...");
        final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        overlapDetector = new CachedOverlapDetector<>(intervals);

        //verify again that intervals do not overlap
        Utils.validateArg(intervals.stream().noneMatch(i -> overlapDetector.getOverlapDetector().getOverlaps(i).size() > 1),
                "Input intervals may not be overlapping.");

        intervals.stream().forEach(interval -> readCountDataMap.put(
                interval, ReadCountDataFactory.getReadCountDataObject(readCountType, interval, covariateBinCollection)));

        // Next we start the traversal:
        logger.log(Level.INFO, "Collecting read counts ...");
    }

    private ReadCountCovariateBinCollection buildCovariateBinCollection() {
        final List<ReadCountCovariateBinningConfiguration> configurations = new ArrayList<>();
        if (includeGCBinning) {
            configurations.add(ReadCountCovariateBinningConfiguration.FRAGMENT_GC_CONTENT.setParameters(GC_MIN_BIN_VALUE, GC_MAX_BIN_VALUE, numGCBins));
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
        final Pair<Integer, Integer> startEndFragmentPositions =
                ReadOrientation.getReadOrientation(read).getReadToStartAndEndFragmentPositionsMapper().apply(read);
        final int windowLeadingBases = readLocation.getStart() - startEndFragmentPositions.getLeft();
        final int windowTrailingBases = startEndFragmentPositions.getRight() - readLocation.getEnd();
        referenceContext.setWindow(windowLeadingBases, windowTrailingBases);
        final Locatable fragmentCenter = ReadOrientation.getFragmentCenter(read);

        final SimpleInterval overlappingInterval = overlapDetector.getOverlap(fragmentCenter);
        //if fragment doesn't overlap any of the provided intervals, do nothing
        if (overlappingInterval == null) {
            return;
        }
        readCountDataMap.get(overlappingInterval).updateReadCount(read, referenceContext);
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
}
