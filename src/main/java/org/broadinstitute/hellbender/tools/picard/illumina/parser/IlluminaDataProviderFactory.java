package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;
import static htsjdk.samtools.util.CollectionUtil.makeSet;
import static htsjdk.samtools.util.Log.getInstance;
import static htsjdk.samtools.util.StringUtil.join;
import static java.util.Arrays.asList;
import static java.util.Collections.unmodifiableSet;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaDataType.*;
import static org.broadinstitute.hellbender.tools.picard.illumina.parser.IlluminaFileUtil.SupportedIlluminaFormat.*;

/**
 * IlluminaDataProviderFactory accepts options for parsing Illumina data files for a lane and creates an
 * IlluminaDataProvider, an iterator over the ClusterData for that lane, which utilizes these options.
 * <p>
 * <p>
 * Note: Since we tend to use IlluminaDataProviderFactory in multithreaded environments (e.g. we call makeDataProvider
 * in a different thread per tile in IlluminaBasecallsToSam).  I've made it essentially immutable.  makeDataProvider/getTiles
 * are now idempotent (well as far as IlluminaDataProviderFactory is concerned, many file handles and other things are
 * opened when makeDataProvider is called).  We may in the future want dataTypes to be provided to the
 * makeDataProvider factory methods so configuration is not done multiple times for the same basecallDirectory in
 * client code.
 *
 * @author jburke@broadinstitute.org
 */
public class IlluminaDataProviderFactory {
    private static final Log log = getInstance(IlluminaDataProviderFactory.class);

    /**
     * A map of data types to a list of file formats in the order in which we prefer those file types (E.g. we would rather parse Bcls before QSeqs, Locs files before Clocs files ...)
     * We try to prefer data types that will be the fastest to parse/smallest in memory
     * NOTE: In the code below, if Qseq is chosen to provide for ANY data type then it is used for ALL its data types (since we'll have to parse the entire line for each Qseq anyways)
     */
    private static final Map<IlluminaDataType, List<SupportedIlluminaFormat>> DATA_TYPE_TO_PREFERRED_FORMATS = new HashMap<>();

    static {
        /** For types found in Qseq, we prefer the NON-Qseq file formats first.  However, if we end up using Qseqs then we use Qseqs for EVERY type it provides,
         * see determineFormats
         */
        DATA_TYPE_TO_PREFERRED_FORMATS.put(BaseCalls, makeList(
                MultiTileBcl, Bcl));
        DATA_TYPE_TO_PREFERRED_FORMATS.put(QualityScores, makeList(
                MultiTileBcl, Bcl));
        DATA_TYPE_TO_PREFERRED_FORMATS.put(PF, makeList(
                MultiTileFilter, Filter));
        DATA_TYPE_TO_PREFERRED_FORMATS.put(Position, makeList(
                MultiTileLocs, Locs, Clocs,
                Pos));

        DATA_TYPE_TO_PREFERRED_FORMATS.put(Barcodes, makeList(Barcode));
    }

    // The following properties must be specified by caller.
    /**
     * basecallDirectory holds QSeqs or bcls *
     */
    private final File basecallDirectory;
    private final int lane;

    /**
     * Whether or not to apply EAMSS filtering if parsing BCLs for the bases and quality scores.
     */
    private boolean applyEamssFiltering = true;

    /**
     * A Map of file formats to the dataTypes they will provide for this run.
     */
    protected final Map<SupportedIlluminaFormat, Set<IlluminaDataType>> formatToDataTypes;

    /**
     * Basecall Directory/lane parameterized util for finding IlluminaFiles
     */
    private final IlluminaFileUtil fileUtil;


    private final List<Integer> availableTiles;

    private final OutputMapping outputMapping;
    private final BclQualityEvaluationStrategy bclQualityEvaluationStrategy;

    /**
     * Create factory with the specified options, one that favors using QSeqs over all other files
     *
     * @param basecallDirectory The baseCalls directory of a complete Illumina directory.  Files are found by searching relative to this folder (some of them higher up in the directory tree).
     * @param lane              Which lane to iterate over.
     * @param readStructure     The read structure to which output clusters will conform.  When not using QSeqs, EAMSS masking(see BclParser) is run on individual reads as found in the readStructure, if
     *                          the readStructure specified does not match the readStructure implied by the sequencer's output than the quality scores output may differ than what would be found
     *                          in a run's QSeq files
     * @param dataTypesArg      Which data types to read
     */
    public IlluminaDataProviderFactory(final File basecallDirectory, final int lane, final ReadStructure readStructure,
                                       final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                                       final IlluminaDataType... dataTypesArg) {
        this(basecallDirectory, null,
                lane, readStructure,
                bclQualityEvaluationStrategy,
                dataTypesArg);
    }

    /**
     * Create factory with the specified options, one that favors using QSeqs over all other files
     *
     * @param basecallDirectory The baseCalls directory of a complete Illumina directory.  Files are found by searching relative to this folder (some of them higher up in the directory tree).
     * @param barcodesDirectory The barcodesDirectory with barcode files extracted by 'ExtractIlluminaBarcodes' (optional, use basecallDirectory if not specified)
     * @param lane              Which lane to iterate over.
     * @param readStructure     The read structure to which output clusters will conform.  When not using QSeqs, EAMSS masking(see BclParser) is run on individual reads as found in the readStructure, if
     *                          the readStructure specified does not match the readStructure implied by the sequencer's output than the quality scores output may differ than what would be found
     *                          in a run's QSeq files
     * @param dataTypesArg      Which data types to read
     */
    public IlluminaDataProviderFactory(final File basecallDirectory, File barcodesDirectory, final int lane,
                                       final ReadStructure readStructure,
                                       final BclQualityEvaluationStrategy bclQualityEvaluationStrategy, final IlluminaDataType... dataTypesArg) {
        this.basecallDirectory = basecallDirectory;
        this.bclQualityEvaluationStrategy = bclQualityEvaluationStrategy;

        this.lane = lane;
        /**
         * The types of data that will be returned by any IlluminaDataProviders created by this factory.
         * Note: In previous version, data of types not specified might be returned if a data type was specified
         * for data residing in QSeqs (since QSeqs span multiple data types).  This is no longer the case, you
         * MUST specify all data types that should be returned.
         */
        final Set<IlluminaDataType> dataTypes = unmodifiableSet(new HashSet<>(asList(dataTypesArg)));

        if (dataTypes.isEmpty()) {
            throw new IlluminaParserException("No data types have been specified for basecall output " + basecallDirectory + ", lane " + lane);
        }

        this.fileUtil = new IlluminaFileUtil(basecallDirectory, barcodesDirectory, lane);

        //find what request IlluminaDataTypes we have files for and select the most preferred file format available for that type
        formatToDataTypes = determineFormats(dataTypes, fileUtil);

        //find if we have any IlluminaDataType with NO available file formats and, if any exist, throw an exception
        final Set<IlluminaDataType> unmatchedDataTypes = findUnmatchedTypes(dataTypes, formatToDataTypes);
        if (unmatchedDataTypes.size() > 0) {
            throw new IlluminaParserException("Could not find a format with available files for the following data types: "
                    + join(", ", new ArrayList<>(unmatchedDataTypes)));
        }

        log.debug("The following file formats will be used by IlluminaDataProvider: " + join("," + formatToDataTypes.keySet()));

        availableTiles = fileUtil.getActualTiles(new ArrayList<>(formatToDataTypes.keySet()));
        if (availableTiles.isEmpty()) {
            throw new IlluminaParserException("No available tiles were found, make sure that " + basecallDirectory.getAbsolutePath() + " has a lane " + lane);
        }

        outputMapping = new OutputMapping(readStructure);
    }

    /**
     * Sometimes (in the case of skipped reads) the logical read structure of the output cluster data is different from the input
     * readStructure
     *
     * @return The ReadStructure describing the output cluster data
     */
    public ReadStructure getOutputReadStructure() {
        return outputMapping.getOutputReadStructure();
    }

    /**
     * Return the list of tiles available for this flowcell and lane.  These are in ascending numerical order.
     *
     * @return List of all tiles available for this flowcell and lane.
     */
    public List<Integer> getAvailableTiles() {
        return availableTiles;
    }

    /**
     * Sets whether or not EAMSS filtering will be applied if parsing BCL files for bases and quality scores.
     */
    public void setApplyEamssFiltering(final boolean applyEamssFiltering) {
        this.applyEamssFiltering = applyEamssFiltering;
    }

    /**
     * Call this method to create a ClusterData iterator over all clusters for all tiles in ascending numeric order.
     *
     * @return An iterator for reading the Illumina basecall output for the lane specified in the ctor.
     */
    public IlluminaDataProvider makeDataProvider() {
        return makeDataProvider(null);
    }

    /**
     * Call this method to create a ClusterData iterator over the specified tiles.
     *
     * @return An iterator for reading the Illumina basecall output for the lane specified in the constructor.
     */
    public IlluminaDataProvider makeDataProvider(List<Integer> requestedTiles) {
        if (requestedTiles == null) {
            requestedTiles = availableTiles;
        } else {
            if (requestedTiles.size() == 0) {
                throw new IlluminaParserException("Zero length tile list supplied to makeDataProvider, you must specify at least 1 tile OR pass NULL to use all available tiles");
            }
        }

        final Map<IlluminaParser, Set<IlluminaDataType>> parsersToDataType = new HashMap<>();
        for (final Map.Entry<SupportedIlluminaFormat, Set<IlluminaDataType>> fmToDt : formatToDataTypes.entrySet()) {
            parsersToDataType.put(makeParser(fmToDt.getKey(), requestedTiles), fmToDt.getValue());
        }

        log.debug("The following parsers will be used by IlluminaDataProvider: " + join("," + parsersToDataType.keySet()));

        return new IlluminaDataProvider(outputMapping, parsersToDataType, basecallDirectory, lane);
    }

    /**
     * Given a set of formats to data types they provide, find any requested data types that do not have a format associated with them and return them
     *
     * @param requestedDataTypes   Data types that need to be provided
     * @param formatToMatchedTypes A map of file formats to data types that will support them
     * @return The data types that go unsupported by the formats found in formatToMatchedTypes
     */
    public static Set<IlluminaDataType> findUnmatchedTypes(final Set<IlluminaDataType> requestedDataTypes, final Map<SupportedIlluminaFormat, Set<IlluminaDataType>> formatToMatchedTypes) {
        final Set<IlluminaDataType> copiedTypes = new HashSet<>(requestedDataTypes);
        for (final Set<IlluminaDataType> matchedTypes : formatToMatchedTypes.values()) {
            copiedTypes.removeAll(matchedTypes);
        }

        return copiedTypes;
    }

    /**
     * For all requestedDataTypes return a map of file format to set of provided data types that covers as many requestedDataTypes as possible and
     * chooses the most preferred available formats possible
     *
     * @param requestedDataTypes Data types to be provided
     * @param fileUtil           A file util for the lane/directory we wish to provide data for
     * @return A Map<Supported file format, Set of data types file format provides>
     */
    public static Map<SupportedIlluminaFormat, Set<IlluminaDataType>> determineFormats(final Set<IlluminaDataType> requestedDataTypes, final IlluminaFileUtil fileUtil) {
        //For predictable ordering and uniqueness only, put the requestedDataTypes into a treeSet
        final SortedSet<IlluminaDataType> toSupport = new TreeSet<>(requestedDataTypes);
        final Map<SupportedIlluminaFormat, Set<IlluminaDataType>> fileTypeToDataTypes = new HashMap<>();
        final Map<IlluminaDataType, SupportedIlluminaFormat> dataTypeToFormat = new HashMap<>();

        for (final IlluminaDataType ts : toSupport) {
            final SupportedIlluminaFormat preferredFormat = findPreferredAvailableFormat(ts, fileUtil);
            if (preferredFormat != null) {
                dataTypeToFormat.put(ts, preferredFormat);
            }
        }

        for (final IlluminaDataType dt : toSupport) {
            final SupportedIlluminaFormat format = dataTypeToFormat.get(dt);

            if (format != null) {
                if (fileTypeToDataTypes.containsKey(format)) {
                    fileTypeToDataTypes.get(format).add(dt);
                } else {
                    fileTypeToDataTypes.put(dataTypeToFormat.get(dt), makeSet(dt));
                }
            }
        }

        return fileTypeToDataTypes;
    }

    /**
     * Given a data type find the most preferred file format that also has files available
     *
     * @param dt       Type of desired data
     * @param fileUtil Util for the lane/directory in which we will find data
     * @return The file format that is "most preferred" (i.e. fastest to parse/smallest in memory)
     */
    private static SupportedIlluminaFormat findPreferredAvailableFormat(final IlluminaDataType dt, final IlluminaFileUtil fileUtil) {
        return findPreferredFormat(dt, fileUtil, true);
    }

    /**
     * Given a data type find the most preferred file format even if files are not available
     *
     * @param dt       Type of desired data
     * @param fileUtil Util for the lane/directory in which we will find data
     * @return The file format that is "most preferred" (i.e. fastest to parse/smallest in memory)
     */
    public static SupportedIlluminaFormat findPreferredFormat(final IlluminaDataType dt, final IlluminaFileUtil fileUtil) {
        return findPreferredFormat(dt, fileUtil, false);
    }

    private static SupportedIlluminaFormat findPreferredFormat(final IlluminaDataType dt, final IlluminaFileUtil fileUtil,
                                                               final boolean checkAvailable) {
        final List<SupportedIlluminaFormat> preferredFormats = DATA_TYPE_TO_PREFERRED_FORMATS.get(dt);
        SupportedIlluminaFormat format = null;
        for (int i = 0; i < preferredFormats.size() && format == null; i++) {
            if (checkAvailable && fileUtil.getUtil(preferredFormats.get(i)).filesAvailable()) {
                format = preferredFormats.get(i);
            } else if (!checkAvailable) {
                format = preferredFormats.get(i);
            }
        }

        return format;
    }

    /**
     * There are multiple parsers for the same IlluminaDataType (e.g. BCLParser and QSeqParser).  Instantiate an instance of the preferred parser for
     * the given data type with the information available and return it.
     *
     * @param format         The type of data we want to parse
     * @param requestedTiles The requestedTiles over which we will be parsing data
     * @return A parser that will parse dataType data over the given requestedTiles and cycles and output it in groupings of the sizes specified in outputLengths
     */
    private IlluminaParser makeParser(final SupportedIlluminaFormat format, final List<Integer> requestedTiles) {
        final IlluminaParser parser;
        switch (format) {
            case Barcode:
                parser = new BarcodeParser(((PerTileFileUtil) fileUtil.getUtil(Barcode)).getFiles(requestedTiles));
                break;

            case Bcl: {
                final CycleIlluminaFileMap bclFileMap = ((PerTilePerCycleFileUtil) fileUtil.getUtil(Bcl))
                        .getFiles(requestedTiles, outputMapping.getOutputCycles());
                bclFileMap.assertValid(requestedTiles, outputMapping.getOutputCycles());
                parser = new BclParser(basecallDirectory, lane, bclFileMap, outputMapping, this.applyEamssFiltering, bclQualityEvaluationStrategy);
                break;
            }

            case Filter:
                final IlluminaFileMap filterFileMap = ((PerTileFileUtil) fileUtil.getUtil(Filter)).getFiles(requestedTiles);
                parser = new FilterParser(filterFileMap);
                break;

            case Locs:
            case Clocs:
            case Pos:
                final PerTileFileUtil fu = (PerTileFileUtil) fileUtil.getUtil(format);
                parser = new PosParser(fu.getFiles(requestedTiles), format);
                break;

            case MultiTileFilter:
                parser = ((MultiTileFilterFileUtil) fileUtil.getUtil(MultiTileFilter)).makeParser(requestedTiles);
                break;

            case MultiTileLocs:
                parser = ((MultiTileLocsFileUtil) fileUtil.getUtil(MultiTileLocs)).makeParser(requestedTiles);
                break;

            case MultiTileBcl: {
                final MultiTileBclFileUtil util = (MultiTileBclFileUtil) fileUtil.getUtil(MultiTileBcl);
                final CycleIlluminaFileMap bclFileMap = util.getFiles(requestedTiles, outputMapping.getOutputCycles());
                bclFileMap.assertValid(requestedTiles, outputMapping.getOutputCycles());
                parser = new MultiTileBclParser(basecallDirectory, lane, bclFileMap, outputMapping,
                        this.applyEamssFiltering, bclQualityEvaluationStrategy, util.tileIndex);
                break;
            }

            default:
                throw new IlluminaParserException("Unrecognized data type(" + format + ") found by IlluminaDataProviderFactory!");
        }

        return parser;
    }
}
