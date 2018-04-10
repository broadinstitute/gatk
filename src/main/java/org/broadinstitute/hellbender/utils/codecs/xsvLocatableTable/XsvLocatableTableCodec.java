package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Codec class to read from XSV (e.g. csv, tsv, etc.) files.
 *
 * Files that can be parsed by the {@link XsvLocatableTableCodec} will have a sibling configuration file of the same
 * name and the `.config` extension.  This file will contain the following keys:
 *      contig
 *      start
 *      end
 *      delimiter
 *      name
 *
 * Alternatively, an arbitrary configuration file can be specified in the {@link XsvLocatableTableCodec#XsvLocatableTableCodec(Path)}.
 *  However, this cannot be used within tribble, such as in a FeatureWalker.
 *
 * These tables are assumed to have preamble lines that start with `#` (comments) or '@' (SAM File Header) and a
 * header that has the names for each column in the table as the top row.  The type of preamble (comments vs.
 * SAM File header) is detected automatically.  Heterogeneous preambles (mix of '#' and '@') are not supported and will
 * cause an exception to be thrown.
 *
 * This class can render a SAMFileHeader from a preamble that is only comments ("#"), but the header will empty
 *   (minus the comment fields (@CO))
 *
 * Two or three columns will specify the location of each row in the data (contig, start, end; start and end can be the same
 * column).
 *
 * contig, start, and end can be specified as a comma-separated list.
 *
 * Created by jonn on 12/4/17.
 */
public final class XsvLocatableTableCodec extends AsciiFeatureCodec<XsvTableFeature> {

    private static final Logger logger = LogManager.getLogger(XsvLocatableTableCodec.class);

    //==================================================================================================================
    // Public Static Members:

    private static final String COMMENT_DELIMITER = "#";
    private static final String SAM_FILE_HEADER_LINE_START = "@";

    public static final String CONFIG_FILE_CONTIG_COLUMN_KEY = "contig_column";
    public static final String CONFIG_FILE_START_COLUMN_KEY = "start_column";
    public static final String CONFIG_FILE_END_COLUMN_KEY = "end_column";
    public static final String CONFIG_FILE_DELIMITER_KEY = "xsv_delimiter";
    public static final String CONFIG_FILE_DATA_SOURCE_NAME_KEY = "name";
    public static final String SAM_FILE_HEADER_START = "@HD\tVN:";


    //==================================================================================================================
    // Private Static Members:


    private static final String CONFIG_FILE_EXTENSION = ".config";

    //==================================================================================================================
    // Private Members:

    /** Column name (or index) from which to get the contig string for each entry.  As specified in the input.*/
    private String inputContigColumn;

    /** Column name (or index) from which to get the start position for each entry.  As specified in the input. */
    private String inputStartColumn;

    /** Column name (or index) from which to get the end position for each entry.  As specified in the input. */
    private String inputEndColumn;

    /** Column name from which to get the contig string for each entry. */
    private String finalContigColumn;

    /** Column name from which to get the start position for each entry. */
    private String finalStartColumn;

    /** Column name from which to get the end position for each entry. */
    private String finalEndColumn;

    /** Delimiter for entries in this XSV Table. */
    private String delimiter;

    /** The name of the data source that is associated with this {@link XsvLocatableTableCodec}. */
    private String dataSourceName;

    /** The XSV Table Header with prepends (datasource name) already applied */
    private List<String> header;

    /** A mapping of the fields in {@link XsvLocatableTableCodec#header} to the corresponding position. */
    private Map<String, Integer> headerToIndex;

    /** The locatable columns, once determined, in contig, start, end order. */
    private List<String> locatableColumns;

    /** The current position in the file that is being read. */
    private long currentLine = 0;

    /** Config file to use instead of a sibling config file.  Null if not using an override.*/
    private Path overrideConfigFile = null;

    /** Comments or SamFileHeader, if any.  Never {@code null}.  */
    private List<String> preamble = new ArrayList<>();

    /** Will hold starting string for preamble lines. */
    private String preambleLineStart;

    private boolean isHeaderInitialized = false;

    private boolean hasFirstRecordBeenRead = false;

    //==================================================================================================================
    // Constructors:

    public XsvLocatableTableCodec() {
        super(XsvTableFeature.class);
    }

    /** Constructor for when a configuration file is specified instead of using a sibling config file.
     *
     * This cannot be used with auto decoding.
     *
     * @param overrideConfigFile {@link Path} to the file to use as a configuration file for the given file.
     */
    public XsvLocatableTableCodec(final Path overrideConfigFile) {
        super(XsvTableFeature.class);
        this.overrideConfigFile = overrideConfigFile;
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public boolean canDecode(final String path) {

        // seg files are handled by a different codec.  This check has to be done, since seg files will return true in
        //  this codec and the AnnotatedIntervalCodec.
        return !path.endsWith(".seg") && canDecodeMinusExtensionChecks(path);
    }

    /**
     * Minus checking the file extension, can this class decode the given path.
     *
     * TODO: This method should be inside an abstract superclass.  {@link XsvLocatableTableCodec} and {@link org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCodec} should inherit.  See https://github.com/broadinstitute/gatk/issues/4580
     *
     * @param path File to check.  Never {@code null}
     * @return true if the file can be decoded.  False otherwise.
     */
    public boolean canDecodeMinusExtensionChecks(final String path) {
        Utils.nonNull(path);

        // Get the paths to our file and the config file:
        final Path inputFilePath = IOUtils.getPath(path);
        final Path configFilePath = (overrideConfigFile != null ?
                overrideConfigFile : getConfigFilePath(inputFilePath));

        // Check that our files are good for eating... I mean reading...
        if ( validateInputDataFile(inputFilePath) && validateInputDataFile(configFilePath) ) {

            // auto-determine the preamble format
            preambleLineStart = determinePreambleLineStart(inputFilePath);

            // Get our metadata and set up our internals so we can read from this file:
            readMetadataFromConfigFile(configFilePath);
            return true;
        }
        else {
            return false;
        }
    }

    @Override
    public XsvTableFeature decode(final String s) {

        // Increment our line counter:
        ++currentLine;

        if (s.startsWith(preambleLineStart)) {
            return null;
        }

        final List<String> split = new ArrayList<>(Arrays.asList(s.split(delimiter)));
        if (split.size() < 1) {
            throw new UserException.BadInput("XSV file has a line with no delimiter at line number: " + currentLine);
        }
        else if ( split.size() < header.size() ) {
            while (split.size() < header.size() ) {
                split.add("");
            }
        }
        else if ( split.size() > header.size() ) {
            logger.warn("WARNING: Line " + currentLine + " does not have the same number of fields as header (" + split.size() + " > " + header.size() + ")!  Truncating fields from end...");
            while (split.size() > header.size() ) {
                split.remove( split.size() - 1 );
            }
        }
        // HACK: For tribble to work properly, we need to know the position (in bytes) of the file pointer.  Currently, that cannot be ascertained.
        // HACK:  this code will just detect if the header is being read again and ignore it.
        else if (!hasFirstRecordBeenRead && ( IntStream.range(0, split.size()).allMatch(i -> header.get(i).equals(split.get(i))))) {
            return null;
        }
        hasFirstRecordBeenRead = true;
        return new XsvTableFeature(headerToIndex.get(finalContigColumn), headerToIndex.get(finalStartColumn),
                headerToIndex.get(finalEndColumn), header, split, dataSourceName);
    }

    /**
     * {@inheritDoc}
     * Read until we get to the header of this xsv
     *
     * Dev note:  We also determine the actual locatable columns here.
     *
     * @param reader iterator of the lines in the file.  Never {@code null}.
     * @return a list of strings that are the header columns.  Throws exception if no valid header line is found.
     */
    @Override
    public List<String> readActualHeader(final LineIterator reader) {

        Utils.nonNull(reader);

        // All leading lines with preamble / header info are headers:
        while ( reader.hasNext() ) {

            final String line = reader.next();
            ++currentLine;

            // Ignore preamble lines:
            if (!isPreambleLine(line)) {

                // The first non-commented line is the column header.
                // Add the data source name to the start of each header row,
                // then add those rows to the header object.
                header = Arrays.stream(line.split(delimiter))
                        .map(x -> determinePrefixForHeader() + x)
                        .collect(Collectors.toCollection(ArrayList::new));
                headerToIndex = IntStream.range(0, header.size()).boxed()
                        .collect(Collectors.toMap(i-> header.get(i), Function.identity()));
                isHeaderInitialized = true;

                finalContigColumn = determineFinalColumn(inputContigColumn);
                finalStartColumn = determineFinalColumn(inputStartColumn);
                finalEndColumn = determineFinalColumn(inputEndColumn);
                validateFinalColumns();

                locatableColumns = Arrays.asList(finalContigColumn, finalStartColumn, finalEndColumn);

                assertLocatableColumnsInHeaderToIndex(locatableColumns, headerToIndex);

                return header;

            } else {
                preamble.add(line.substring(preambleLineStart.length()));
            }
        }

        throw new UserException.BadInput("Given file is malformed - does not contain a header!");
    }

    private void validateFinalColumns() {
        if (finalContigColumn.equals(finalStartColumn) || finalContigColumn.equals(finalEndColumn)) {
            throw new UserException.BadInput("Contig column: " + finalContigColumn +
                    " is the same as start or end column.  Start: " + finalStartColumn +
                    "  End: " + finalEndColumn);
        }
    }


    @VisibleForTesting
    String determineFinalColumn(final String rawInputListOrIndex) {
        return StringUtils.isNumeric(rawInputListOrIndex) ? header.get(Integer.valueOf(rawInputListOrIndex))
                : determinePrefixForHeader() + determineColumnNameToUse(rawInputListOrIndex);
    }

    private List<String> getRawHeaders() {
        assertHeaderInitialized();
        return header.stream().map(h -> getHeaderWithoutPrefix(h)).collect(Collectors.toList());
    }

    private String determineColumnNameToUse(final String rawInputListAsString) {
        final List<String> rawHeaders = getRawHeaders();
        final String[] rawInputList = StringUtils.split(rawInputListAsString, ",");
        final Optional<String> candidateName = Stream.of(rawInputList).filter(c -> rawHeaders.contains(c)).findFirst();
        if (!candidateName.isPresent()) {
            throw new UserException.BadInput("Input did not contain any headers from the list: " + rawInputListAsString);
        } else {
            return candidateName.get();
        }
    }

    private String getHeaderWithoutPrefix(final String headerField) {
        if (StringUtils.isEmpty(determinePrefixForHeader())) {
            return headerField;
        } else {
            return StringUtils.replace(headerField, determinePrefixForHeader(), "", 1);
        }
    }

    private void assertLocatableColumnsInHeaderToIndex(final List<String> locatableColumns, final Map<String, Integer> headerToIndex) {
        final List<String> missingColumns =
                locatableColumns.stream().filter(c -> headerToIndex.get(c) == null)
                        .map(c -> getHeaderWithoutPrefix(c))
                        .collect(Collectors.toList());

        if (missingColumns.size() > 0) {
            final String missingColumnsString = StringUtil.join(", ", missingColumns);
            throw new UserException.BadInput("Error in input file: cannot find the locatable column(s): " + missingColumnsString + ", though these were specified in the parsing configuration.  Do those columns need to be added to the input file?  Do you have a heterogenous preamble (e.g. lines that start with both '#' and '@') before the headers?  Does each line of your preamble start with the correct string ('" + preambleLineStart + "')?");
        }
    }

    private void assertHeaderInitialized() {
        if (!isHeaderInitialized) {
            throw new GATKException.ShouldNeverReachHereException("Method that needs an initialized header was called before header was initialized.");
        }
    }

    private String determinePrefixForHeader() {
        return (StringUtils.isEmpty(dataSourceName) ? "" : dataSourceName + "_");
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Get the properties from the given {@code configFilePath}, validate that all required properties are present,
     * and return the property map.
     * @param configFilePath {@link Path} to the configuration file.
     * @return The {@link Properties} as contained in the given {@code configFilePath}.
     */
    public static Properties getAndValidateConfigFileContents(final Path configFilePath) {

        Utils.nonNull(configFilePath);

        // Read in the contents of the config file:
        final Properties configFileContents = new Properties();
        try ( final InputStream inputStream = Files.newInputStream(configFilePath, StandardOpenOption.READ) ) {
            configFileContents.load(inputStream);
        }
        catch (final Exception ex) {
            throw new UserException.BadInput("Unable to read from XSV config file: " + configFilePath.toUri().toString(), ex);
        }

        // Validate that it has the right keys:
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_CONTIG_COLUMN_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_START_COLUMN_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_END_COLUMN_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_DELIMITER_KEY, configFilePath);
        assertConfigPropertiesContainsKey(configFileContents, CONFIG_FILE_DATA_SOURCE_NAME_KEY, configFilePath);

        return configFileContents;
    }

    private boolean isPreambleLine(final String line) {
        return line.startsWith(preambleLineStart);
    }

    /**
     * Gets the path to the corresponding configuration file for the given {@code inputFilePath}.
     * The resulting path may or may not exist.
     * @param inputFilePath The data file {@link Path} from which to construct the path to the configuration file.
     * @return The {@link Path} for the configuration file associated with {@code inputFilePath}.
     */
    public static Path getConfigFilePath(final Path inputFilePath) {
        // Check for a sibling config file with the same name, .config as extension
        final String configFilePath = IOUtils.replaceExtension( inputFilePath.toString(), CONFIG_FILE_EXTENSION );
        return Paths.get(configFilePath);
    }

    /**
     * Ensures that the given {@link Properties} contain the given key.
     * @param configProperties The {@link Properties} in which to look for the given key.
     * @param key The value to find in the given {@link Properties}.
     * @param configFilePath The {@link Path} for the config file from which {@link Properties} were derived.  Used for printing output only.
     */
    private static void assertConfigPropertiesContainsKey(final Properties configProperties, final String key, final Path configFilePath) {
        if ( !configProperties.stringPropertyNames().contains(key) ) {
            throw new UserException.BadInput("Config file for datasource (" + configFilePath.toUri().toString() + ") does not contain required key: " + key);
        }
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * Asserts that the given {@code filePath} is a valid file from which to read.
     * @param filePath The {@link Path} to the data file to validate.
     * @return {@code true} if the given {@code filePath} is valid; {@code false} otherwise.
     */
    private boolean validateInputDataFile(final Path filePath) {
        return Files.exists(filePath) && Files.isReadable(filePath) && !Files.isDirectory(filePath);
    }

    /**
     * Reads the metadata required for parsing from the given {@code configFilePath}.
     * @param configFilePath {@link Path} to the configuration file from which to read in and setup metadata values.
     */
    private void readMetadataFromConfigFile(final Path configFilePath) {

        final Properties configProperties = getAndValidateConfigFileContents(configFilePath);

        // Get the properties and remove the leading/trailing whitespace if there is any:
        inputContigColumn = configProperties.getProperty(CONFIG_FILE_CONTIG_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");
        inputStartColumn  = configProperties.getProperty(CONFIG_FILE_START_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");
        inputEndColumn    = configProperties.getProperty(CONFIG_FILE_END_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");
        dataSourceName    = configProperties.getProperty(CONFIG_FILE_DATA_SOURCE_NAME_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");

        // Get the delimiter - we do NOT remove whitespace here on purpose:
        delimiter         = configProperties.getProperty(CONFIG_FILE_DELIMITER_KEY);

        // Process delimiter just in case it is a tab escape character:
        if ( delimiter.equals("\\t") ) {
            delimiter = "\t";
        }
    }

    /**
     * Creates a copy of the internal preamble upon each invocation.
     * {@link #readActualHeader(LineIterator)} must have been called or this method will throw an exception.
     *
     * @return an immutable list of all of the preamble lines in the xsv
     */
    private ImmutableList<String> getPreamble() {
        assertHeaderInitialized();
        return ImmutableList.copyOf(preamble);
    }

    /**
     * Get the header from this {@link XsvLocatableTableCodec} without the columns that contain location information.
     *
     * {@link #readActualHeader(LineIterator)} must have been called before this method.  Otherwise, exception will be thrown.
     *
     * Specifically the columns specified by the following fields are not included:
     *  {@link XsvLocatableTableCodec#inputContigColumn}
     *  {@link XsvLocatableTableCodec#inputStartColumn}
     *  {@link XsvLocatableTableCodec#inputEndColumn}
     * @return The header for this {@link XsvLocatableTableCodec} without location columns.
     */
    public List<String> getHeaderWithoutLocationColumns() {
        assertHeaderInitialized();
        return header.stream().filter(h -> !locatableColumns.contains(h))
                .collect(Collectors.toList());
    }

    /**
     * Throw an exception if the given column name cannot be used for one of the locatable columns.
     * @param columnName candidate column name for one of the locatable fields (contig, start, or end)
     */
    public static void validateLocatableColumnName(final String columnName) {
        Utils.validateArg(!StringUtils.isEmpty(columnName), "column header is blank.");
        Utils.validateArg(!NumberUtils.isCreatable(columnName), "column header cannot be a number: " + columnName);
    }

    private String determinePreambleLineStart(final Path path) {

        try (final InputStream stream = Files.newInputStream(path)){

            byte[] buff = new byte[SAM_FILE_HEADER_START.length()];
            stream.read(buff, 0, SAM_FILE_HEADER_START.length());
            final boolean eq = Arrays.equals(buff, SAM_FILE_HEADER_START.getBytes());

            if (eq) {
                return SAM_FILE_HEADER_LINE_START;
            } else {
                return COMMENT_DELIMITER;
            }
        } catch ( final IOException e ) {
            throw new UserException.CouldNotReadInputFile("Could not read file: " + path.toString(), e);
        }
    }

    /**
     * Create a SAM File Header from the preamble.
     *
     * Must be called after {@link #readActualHeader(LineIterator)} or an exception will be thrown
     * This method will return an empty SAMFileHeader if the preamble was empty.
     * @return Always returns a SAMFileHeader, even if it is empty or comments only.  Never {@code null}
     */
    public SAMFileHeader renderSamFileHeader() {
        assertHeaderInitialized();
        final ImmutableList<String> preamble = getPreamble();
        if (isPreambleSamFileHeader(preamble)) {
            final List<String> samHeaderAsString = preamble.stream().map(p -> preambleLineStart + p).collect(Collectors.toList());
            return createSamFileHeader(samHeaderAsString);
        } else {
            return createSamFileHeaderWithCommentsOnly(preamble);
        }
    }

    private boolean isPreambleSamFileHeader(final List<String> preamble) {
        if ((preamble == null) || (preamble.size() == 0)) {
            return false;
        }

        final String testPreambleFirstLine =  preambleLineStart + preamble.get(0);
        return testPreambleFirstLine.startsWith(XsvLocatableTableCodec.SAM_FILE_HEADER_START);
    }

    /**
     * @return copy of the sam file header created from the given list of strings.  {@code null} is not possible
     */
    @VisibleForTesting
    static SAMFileHeader createSamFileHeader(final List<String> samFileHeaderAsStrings) {

        final LineReader reader = BufferedLineReader.fromString(StringUtils.join(samFileHeaderAsStrings, "\n"));
        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        return codec.decode(reader, null);
    }

    @VisibleForTesting
    static SAMFileHeader createSamFileHeaderWithCommentsOnly(final List<String> comments) {
        final LineReader reader = BufferedLineReader.fromString("");
        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        final SAMFileHeader result = codec.decode(reader, null);

        // We grab all of the comments from the first decode, just in case the codec has some default comments.
        final List<String> finalComments = new ArrayList<>();
        finalComments.addAll(result.getComments());
        finalComments.addAll(comments);
        result.setComments(finalComments);

        return result;
    }

    /** Do not call before {@link #readActualHeader(LineIterator)} or exception will be thrown..
     *
     * @return the name of the contig column.  Never {@code null}.
     */
    public String getContigColumn() {
        assertHeaderInitialized();
        return finalContigColumn;
    }

    /** Do not call before {@link #readActualHeader(LineIterator)} or exception will be thrown..
     *
     * @return the name of the start column.  Never {@code null}.
     */
    public String getStartColumn() {
        assertHeaderInitialized();
        return finalStartColumn;
    }

    /** Do not call before {@link #readActualHeader(LineIterator)} or exception will be thrown..
     *
     * @return the name of the end column.  Never {@code null}.
     */
    public String getEndColumn() {
        assertHeaderInitialized();
        return finalEndColumn;
    }
}
