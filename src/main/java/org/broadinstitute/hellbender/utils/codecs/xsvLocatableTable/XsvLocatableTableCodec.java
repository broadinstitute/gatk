package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.InputStream;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.stream.Collectors;

/**
 * Codec class to read from XSV (e.g. csv, tsv, etc.) files.
 * Designed specifically with use by {@link org.broadinstitute.hellbender.tools.funcotator.Funcotator} in mind.
 *
 * Files that can be parsed by the {@link XsvLocatableTableCodec} will have a sibling configuration file of the same
 * name and the `.config` extension.  This file will contain the following keys:
 *      contig
 *      start
 *      end
 *      delimiter
 *      name
 *
 * These tables are assumed to have comment lines that start with `#` and a header that has the names for each
 * column in the table as the top row.
 *
 * Two or three columns will specify the location of each row in the data (contig, start, end; start and end can be the same
 * column).
 *
 * Created by jonn on 12/4/17.
 */
public final class XsvLocatableTableCodec extends AsciiFeatureCodec<XsvTableFeature> {

    private static final Logger logger = LogManager.getLogger(XsvLocatableTableCodec.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String COMMENT_DELIMITER = "#";

    public static final String CONFIG_FILE_CONTIG_COLUMN_KEY = "contig_column";
    public static final String CONFIG_FILE_START_COLUMN_KEY = "start_column";
    public static final String CONFIG_FILE_END_COLUMN_KEY = "end_column";
    public static final String CONFIG_FILE_DELIMITER_KEY = "xsv_delimiter";
    public static final String CONFIG_FILE_DATA_SOURCE_NAME_KEY = "name";

    //==================================================================================================================
    // Private Static Members:


    private static final String CONFIG_FILE_EXTENSION = ".config";

    //==================================================================================================================
    // Private Members:

    /** Column number from which to get the contig string for each entry. */
    private int contigColumn;

    /** Column number from which to get the start position for each entry. */
    private int startColumn;

    /** Column number from which to get the end position for each entry. */
    private int endColumn;

    /** Delimiter for entries in this XSV Table. */
    private String delimiter;

    /** The name of the data source that is associated with this {@link XsvLocatableTableCodec}. */
    private String dataSourceName;

    /** The XSV Table Header */
    private List<String> header;

    /** The current position in the file that is being read. */
    private long currentLine = 0;

    //==================================================================================================================
    // Constructors:

    public XsvLocatableTableCodec() {
        super(XsvTableFeature.class);
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public boolean canDecode(final String path) {

        // Get the paths to our file and the config file:
        final Path inputFilePath = IOUtils.getPath(path);
        final Path configFilePath = getConfigFilePath(inputFilePath);

        // Check that our files are good for eating... I mean reading...
        if ( validateInputDataFile(inputFilePath) && validateInputDataFile(configFilePath) ) {

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

        if (s.startsWith(COMMENT_DELIMITER)) {
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

        return new XsvTableFeature(contigColumn, startColumn, endColumn, header, split, dataSourceName);
    }

    @Override
    public List<String> readActualHeader(final LineIterator reader) {
        // All leading lines with comments / header info are headers:
        while ( reader.hasNext() ) {

            final String line = reader.next();
            ++currentLine;

            // Ignore commented out lines:
            if ( !line.startsWith(COMMENT_DELIMITER) ) {

                // The first non-commented line is the column header.
                // Add the data source name to teh start of each header row,
                // then add those rows to the header object.
                header = Arrays.stream(line.split(delimiter))
                        .map(x -> dataSourceName + "_" + x)
                        .collect(Collectors.toCollection(ArrayList::new));

                return header;
            }
        }

        throw new UserException.BadInput("Given file is malformed - does not contain a header!");
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
        contigColumn   = Integer.valueOf(configProperties.getProperty(CONFIG_FILE_CONTIG_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", ""));
        startColumn    = Integer.valueOf(configProperties.getProperty(CONFIG_FILE_START_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", ""));
        endColumn      = Integer.valueOf(configProperties.getProperty(CONFIG_FILE_END_COLUMN_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", ""));
        dataSourceName = configProperties.getProperty(CONFIG_FILE_DATA_SOURCE_NAME_KEY).replaceAll("^\\s+", "").replaceAll("\\s+$", "");

        // Get the delimiter - we do NOT remove whitespace here on purpose:
        delimiter      = configProperties.getProperty(CONFIG_FILE_DELIMITER_KEY);

        // Process delimiter just in case it is a tab escape character:
        if ( delimiter.equals("\\t") ) {
            delimiter = "\t";
        }
    }



    //==================================================================================================================
    // Helper Data Types:

}
