package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.TranscriptSelectionMode;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.LocatableXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Utilities for reading / working with / manipulating Data Sources.
 * Designed to be a static utility class with no state.
 * Created by jonn on 3/8/18.
 */
final public class DataSourceUtils {

    // Private constructor.  No makey, no brakey.
    private DataSourceUtils() {}

    //==================================================================================================================
    // Private Static Members:

    private static final Logger logger = LogManager.getLogger(DataSourceUtils.class);

    private static final PathMatcher configFileMatcher =
            FileSystems.getDefault().getPathMatcher("glob:**/*.config");

    //==================================================================================================================
    // Public Static Members:

    public static final String CONFIG_FILE_FIELD_NAME_NAME                 = "name";
    public static final String CONFIG_FILE_FIELD_NAME_VERSION              = "version";
    public static final String CONFIG_FILE_FIELD_NAME_SRC_FILE             = "src_file";
    public static final String CONFIG_FILE_FIELD_NAME_ORIGIN_LOCATION      = "origin_location";
    public static final String CONFIG_FILE_FIELD_NAME_PREPROCESSING_SCRIPT = "preprocessing_script";
    public static final String CONFIG_FILE_FIELD_NAME_TYPE                 = "type";
    public static final String CONFIG_FILE_FIELD_NAME_GENCODE_FASTA_PATH   = "gencode_fasta_path";
    public static final String CONFIG_FILE_FIELD_NAME_XSV_KEY              = "xsv_key";
    public static final String CONFIG_FILE_FIELD_NAME_XSV_KEY_COLUMN       = "xsv_key_column";
    public static final String CONFIG_FILE_FIELD_NAME_XSV_DELIMITER        = "xsv_delimiter";
    public static final String CONFIG_FILE_FIELD_NAME_XSV_PERMISSIVE_COLS  = "xsv_permissive_cols";
    public static final String CONFIG_FILE_FIELD_NAME_CONTIG_COLUMN        = "contig_column";
    public static final String CONFIG_FILE_FIELD_NAME_START_COLUMN         = "start_column";
    public static final String CONFIG_FILE_FIELD_NAME_END_COLUMN           = "end_column";

    //==================================================================================================================
    // Public Static Methods:

    /**
     * Initializes the data sources for {@link Funcotator}.
     * @param refVersion The version of the reference we're using to create annotations.  Must not be {@code null}.
     * @param dataSourceDirectories A {@link List} of {@link Path} to the directories containing our data sources.  Must not be {@code null}.
     * @return The contents of the config files for each of the data sources found in the given {@code dataSourceDirectories}.
     */
    public static Map<Path, Properties> getAndValidateDataSourcesFromPaths(final String refVersion,
                                                                            final List<String> dataSourceDirectories) {
        Utils.nonNull(refVersion);
        Utils.nonNull(dataSourceDirectories);

        final Map<Path, Properties> metaData = new LinkedHashMap<>();

        boolean hasGencodeDataSource = false;

        // Go through our directories:
        final Set<String> names = new LinkedHashSet<>();
        for ( final String pathString : dataSourceDirectories ) {
            final Path p = IOUtils.getPath(pathString);
            if ( !isValidDirectory(p) ) {
                logger.warn("WARNING: Given path is not a valid directory: " + p.toUri().toString());
                continue;
            }

            // Now that we have a valid directory, we need to grab a list of sub-directories in it:
            try {
                for ( final Path dataSourceTopDir : Files.list(p).filter(DataSourceUtils::isValidDirectory).collect(Collectors.toSet()) ) {

                    // Get the path that corresponds to our reference version:
                    final Path dataSourceDir = dataSourceTopDir.resolve(refVersion);

                    // Make sure that we have a good data source directory:
                    if ( isValidDirectory(dataSourceDir) ) {

                        // Get the config file path:
                        final Path configFilePath = getConfigfile(dataSourceDir);

                        // Read the config file into Properties:
                        final Properties configFileProperties = readConfigFileProperties(configFilePath);

                        // Validate the config file properties:
                        assertConfigFilePropertiesAreValid(configFileProperties, configFilePath);

                        // Make sure we only have 1 of this data source:
                        if ( names.contains(configFileProperties.getProperty("name")) ) {
                            throw new UserException.BadInput("Error: contains more than one dataset of name: "
                                    + configFileProperties.getProperty("name") + " - one is: " + configFilePath.toUri().toString());
                        }
                        else {
                            names.add( configFileProperties.getProperty("name") );

                            if ( FuncotatorArgumentDefinitions.DataSourceType.getEnum(configFileProperties.getProperty("type")) ==
                                    FuncotatorArgumentDefinitions.DataSourceType.GENCODE ) {
                                hasGencodeDataSource = true;
                            }

                            // Add our config file to the properties list:
                            metaData.put( configFilePath, configFileProperties );
                        }
                    }
                }
            }
            catch (final IOException ex) {
                throw new GATKException("Unable to read contents of: " + p.toUri().toString(), ex);
            }
        }

        // Sanity checks to make sure we actually found our data sources:
        if ( metaData.size() == 0 ) {
            throw new UserException("ERROR: Could not find any data sources for given reference: " + refVersion);
        }
        else if ( !hasGencodeDataSource ) {
            throw new UserException("ERROR: a Gencode datasource is required!");
        }

        return metaData;
    }

    /**
     * Gets the {@link Path} config file in the given directory.
     * Will throw a {@link UserException}if no config file is found.
     * @param directory The {@link Path} of the directory in which to search for a config file.  Must not be {@code null}.
     * @return The {@link Path} to the config file found in the given {@code directory}.
     */
    public static Path getConfigfile(final Path directory) {

        Utils.nonNull(directory);

        final List<Path> configFileSet;

        try {
            configFileSet = Files.list(directory)
                    .filter(x -> configFileMatcher.matches(x))
                    .filter(Files::exists)
                    .filter(Files::isReadable)
                    .filter(Files::isRegularFile)
                    .collect(Collectors.toList());

            if (configFileSet.size() > 1) {
                throw new UserException("ERROR: Directory contains more than one config file: " + directory.toUri().toString());
            }
            else if ( configFileSet.size() == 0 ) {
                throw new UserException("ERROR: Directory does not contain a config file: " + directory.toUri().toString());
            }

            return configFileSet.get(0);
        }
        catch (final IOException ex) {
            throw new UserException("Unable to read contents of: " + directory.toUri().toString(), ex);
        }
    }

    /** @return {@code true} if the given {@link Path} exists, is readable, and is a directory; {@code false} otherwise. */
    public static boolean isValidDirectory(final Path p) {
        Utils.nonNull(p);
        return Files.exists(p) && Files.isReadable(p) && Files.isDirectory(p);
    }

    /**
     * Create a {@link List} of {@link DataSourceFuncotationFactory} based on meta data on the data sources, overrides, and transcript reporting priority information.
     * @param dataSourceMetaData {@link Map} of {@link Path}->{@link Properties} containing metadata about each data source.  Must not be {@code null}.
     * @param annotationOverridesMap {@link LinkedHashMap} of {@link String}->{@link String} containing any annotation overrides to include in data sources.  Must not be {@code null}.
     * @param transcriptSelectionMode {@link TranscriptSelectionMode} to use when choosing the transcript for detailed reporting.  Must not be {@code null}.
     * @param userTranscriptIdSet {@link Set} of {@link String}s containing transcript IDs of interest to be selected for first.  Must not be {@code null}.
     * @return A {@link List} of {@link DataSourceFuncotationFactory} given the data source metadata, overrides, and transcript reporting priority information.
     */
    public static List<DataSourceFuncotationFactory> createDataSourceFuncotationFactoriesForDataSources(final Map<Path, Properties> dataSourceMetaData,
                                                                                               final LinkedHashMap<String, String> annotationOverridesMap,
                                                                                               final TranscriptSelectionMode transcriptSelectionMode,
                                                                                               final Set<String> userTranscriptIdSet) {

        Utils.nonNull(dataSourceMetaData);
        Utils.nonNull(annotationOverridesMap);
        Utils.nonNull(transcriptSelectionMode);
        Utils.nonNull(userTranscriptIdSet);

        final List<DataSourceFuncotationFactory> dataSourceFactories = new ArrayList<>(dataSourceMetaData.size());

        // Now we know we have unique and valid data.
        // Now we must instantiate our data sources:
        for ( final Map.Entry<Path, Properties> entry : dataSourceMetaData.entrySet() ) {

            logger.debug("Creating Funcotation Factory for " + entry.getValue().getProperty("name") + " ...");

            final Path path = entry.getKey();
            final Properties properties = entry.getValue();

            final DataSourceFuncotationFactory funcotationFactory;

            // Note: we need no default case since we know these are valid:
            final String stringType = properties.getProperty("type");
            switch ( FuncotatorArgumentDefinitions.DataSourceType.getEnum(stringType) ) {
                case LOCATABLE_XSV:
                    funcotationFactory = DataSourceUtils.createLocatableXsvDataSource(path, properties, annotationOverridesMap);
                    break;
                case SIMPLE_XSV:
                    funcotationFactory = DataSourceUtils.createSimpleXsvDataSource(path, properties, annotationOverridesMap);
                    break;
                case COSMIC:
                    funcotationFactory = DataSourceUtils.createCosmicDataSource(path, properties, annotationOverridesMap);
                    break;
                case GENCODE:
                    funcotationFactory = DataSourceUtils.createGencodeDataSource(path, properties, annotationOverridesMap, transcriptSelectionMode, userTranscriptIdSet);
                    break;
                default:
                    throw new GATKException("Unknown type of DataSourceFuncotationFactory encountered: " + stringType );
            }

            // Add in our factory:
            dataSourceFactories.add(funcotationFactory);
        }

        logger.debug("All Data Sources have been created.");
        return dataSourceFactories;
    }

    /**
     * Create a {@link LocatableXsvFuncotationFactory} from filesystem resources and field overrides.
     * @param dataSourceFile {@link Path} to the data source file.  Must not be {@code null}.
     * @param dataSourceProperties {@link Properties} consisting of the contents of the config file for the data source.  Must not be {@code null}.
     * @param annotationOverridesMap {@link LinkedHashMap}{@code <String->String>} containing any annotation overrides to be included in the resulting data source.  Must not be {@code null}.
     * @return A new {@link LocatableXsvFuncotationFactory} based on the given data source file information and field overrides map.
     */
    public static LocatableXsvFuncotationFactory createLocatableXsvDataSource(final Path dataSourceFile,
                                                                              final Properties dataSourceProperties,
                                                                              final LinkedHashMap<String, String> annotationOverridesMap) {
        Utils.nonNull(dataSourceFile);
        Utils.nonNull(dataSourceProperties);
        Utils.nonNull(annotationOverridesMap);

        final String name      = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_NAME);
        final String version   = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_VERSION);

        // Create a locatable XSV feature reader to handle XSV Locatable features:
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(name, version, annotationOverridesMap);

        // Set the supported fields by the LocatableXsvFuncotationFactory:
        locatableXsvFuncotationFactory.setSupportedFuncotationFields(
                new ArrayList<>(
                        Collections.singletonList(
                                dataSourceFile.resolveSibling(
                                        IOUtils.getPath( dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_SRC_FILE) )
                                )
                        )
                )
        );

        return locatableXsvFuncotationFactory;
    }

    /**
     * Create a {@link SimpleKeyXsvFuncotationFactory} from filesystem resources and field overrides.
     * @param dataSourceFile {@link Path} to the data source file.  Must not be {@code null}.
     * @param dataSourceProperties {@link Properties} consisting of the contents of the config file for the data source.  Must not be {@code null}.
     * @param annotationOverridesMap {@link LinkedHashMap}{@code <String->String>} containing any annotation overrides to be included in the resulting data source.  Must not be {@code null}.
     * @return A new {@link SimpleKeyXsvFuncotationFactory} based on the given data source file information and field overrides map.
     */
    public static SimpleKeyXsvFuncotationFactory createSimpleXsvDataSource(final Path dataSourceFile,
                                                                   final Properties dataSourceProperties,
                                                                   final LinkedHashMap<String, String> annotationOverridesMap) {

        Utils.nonNull(dataSourceFile);
        Utils.nonNull(dataSourceProperties);
        Utils.nonNull(annotationOverridesMap);

        // Create our SimpleKeyXsvFuncotationFactory:
        return new SimpleKeyXsvFuncotationFactory(
                        dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_NAME),
                        dataSourceFile.resolveSibling(IOUtils.getPath(dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_SRC_FILE))),
                        dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_VERSION),
                        dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_XSV_DELIMITER),
                        Integer.valueOf(dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_XSV_KEY_COLUMN)),
                        SimpleKeyXsvFuncotationFactory.XsvDataKeyType.valueOf(dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_XSV_KEY)),
                        annotationOverridesMap,
                        0,
                        Boolean.valueOf(dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_XSV_PERMISSIVE_COLS))
                );
    }

    /**
     * Create a {@link CosmicFuncotationFactory} from filesystem resources and field overrides.
     * @param dataSourceFile {@link Path} to the data source file.  Must not be {@code null}.
     * @param dataSourceProperties {@link Properties} consisting of the contents of the config file for the data source.  Must not be {@code null}.
     * @param annotationOverridesMap {@link LinkedHashMap}{@code <String->String>} containing any annotation overrides to be included in the resulting data source.  Must not be {@code null}.
     * @return A new {@link CosmicFuncotationFactory} based on the given data source file information and field overrides map.
     */
    public static CosmicFuncotationFactory createCosmicDataSource(final Path dataSourceFile,
                                                                final Properties dataSourceProperties,
                                                                final LinkedHashMap<String, String> annotationOverridesMap) {
        Utils.nonNull(dataSourceFile);
        Utils.nonNull(dataSourceProperties);
        Utils.nonNull(annotationOverridesMap);

        final String version   = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_VERSION);

        return new CosmicFuncotationFactory(
                        dataSourceFile.resolveSibling(IOUtils.getPath(dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_SRC_FILE))),
                        annotationOverridesMap,
                        version
                );
    }

    /**
     * Create a {@link GencodeFuncotationFactory} from filesystem resources and field overrides.
     * @param dataSourceFile {@link Path} to the data source file.  Must not be {@code null}.
     * @param dataSourceProperties {@link Properties} consisting of the contents of the config file for the data source.  Must not be {@code null}.
     * @param annotationOverridesMap {@link LinkedHashMap}{@code <String->String>} containing any annotation overrides to be included in the resulting data source.  Must not be {@code null}.
     * @param transcriptSelectionMode {@link TranscriptSelectionMode} to use when choosing the transcript for detailed reporting.  Must not be {@code null}.
     * @param userTranscriptIdSet {@link Set} of {@link String}s containing transcript IDs of interest to be selected for first.  Must not be {@code null}.
     * @return A new {@link GencodeFuncotationFactory} based on the given data source file information, field overrides map, and transcript information.
     */
    public static GencodeFuncotationFactory createGencodeDataSource(final Path dataSourceFile,
                                                                 final Properties dataSourceProperties,
                                                                 final LinkedHashMap<String, String> annotationOverridesMap,
                                                                 final TranscriptSelectionMode transcriptSelectionMode,
                                                                 final Set<String> userTranscriptIdSet) {

        Utils.nonNull(dataSourceFile);
        Utils.nonNull(dataSourceProperties);
        Utils.nonNull(annotationOverridesMap);
        Utils.nonNull(transcriptSelectionMode);
        Utils.nonNull(userTranscriptIdSet);

        // Get some metadata:
        final String fastaPath = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_GENCODE_FASTA_PATH);
        final String version   = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_VERSION);

        // Create our gencode factory:
        return new GencodeFuncotationFactory(dataSourceFile.resolveSibling(fastaPath),
                        version,
                        transcriptSelectionMode,
                        userTranscriptIdSet,
                        annotationOverridesMap
                );
    }

    //==================================================================================================================
    // Private Static Methods:

    private static Properties readConfigFileProperties(final Path configFilePath) {
        final Properties configProperties = new Properties();
        try ( final InputStream inputStream = Files.newInputStream(configFilePath, StandardOpenOption.READ) ) {
            configProperties.load(inputStream);
        }
        catch (final Exception ex) {
            throw new UserException.BadInput("Unable to read from data source config file: " + configFilePath.toUri().toString(), ex);
        }

        return configProperties;
    }

    // ========================================================================================================
    // Config file static helper methods:
    // ------------------------------------------------

    private static void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {

        // Universally required config properties:
        assertConfigPropertiesContainsKey(CONFIG_FILE_FIELD_NAME_NAME, configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey(CONFIG_FILE_FIELD_NAME_VERSION, configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey(CONFIG_FILE_FIELD_NAME_SRC_FILE, configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey(CONFIG_FILE_FIELD_NAME_ORIGIN_LOCATION, configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey(CONFIG_FILE_FIELD_NAME_PREPROCESSING_SCRIPT, configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey(CONFIG_FILE_FIELD_NAME_TYPE, configFileProperties, configFilePath);

        // Validate our source file:
        assertPathFilePropertiesField( configFileProperties, CONFIG_FILE_FIELD_NAME_SRC_FILE, configFilePath);

        // Validate our type:
        final String stringType = configFileProperties.getProperty(CONFIG_FILE_FIELD_NAME_TYPE);
        final FuncotatorArgumentDefinitions.DataSourceType type;
        try {
            type = FuncotatorArgumentDefinitions.DataSourceType.getEnum(stringType);
        }
        catch (final IllegalArgumentException ex) {
            throw new UserException.BadInput("ERROR in config file: " + configFilePath.toUri().toString() +
                    " - Invalid value in \"" + CONFIG_FILE_FIELD_NAME_TYPE + "\" field: " + stringType, ex);
        }

        // Validate remaining values based on type:
        type.assertConfigFilePropertiesAreValid(configFileProperties, configFilePath);
    }

    /**
     * Asserts that the given {@code key} is contained in the given {@code configProperties}.
     * @param key {@link String} name of the key, the existence of which will be confirmed in {@code configProperties}.
     * @param configProperties {@link Properties} corresponding to the given {@code configFilePath} in which to check for the existence of {@code key}.
     * @param configFilePath {@link Path} to config file.  For output purposes only.
     */
    public static void assertConfigPropertiesContainsKey(final String key, final Properties configProperties, final Path configFilePath) {
        if ( !configProperties.stringPropertyNames().contains(key) ) {
            throw new UserException.BadInput("Config file for datasource (" + configFilePath.toUri().toString() + ") does not contain required key: \"" + key + "\"");
        }
    }

    /**
     * Asserts that the given {@code field} is contained in the given {@code props} and is an integer value.
     * @param props {@link Properties} corresponding to the given {@code filePath} in which to check for the validity of {@code field}.
     * @param field {@link String} name of the field, the existence and correct type of which will be confirmed in {@code props}.
     * @param filePath {@link Path} to config file.  For output purposes only.
     */
    public static void assertIntegerPropertiesField(final Properties props, final String field, final Path filePath) {
        try {
            Integer.valueOf(props.getProperty(field));
        }
        catch (final NumberFormatException ex) {
            throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                    " - Invalid value in \"" + field + "\" field: " + props.getProperty(field));
        }
    }

    /**
     * Asserts that the given {@code field} is contained in the given {@code props} and is a boolean value.
     * @param props {@link Properties} corresponding to the given {@code filePath} in which to check for the validity of {@code field}.
     * @param field {@link String} name of the field, the existence and correct type of which will be confirmed in {@code props}.
     * @param filePath {@link Path} to config file.  For output purposes only.
     */
    public static void assertBooleanPropertiesField(final Properties props, final String field, final Path filePath) {

            final String comparableValue = props.getProperty(field).trim().toLowerCase();

            if ( !comparableValue.equals("true") && !comparableValue.equals("false") ) {
                throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                        " - Invalid value in \"" + field + "\" field: " + props.getProperty(field));
            }
    }

    /**
     * Asserts that the given {@code field} is contained in the given {@code props} and is a file path.
     * @param props {@link Properties} corresponding to the given {@code filePath} in which to check for the validity of {@code field}.
     * @param field {@link String} name of the field, the existence and correct type of which will be confirmed in {@code props}.
     * @param filePath {@link Path} to config file.  For output purposes only.
     */
    public static void assertPathFilePropertiesField(final Properties props, final String field, final Path filePath) {
        final Path sourceFilePath = filePath.resolveSibling(props.getProperty(field));
        if ( !Files.exists(sourceFilePath) ) {
            throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                    " - " + field + " does not exist: " + sourceFilePath);
        }
        else if ( !Files.isRegularFile(sourceFilePath) ) {
            throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                    " -  " + field + " is not a regular file: " + sourceFilePath);
        }
        else if ( !Files.isReadable(sourceFilePath) ) {
            throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                    " - " + field + " is not readable: " + sourceFilePath);
        }
    }
}
