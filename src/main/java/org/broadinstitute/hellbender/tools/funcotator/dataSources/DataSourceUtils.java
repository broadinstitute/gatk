package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.TranscriptSelectionMode;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf.VcfFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.LocatableXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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

    @VisibleForTesting
    static final         String  MANIFEST_VERSION_LINE_START    = "Version:";
    private static final String  MANIFEST_SOURCE_LINE_START     = "Source:";
    private static final String  MANIFEST_ALT_SOURCE_LINE_START = "Alternate Source:";
    @VisibleForTesting
    static final Pattern VERSION_PATTERN                = Pattern.compile(MANIFEST_VERSION_LINE_START + "\\s+(\\d+)\\.(\\d+)\\.(\\d\\d\\d\\d)(\\d\\d)(\\d\\d)(.*)");
    private static final Pattern SOURCE_PATTERN                 = Pattern.compile(MANIFEST_SOURCE_LINE_START + "\\s+(ftp.*)");
    private static final Pattern ALT_SOURCE_PATTERN             = Pattern.compile(MANIFEST_ALT_SOURCE_LINE_START + "\\s+(gs.*)");

    // Track our minimum version number here:
    @VisibleForTesting
    static final int MIN_MAJOR_VERSION_NUMBER = 1;
    @VisibleForTesting
    static final int MIN_MINOR_VERSION_NUMBER = 4;
    @VisibleForTesting
    static final int MIN_YEAR_RELEASED        = 2018;
    @VisibleForTesting
    static final int MIN_MONTH_RELEASED       = 8;
    @VisibleForTesting
    static final int MIN_DAY_RELEASED         = 29;

    //==================================================================================================================
    // Public Static Members:

    /** The minimum version of the data sources required for funcotator to run.  */
    public static final String CURRENT_MINIMUM_DATA_SOURCE_VERSION         = String.format("v%d.%d.%d%02d%02d", MIN_MAJOR_VERSION_NUMBER, MIN_MINOR_VERSION_NUMBER, MIN_YEAR_RELEASED, MIN_MONTH_RELEASED, MIN_DAY_RELEASED);
    public static final String MANIFEST_FILE_NAME                          = "MANIFEST.txt";
    public static final String DATA_SOURCES_FTP_PATH                       = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/";
    public static final String DATA_SOURCES_BUCKET_PATH                    = "gs://broad-public-datasets/funcotator/";
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
    public static Map<Path, Properties> getAndValidateDataSourcesFromPaths( final String refVersion,
                                                                            final List<String> dataSourceDirectories) {
        Utils.nonNull(refVersion);
        Utils.nonNull(dataSourceDirectories);

        final Map<Path, Properties> metaData = new LinkedHashMap<>();

        boolean hasGencodeDataSource = false;

        // Go through our directories:
        final Set<String> names = new LinkedHashSet<>();
        for ( final String pathString : dataSourceDirectories ) {

            logger.info("Initializing data sources from directory: " + pathString);

            final Path pathToDatasources = IOUtils.getPath(pathString);
            if ( !isValidDirectory(pathToDatasources) ) {
                throw new UserException("ERROR: Given data source path is not a valid directory: " + pathToDatasources.toUri());
            }

            // Log information from the datasources directory so we can have a record of what we're using:
            final boolean isGoodVersionOfDataSources = logDataSourcesInfo(pathToDatasources);

            if ( !isGoodVersionOfDataSources ) {
                continue;
            }

            // Now that we have a valid directory, we need to grab a list of sub-directories in it:
            try {
                for ( final Path dataSourceTopDir : Files.list(pathToDatasources).filter(DataSourceUtils::isValidDirectory).collect(Collectors.toSet()) ) {

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
                throw new GATKException("Unable to read contents of: " + pathToDatasources.toUri().toString(), ex);
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
                case VCF:
                    funcotationFactory = DataSourceUtils.createVcfDataSource(path, properties, annotationOverridesMap, transcriptSelectionMode, userTranscriptIdSet);
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
        final String name      = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_NAME);

        // Create our gencode factory:
        return new GencodeFuncotationFactory(dataSourceFile.resolveSibling(fastaPath),
                        version,
                        name,
                        transcriptSelectionMode,
                        userTranscriptIdSet,
                        annotationOverridesMap
                );
    }

    /**
     * Create a {@link VcfFuncotationFactory} from filesystem resources and field overrides.
     * @param dataSourceFile {@link Path} to the data source file.  Must not be {@code null}.
     * @param dataSourceProperties {@link Properties} consisting of the contents of the config file for the data source.  Must not be {@code null}.
     * @param annotationOverridesMap {@link LinkedHashMap}{@code <String->String>} containing any annotation overrides to be included in the resulting data source.  Must not be {@code null}.
     * @param transcriptSelectionMode {@link TranscriptSelectionMode} to use when choosing the transcript for detailed reporting.  Must not be {@code null}.
     * @param userTranscriptIdSet {@link Set} of {@link String}s containing transcript IDs of interest to be selected for first.  Must not be {@code null}.
     * @return A new {@link GencodeFuncotationFactory} based on the given data source file information, field overrides map, and transcript information.
     */
    public static VcfFuncotationFactory createVcfDataSource(final Path dataSourceFile,
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
        final String name      = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_NAME);
        final String srcFile   = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_SRC_FILE);
        final String version   = dataSourceProperties.getProperty(CONFIG_FILE_FIELD_NAME_VERSION);

        // Create our VCF factory:

        return new VcfFuncotationFactory(
                name,
                version,
                dataSourceFile.resolveSibling(srcFile).toAbsolutePath(),
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

    /**
     * Logs the meta data about the data sources from the given {@code dataSourcesPath}.
     * We assume the data sources path is OK in the case that the version information cannot be read because the
     * user can create their own data sources directory, which may not contain the metadata we seek.
     *
     * NOTE: The MANIFEST file in a Data Sources directory is assumed to have the following properties:
     *       - Its name must be {@link #MANIFEST_FILE_NAME}
     *       - It must contain a line starting with {@link #MANIFEST_VERSION_LINE_START} containing an alphanumeric string containing the version number information.
     *           - This version information takes the form of:
     *               [MAJOR_VERSION].[MINOR_VERSION].[RELEASE_YEAR][RELEASE_MONTH][RELEASE_DAY][VERSION_DECORATOR]?
     *               e.g.
     *               1.1.20180204        (version 1.1 released Feb. 2, 2018)
     *               4.2.20480608somatic (version 4.2 released June 6, 2048 - somatic data sources)
     *               1.7.20190918X       (version 1.7 released Sept. 18, 2048 - X data sources)
     *
     *
     * @param dataSourcesPath {@link Path} to a Data Sources directory to check.
     * @return {@code True} if the given {@code dataSourcesPath} is equal to or newer than the minumum version, or if the version is unreadable.  {@code False} otherwise.
     */
    private static boolean logDataSourcesInfo(final Path dataSourcesPath) {

        boolean dataSourcesPathIsAcceptable = true;

        final Path manifestPath = dataSourcesPath.resolve(IOUtils.getPath(MANIFEST_FILE_NAME));

        String version = null;

        if ( Files.exists(manifestPath) && Files.isRegularFile(manifestPath) && Files.isReadable(manifestPath) ) {

            try ( final BufferedReader reader = Files.newBufferedReader(manifestPath) ) {

                Integer versionMajor     = null;
                Integer versionMinor     = null;
                Integer versionYear      = null;
                Integer versionMonth     = null;
                Integer versionDay       = null;
                String  versionDecorator = null;
                String  source           = null;
                String  alternateSource  = null;

                // Get the info from our README file:
                String line = reader.readLine();
                while ((line != null) && ((version == null) || (source == null) || (alternateSource == null))) {

                    if (version == null && line.startsWith(MANIFEST_VERSION_LINE_START)) {
                        final Matcher matcher = VERSION_PATTERN.matcher(line);
                        if ( matcher.matches() ) {
                            versionMajor     = Integer.valueOf(matcher.group(1));
                            versionMinor     = Integer.valueOf(matcher.group(2));
                            versionYear      = Integer.valueOf(matcher.group(3));
                            versionMonth     = Integer.valueOf(matcher.group(4));
                            versionDay       = Integer.valueOf(matcher.group(5));
                            versionDecorator = matcher.group(6);

                            version = versionMajor + "." + versionMinor + "." + versionYear + "" + versionMonth + "" + versionDay;
                        }
                        else {
                            logger.warn("README file has improperly formatted version string: " + line);
                        }
                    }

                    if (source == null && line.startsWith(MANIFEST_SOURCE_LINE_START)) {
                        final Matcher m = SOURCE_PATTERN.matcher(line);
                        if ( m.matches() ) {
                            source = m.group(1);
                        }
                        else {
                            logger.warn("README file has improperly formatted source string: " + line);
                        }
                    }

                    if (alternateSource == null && line.startsWith(MANIFEST_ALT_SOURCE_LINE_START)) {
                        final Matcher m = ALT_SOURCE_PATTERN.matcher(line);
                        if ( m.matches() ) {
                            alternateSource = m.group(1);
                        }
                        else {
                            logger.warn("README file has improperly formatted alternate source string: " + line);
                        }
                    }

                    line = reader.readLine();
                }

                // Make sure we have good info:
                if ( version == null ) {
                    logger.warn("Unable to read version information from data sources info/readme file: " + manifestPath.toUri().toString());
                }
                else {
                    logger.info("Data sources version: " + version);

                    // Make sure our version is OK:
                    dataSourcesPathIsAcceptable = validateVersionInformation(versionMajor, versionMinor, versionYear, versionMonth, versionDay);
                }

                if ( source == null ) {
                    logger.warn("Unable to read source information from data sources info/readme file: " + manifestPath.toUri().toString());
                }
                else {
                    logger.info("Data sources source: " + source);
                }

                if ( alternateSource == null ) {
                    logger.warn("Unable to read alternate source information from data sources info/readme file: " + manifestPath.toUri().toString());
                }
                else {
                    logger.info("Data sources alternate source: " + alternateSource);
                }
            }
            catch (final Exception ex) {
                logger.warn("Could not read " + MANIFEST_FILE_NAME + ": unable to log data sources version information.", ex);
            }
        }
        else {
            logger.warn("Could not read " + MANIFEST_FILE_NAME + ": unable to log data sources version information.");
        }

        // Warn the user if they need newer stuff.
        if ( !dataSourcesPathIsAcceptable ) {

            String message = "";
            message = message + "ERROR: Given data source path is too old!  Minimum required version is: " + CURRENT_MINIMUM_DATA_SOURCE_VERSION + " (yours: " + version + ")\n";
            message = message + "       You must download a newer version of the data sources from the Broad Institute FTP site: " + DATA_SOURCES_FTP_PATH + "\n";
            message = message + "       or the Broad Institute Google Bucket: " + DATA_SOURCES_BUCKET_PATH + "\n";
            throw new UserException( message );
        }

        return dataSourcesPathIsAcceptable;
    }

    @VisibleForTesting
    static boolean validateVersionInformation(final int major, final int minor, final int year, final int month, final int day) {

        // Compare from largest to smallest differences:

        if ( major < MIN_MAJOR_VERSION_NUMBER ) {
            return false;
        }

        if ( minor <  MIN_MINOR_VERSION_NUMBER ) {
            return false;
        }

        if ( year < MIN_YEAR_RELEASED ) {
            return false;
        }

        else if ( year == MIN_YEAR_RELEASED ) {
            if ( month < MIN_MONTH_RELEASED ) {
                return false;
            }
            else if ( month == MIN_MONTH_RELEASED ) {
                if ( day < MIN_DAY_RELEASED ) {
                    return false;
                }
            }
        }

        return true;
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

    /**
     * Sort by putting gencode datasources first and then non-gencode in alphabetical order.
     *
     * Note that this must match how the datasources are being rendered in the {@link Funcotator#apply(VariantContext, ReadsContext, ReferenceContext, FeatureContext)}
     * @param f1 datasource to compare.  Never {@code null}
     * @param f2 datasource to compare.  Never {@code null}
     * @return -1 if f1 should be put first, 1 if f2 should be put first.
     */
    public static int datasourceComparator(final DataSourceFuncotationFactory f1, final DataSourceFuncotationFactory f2) {
        Utils.nonNull(f1);
        Utils.nonNull(f2);
        final boolean isF1Gencode = f1.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE);
        final boolean isF2Gencode = f2.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE);
        if (isF1Gencode == isF2Gencode) {
            return f1.getInfoString().compareTo(f2.getInfoString());
        } else {
            if (isF1Gencode) {
                return -1;
            } else {
                return 1;
            }
        }
    }
}
