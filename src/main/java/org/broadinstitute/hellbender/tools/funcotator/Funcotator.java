package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.LocatableXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Funcotator (FUNCtional annOTATOR) performs functional analysis on given variants
 * and reports output in a specified output file.
 *
 * This tool is the GATK analog of the Oncotator.
 *
 * Created by jonn on 8/22/17.
 */
@CommandLineProgramProperties(
        summary = "Create functional annotations on given variants cross-referenced by a given database.\n" +
                "A GATK version of the Oncotator.",
        oneLineSummary = "(Experimental) Functional Annotator",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class Funcotator extends VariantWalker {
    private static final Logger logger = LogManager.getLogger(Funcotator.class);

    private static final PathMatcher configFileMatcher =
            FileSystems.getDefault().getPathMatcher("glob:**/*.config");

    //==================================================================================================================
    // Arguments:

    //-----------------------------------------------------
    // Required args:

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF File.")
    protected File outputFile;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME,
            doc = "The version of the Human Genome reference to use (hg19 or Hhg38)."
    )
    protected FuncotatorArgumentDefinitions.ReferenceVersionType referenceVersion;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME,
            doc = "The paths to any datasource directories for Funcotator."
    )
    protected List<String> dataSourceDirectories;

    //-----------------------------------------------------
    // Optional args:

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME,
            optional = true,
            doc = "Method of detailed transcript selection."
    )
    protected FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME,
            optional = true,
            doc = "List of transcripts to use for annotation to override selected transcript."
    )
    protected Set<String> transcriptList = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_DEFAULT_VALUE;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME,
            optional = true,
            doc = "Default values for annotations that are not added (in the format <ANNOTATION>:<VALUE>)"
    )
    protected List<String> annotationDefaults = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_DEFAULT_VALUE;

    @Argument(
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME,
            optional = true,
            doc = "Override values for annotations.  Replaces existing matchihng annotations with given values (in the format <ANNOTATION>:<VALUE>)"
    )
    protected List<String> annotationOverrides = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_DEFAULT_VALUE;

    //==================================================================================================================

    private OutputRenderer outputRenderer;
    private final List<DataSourceFuncotationFactory> dataSourceFactories = new ArrayList<>();
    private List<GencodeFuncotationFactory> gencodeFuncotationFactories = new ArrayList<>();

    private List<FeatureInput<? extends Feature>> manualFeatureInputs = new ArrayList<>();

    //==================================================================================================================

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final LinkedHashMap<String, String> annotationDefaultsMap = splitAnnotationArgsIntoMap(annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = splitAnnotationArgsIntoMap(annotationOverrides);

        // Initialize all of our data sources:
        final Map<Path, Properties> configData = getAndValidateDataSourcesFromPaths(referenceVersion, dataSourceDirectories);
        initializeDataSources(configData, annotationOverridesMap);

        // Set up our other data source factories:

        // Determine which annotations are accounted for (by the funcotation factories) and which are not.
        final LinkedHashMap<String, String> unaccountedForDefaultAnnotations = getUnaccountedForAnnotations( dataSourceFactories, annotationDefaultsMap );
        final LinkedHashMap<String, String> unaccountedForOverrideAnnotations = getUnaccountedForAnnotations( dataSourceFactories, annotationOverridesMap );

        // Set up our output renderer:
        // TODO: in the future this should be encapsulated into a factory for output renderers based on an input argument.
        outputRenderer = new VcfOutputRenderer(getHeaderForVariants(),
                                               createVCFWriter(outputFile),
                                               dataSourceFactories,
                                               unaccountedForDefaultAnnotations,
                                               unaccountedForOverrideAnnotations);

        outputRenderer.open();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        if ( !referenceContext.hasBackingDataSource() ) {
            throw new GATKException("No reference context for variant.  Cannot annotate!");
        }

        // Place the variant on our queue to be funcotated:
        enqueueAndHandleVariant(variant, referenceContext, featureContext);
    }

    @Override
    public Object onTraversalSuccess() {
        return true;
    }

    @Override
    public void closeTool() {

        for(final DataSourceFuncotationFactory factory : dataSourceFactories) {
            factory.close();
        }
        outputRenderer.close();

    }

    //==================================================================================================================

    /**
     * Creates a {@link LinkedHashMap} of annotations in the given {@code annotationMap} that do not occur in the given {@code dataSourceFactories}.
     * @param dataSourceFactories {@link List} of {@link DataSourceFuncotationFactory} to check for whether each annotation in the {@code annotationMap} is handled.
     * @param annotationMap {@link Map} (of ANNOTATION_NAME : ANNOTATION_VALUE) to check
     * @return A {@link LinkedHashMap} of annotations in the given {@code annotationMap} that do not occur in the given {@code dataSourceFactories}.
     */
    private LinkedHashMap<String, String> getUnaccountedForAnnotations( final List<DataSourceFuncotationFactory> dataSourceFactories,
                                                                        final Map<String, String> annotationMap ) {
        final LinkedHashMap<String, String> outAnnotations = new LinkedHashMap<>();

        // Check each field in each factory:
        for ( final String field : annotationMap.keySet() ) {
            boolean accountedFor = false;
            for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

                if ( funcotationFactory.getSupportedFuncotationFields().contains(field) ) {
                    accountedFor = true;
                    break;
                }
            }
            if ( !accountedFor ) {
                outAnnotations.put(field, annotationMap.get(field));
            }
        }

        return outAnnotations;
    }

    /**
     * Creates an annotation on the given {@code variant} or enqueues it to be processed during a later call to this method.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variant}.
     */
    private void enqueueAndHandleVariant(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // Get our feature inputs:
        final List<Feature> featureList = new ArrayList<>();
        for ( final FeatureInput<? extends Feature> featureInput : manualFeatureInputs ) {
            featureList.addAll( featureContext.getValues(featureInput) );
        }

        // Create a place to keep our funcotations:
        final List<Funcotation> funcotations = new ArrayList<>();

        // Annotate with Gencode first:

        // Create a list of GencodeFuncotation to use for other Data Sources:
        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>();

        for ( final GencodeFuncotationFactory factory : gencodeFuncotationFactories ) {
            final List<Funcotation> funcotationsFromGencodeFactory = factory.createFuncotations(variant, referenceContext, featureList);
            funcotations.addAll(funcotationsFromGencodeFactory);
            gencodeFuncotations.addAll(
                    funcotationsFromGencodeFactory.stream()
                    .map(x -> (GencodeFuncotation)x)
                    .collect(Collectors.toList())
            );
        }

        // Annotate with the rest of the data sources:
        for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

            // Make sure we don't double up on the Gencodes:
            if ( funcotationFactory.getType().equals(FuncotatorArgumentDefinitions.DataSourceType.GENCODE) ) {
                continue;
            }

            funcotations.addAll( funcotationFactory.createFuncotations(variant, referenceContext, featureList, gencodeFuncotations) );
        }
        outputRenderer.write(variant, funcotations);
    }

    /**
     * Split each element of the given {@link List} into a key and value.
     * Assumes each element of the given {@link List} is formatted as follows:
     *     KEY:VALUE
     * @param annotationArgs {@link List} of strings formatted KEY:VALUE to turn into a {@link Map}.
     * @return A {@link LinkedHashMap} of KEY:VALUE pairs corresponding to entries in the given list.
     */
    private LinkedHashMap<String, String> splitAnnotationArgsIntoMap( final List<String> annotationArgs ) {

        final LinkedHashMap<String, String> annotationMap = new LinkedHashMap<>();

        for ( final String s : annotationArgs ) {
            final List<String> keyVal = ParsingUtils.split(s, FuncotatorArgumentDefinitions.MAP_NAME_VALUE_DELIMITER);
            if ( keyVal.size() != 2) {
                throw new UserException.BadInput( "Argument annotation incorrectly formatted: " + s );
            }

            annotationMap.put( keyVal.get(0), keyVal.get(1) );
        }

        return annotationMap;
    }

    /** @return {@code true} if the given {@link Path} exists, is readable, and is a directory; {@code false} otherwise. */
    private static boolean isValidDirectory(final Path p) {
        return Files.exists(p) && Files.isReadable(p) && Files.isDirectory(p);
    }

    /**
     * Gets the {@link Path} config file in the given directory.
     * Will throw a {@link UserException}if no config file is found.
     * @param directory The {@link Path} of the directory in which to search for a config file.
     * @return The {@link Path} to the config file found in the given {@code directory}.
     */
    private static Path getConfigfile(final Path directory) {
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
            throw new GATKException("Unable to read contents of: " + directory.toUri().toString());
        }
    }

    /**
     * Initializes the data sources for {@link Funcotator}.
     * @param refVersion The version of the reference we're using to create annotations.
     * @param dataSourceDirectories A {@link List} of {@link Path} to the directories containing our data sources.
     * @return The contents of the config files for each of the data sources found in the given {@code dataSourceDirectories}.
     */
    private Map<Path, Properties> getAndValidateDataSourcesFromPaths(final FuncotatorArgumentDefinitions.ReferenceVersionType refVersion,
                                                                     final List<String> dataSourceDirectories) {
        final Map<Path, Properties> metaData = new LinkedHashMap<>();

        boolean hasGencodeDataSource = false;

        // Go through our directories:
        final Set<String> names = new HashSet<>();
        for ( final String pathString : dataSourceDirectories ) {
            final Path p = IOUtils.getPath(pathString);
            if ( !isValidDirectory(p) ) {
                logger.warn("WARNING: Given path is not a valid directory: " + p.toUri().toString());
                continue;
            }

            // Now that we have a valid directory, we need to grab a list of sub-directories in it:
            try {
                for ( final Path dataSourceTopDir : Files.list(p).filter(Funcotator::isValidDirectory).collect(Collectors.toSet()) ) {

                    // Get the path that corresponds to our reference version:
                    final Path dataSourceDir = dataSourceTopDir.resolve(refVersion.toString());

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
                throw new GATKException("Unable to read contents of: " + p.toUri().toString());
            }
        }

        if ( !hasGencodeDataSource ) {
            throw new UserException("ERROR: a Gencode datasource is required!");
        }

        return metaData;
    }

    private void initializeDataSources(final Map<Path, Properties> metaData,
                                       final LinkedHashMap<String, String> annotationOverridesMap) {

        // Now we know we have unique and valid data.
        // Now we must instantiate our data sources:
        for ( final Map.Entry<Path, Properties> entry : metaData.entrySet() ) {

            logger.debug("Initializing " + entry.getValue().getProperty("name") + " ...");

            // Note: we need no default case since we know these are valid:
            final String stringType = entry.getValue().getProperty("type");
            switch ( FuncotatorArgumentDefinitions.DataSourceType.getEnum(stringType) ) {
                case LOCATABLE_XSV: createLocatableXsvDataSource(entry.getKey(), entry.getValue(), annotationOverridesMap); break;
                case SIMPLE_XSV:    createSimpleXsvDataSource(entry.getKey(), entry.getValue(), annotationOverridesMap); break;
                case COSMIC:        createCosmicDataSource(entry.getKey(), entry.getValue(), annotationOverridesMap); break;
                case GENCODE:       createGencodeDataSource(entry.getKey(), entry.getValue(), annotationOverridesMap); break;
            }
        }

        logger.debug("All Data Sources are Initialized.");
    }

    private void createLocatableXsvDataSource(final Path dataSourceFile,
                                              final Properties dataSourceProperties,
                                              final LinkedHashMap<String, String> annotationOverridesMap) {
        final String name   = dataSourceProperties.getProperty("name");

        // Inject our features into our list of feature data sources:
        final FeatureInput<? extends Feature> featureInput = addFeatureInputsAfterInitialization(
                dataSourceFile.resolveSibling(
                    IOUtils.getPath( dataSourceProperties.getProperty("src_file") )
                ).toUri().toString(),
                name,
                XsvTableFeature.class
        );

        // Add our feature input to our list of manual inputs:
        manualFeatureInputs.add(featureInput);

        // Create a locatable XSV feature reader to handle XSV Locatable features:
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory();

        // Set the supported fields by the LocatableXsvFuncotationFactory:
        locatableXsvFuncotationFactory.setSupportedFuncotationFields(
                new ArrayList<>(
                        Collections.singletonList(
                                dataSourceFile.resolveSibling(
                                    IOUtils.getPath( dataSourceProperties.getProperty("src_file") )
                                )
                        )
                )
        );
        dataSourceFactories.add( locatableXsvFuncotationFactory );
    }

    private void createSimpleXsvDataSource(final Path dataSourceFile,
                                           final Properties dataSourceProperties,
                                           final LinkedHashMap<String, String> annotationOverridesMap) {
        // Create our SimpleKeyXsvFuncotationFactory:
        final SimpleKeyXsvFuncotationFactory factory =
                //final String name, final Path filePath, final String delim, final int keyColumn, final XsvDataKeyType keyType
                new SimpleKeyXsvFuncotationFactory(
                    dataSourceProperties.getProperty("name"),
                    dataSourceFile.resolveSibling(IOUtils.getPath(dataSourceProperties.getProperty("src_file"))),
                    dataSourceProperties.getProperty("version"),
                    dataSourceProperties.getProperty("xsv_delimiter"),
                    Integer.valueOf(dataSourceProperties.getProperty("xsv_key_column")),
                    SimpleKeyXsvFuncotationFactory.XsvDataKeyType.valueOf(dataSourceProperties.getProperty("xsv_key")),
                    annotationOverridesMap,
                    0,
                    Boolean.valueOf(dataSourceProperties.getProperty("xsv_permissive_cols"))
                );

        // Add it to our sources:
        dataSourceFactories.add( factory );

    }

    private void createCosmicDataSource(final Path dataSourceFile,
                                        final Properties dataSourceProperties,
                                        final LinkedHashMap<String, String> annotationOverridesMap) {

        final CosmicFuncotationFactory cosmicFuncotationFactory =
                new CosmicFuncotationFactory(
                  dataSourceFile.resolveSibling(IOUtils.getPath(dataSourceProperties.getProperty("src_file"))),
                  annotationOverridesMap
                );

        // Add our factory to our factory list:
        dataSourceFactories.add( cosmicFuncotationFactory );
    }

    private void createGencodeDataSource(final Path dataSourceFile,
                                         final Properties dataSourceProperties,
                                         final LinkedHashMap<String, String> annotationOverridesMap) {

        // Get some metadata:
        final String fastaPath = dataSourceProperties.getProperty("gencode_fasta_path");
        final String version   = dataSourceProperties.getProperty("version");
        final String name   = dataSourceProperties.getProperty("name");

        // Inject our Gencode GTF features into our list of feature data sources:
        final FeatureInput<? extends Feature> featureInput = addFeatureInputsAfterInitialization(
                dataSourceFile.resolveSibling(
                        IOUtils.getPath( dataSourceProperties.getProperty("src_file") )
                ).toUri().toString(),
                name,
                GencodeGtfFeature.class
        );

        // Add our feature input to our list of manual inputs:
        manualFeatureInputs.add(featureInput);

        // Create our gencode factory:
        final GencodeFuncotationFactory gencodeFuncotationFactory =
            new GencodeFuncotationFactory(dataSourceFile.resolveSibling(fastaPath),
                    version,
                    transcriptSelectionMode,
                    transcriptList,
                    annotationOverridesMap
            );

        // Add our factory to our factory lists:
        gencodeFuncotationFactories.add( gencodeFuncotationFactory );
        dataSourceFactories.add( gencodeFuncotationFactory );
    }

    // ========================================================================================================
    // Config file static helper methods:
    // ------------------------------------------------

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

    private static void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {

        // Universally required config properties:
        assertConfigPropertiesContainsKey("name", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("version", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("src_file", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("origin_location", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("preprocessing_script", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("type", configFileProperties, configFilePath);

        // Validate our source file:
        assertPathFilePropertiesField( configFileProperties, "src_file", configFilePath);

        // Validate our type:
        final String stringType = configFileProperties.getProperty("type");
        final FuncotatorArgumentDefinitions.DataSourceType type;
        try {
             type = FuncotatorArgumentDefinitions.DataSourceType.getEnum(stringType);
        }
        catch (final IllegalArgumentException ex) {
            throw new UserException.BadInput("ERROR in config file: " + configFilePath.toUri().toString() +
                    " - Invalid value in \"type\" field: " + stringType, ex);
        }

        // Validate remaining values based on type:
        switch (type) {
            case LOCATABLE_XSV: assertLocatableXsvConfigFilePropertiesAreValid(configFileProperties, configFilePath); break;
            case SIMPLE_XSV:    assertSimpleXsvConfigFilePropertiesAreValid(configFileProperties, configFilePath); break;
            case GENCODE:       assertGencodeConfigFilePropertiesAreValid(configFileProperties, configFilePath); break;
            case COSMIC:        /* No special case for COSMIC needed here. */ break;
            default:
                //Note: this should never happen:
                throw new UserException.BadInput("ERROR in config file: " + configFilePath.toUri().toString() +
                        " - Invalid value in \"type\" field: " + stringType);
        }
    }

    private static void assertConfigPropertiesContainsKey(final String key, final Properties configProperties, final Path configFilePath) {
        if ( !configProperties.stringPropertyNames().contains(key) ) {
            throw new UserException.BadInput("Config file for datasource (" + configFilePath.toUri().toString() + ") does not contain required key: " + key);
        }
    }

    private static void assertLocatableXsvConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
        assertConfigPropertiesContainsKey("xsv_delimiter", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("contig_column", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("start_column", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("end_column", configFileProperties, configFilePath);

        // Ensure typed values:
        assertNumericalPropertiesField(configFileProperties, "contig_column", configFilePath);
        assertNumericalPropertiesField(configFileProperties, "start_column", configFilePath);
        assertNumericalPropertiesField(configFileProperties, "end_column", configFilePath);
    }

    private static void assertSimpleXsvConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
        assertConfigPropertiesContainsKey("xsv_delimiter", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("xsv_key", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("xsv_key_column", configFileProperties, configFilePath);
        assertConfigPropertiesContainsKey("xsv_permissive_cols", configFileProperties, configFilePath);

        // Ensure typed values:
        assertNumericalPropertiesField(configFileProperties, "xsv_key_column", configFilePath);
        assertBooleanPropertiesField(configFileProperties, "xsv_permissive_cols", configFilePath);

        // Validate our xsv_key:
        final String stringXsvKey = configFileProperties.getProperty("xsv_key");
        try {
            SimpleKeyXsvFuncotationFactory.XsvDataKeyType.valueOf(stringXsvKey);
        }
        catch (final IllegalArgumentException ex) {
            throw new UserException.BadInput("ERROR in config file: " + configFilePath.toUri().toString() +
                    " - Invalid value in \"xsv_key\" field: " + stringXsvKey, ex);
        }
    }

    private static void assertGencodeConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath){
        assertConfigPropertiesContainsKey("gencode_fasta_path", configFileProperties, configFilePath);

        // Assert that the path is good:
        assertPathFilePropertiesField(configFileProperties, "gencode_fasta_path", configFilePath);
    }

    private static void assertNumericalPropertiesField(final Properties props, final String field, final Path filePath) {
        try {
            Integer.valueOf(props.getProperty(field));
        }
        catch (final NumberFormatException ex) {
            throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                    " - Invalid value in \"" + field + "\" field: " + props.getProperty(field));
        }
    }

    private static void assertBooleanPropertiesField(final Properties props, final String field, final Path filePath) {
        try {
            Boolean.valueOf(props.getProperty(field));
        }
        catch (final NumberFormatException ex) {
            throw new UserException.BadInput("ERROR in config file: " + filePath.toUri().toString() +
                    " - Invalid value in \"" + field + "\" field: " + props.getProperty(field));
        }
    }

    private static void assertPathFilePropertiesField(final Properties props, final String field, final Path filePath) {
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
