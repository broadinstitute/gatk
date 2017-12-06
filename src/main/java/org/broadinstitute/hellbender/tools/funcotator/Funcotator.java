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
import org.broadinstitute.hellbender.tools.funcotator.dataSources.XSV.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.codecs.GENCODE.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
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
            shortName = FuncotatorArgumentDefinitions.GENCODE_FASTA_ARG_NAME,
            fullName  = FuncotatorArgumentDefinitions.GENCODE_FASTA_ARG_NAME,
            doc = "GENCODE Transcript FASTA File.")
    protected File gencodeTranscriptFastaFile;

    @Argument(
            fullName =  FuncotatorArgumentDefinitions.GTF_FILE_ARG_LONG_NAME,
            shortName = FuncotatorArgumentDefinitions.GTF_FILE_ARG_SHORT_NAME,
            doc = "A GENCODE GTF file containing annotated genes."
    )
    protected FeatureInput<GencodeGtfFeature> gtfVariants;

    //-----------------------------------------------------
    // Optional args:

    @Argument(
            shortName = FuncotatorArgumentDefinitions.XSV_INPUT_ARG_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.XSV_INPUT_ARG_LONG_NAME,
            optional = true,
            doc = "Additional data source file (csv/tsv/<DELIMITER>separated value file).  " +
                  "Delimiters for files are specified with the " + FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME + " argument.  " +
                    "If no delimiters are specified, files are assumed to be in Comma Separated Value (CSV) format.  " +
                    "If delimiters are specified, exactly one delimiter must be specified for each additional data source in the order of the data sources."
    )
    protected List<String> xsvInputPaths = FuncotatorArgumentDefinitions.XSV_INPUT_ARG_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_LONG_NAME,
            fullName  = FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_SHORT_NAME,
            optional = true,
            doc = "Delimiter for additional data source file (csv/tsv/<DELIMITER>separated value file)." +
                    "If no delimiters are specified, files are assumed to be in Comma Separated Value (CSV) format."
    )
    protected List<String> xsvInputDelimiters = FuncotatorArgumentDefinitions.XSV_DELIMITER_ARG_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_LONG_NAME,
            optional = true,
            doc = "The column number (0-based) in the XSV file(s) to use as the matching field.  "
    )
    protected List<Integer> xsvKeyColumns = FuncotatorArgumentDefinitions.XSV_KEY_COLUMN_ARG_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_LONG_NAME,
            optional = true,
            doc = "The type of match to perform on each XSV file."
    )
    protected List<SimpleKeyXsvFuncotationFactory.XsvDataKeyType> xsvFileTypes = FuncotatorArgumentDefinitions.XSV_FILE_TYPE_ARG_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.XSV_NAME_ARG_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.XSV_NAME_ARG_LONG_NAME,
            optional = true,
            doc = "The name for each XSV Data Source."
    )
    protected List<String> xsvDataSourceNames = FuncotatorArgumentDefinitions.XSV_NAME_ARG_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_LONG_NAME,
            optional = true,
            doc = "Whether to permissively match the number of columns in the header and data rows."
    )
    protected List<Boolean> xsvPermissiveColumns = FuncotatorArgumentDefinitions.XSV_PERMISSIVE_COLS_ARG_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME,
            optional = true,
            doc = "Method of detailed transcript selection."
    )
    protected FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME,
            optional = true,
            doc = "List of transcripts to use for annotation."
    )
    protected Set<String> transcriptList = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_LONG_NAME,
            optional = true,
            doc = "Default values for annotations that are not added (in the format <ANNOTATION>:<VALUE>)"
    )
    protected List<String> annotationDefaults = FuncotatorArgumentDefinitions.ANNOTATION_DEFAULTS_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_SHORT_NAME,
            fullName  = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_LONG_NAME,
            optional = true,
            doc = "Override values for annotations.  Replaces existing matchihng annotations with given values (in the format <ANNOTATION>:<VALUE>)"
    )
    protected List<String> annotationOverrides = FuncotatorArgumentDefinitions.ANNOTATION_OVERRIDES_DEFAULT_VALUE;

    //==================================================================================================================

    private OutputRenderer outputRenderer;
    private final List<DataSourceFuncotationFactory> dataSourceFactories = new ArrayList<>();
    private GencodeFuncotationFactory gencodeFuncotationFactory;

    //==================================================================================================================

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        logger.info("Beginning traversal");

        final LinkedHashMap<String, String> annotationDefaultsMap = splitAnnotationArgsIntoMap(annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = splitAnnotationArgsIntoMap(annotationOverrides);

        // Set up and add our gencode factory:
        gencodeFuncotationFactory = new GencodeFuncotationFactory(gencodeTranscriptFastaFile,
                                                                 transcriptSelectionMode,
                                                                 transcriptList,
                                                                 annotationOverridesMap);
        dataSourceFactories.add( gencodeFuncotationFactory );

        // Set up our other data source factories:
        // TODO: Set up ALL datasource factories based on config files / directory:
        assertXsvInputsAreValid();

        // Go through and set up each XSV factory:
        for ( int i = 0 ; i < xsvInputPaths.size() ; ++i ) {

            // Create our SimpleKeyXsvFuncotationFactory:
            final SimpleKeyXsvFuncotationFactory factory =
                    //final String name, final Path filePath, final String delim, final int keyColumn, final XsvDataKeyType keyType
                    new SimpleKeyXsvFuncotationFactory(
                            xsvDataSourceNames.get(i),
                            IOUtils.getPath(xsvInputPaths.get(i)),
                            xsvInputDelimiters.get(i),
                            xsvKeyColumns.get(i),
                            xsvFileTypes.get(i),
                            annotationOverridesMap,
                            0,
                            xsvPermissiveColumns.get(i)
                    );

            // Add it to our sources:
            dataSourceFactories.add( factory );
        }

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
        logger.info("Traversal complete!");
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
     * Verifies that the given inputs for XSV files are valid.
     * If they are valid, this method does nothing.
     * If they are invalid, this method throws a {@link UserException.BadInput}.
     */
    private void assertXsvInputsAreValid() {

        if ( !(xsvInputPaths.isEmpty() && xsvInputDelimiters.isEmpty() && xsvDataSourceNames.isEmpty() && xsvFileTypes.isEmpty() && xsvKeyColumns.isEmpty() && xsvPermissiveColumns.isEmpty()) ) {
            if ((xsvInputPaths.size() != xsvInputDelimiters.size() && (!xsvInputDelimiters.isEmpty()))) {
                throw new UserException.BadInput("Must specify the same number of XSV input files and XSV delimiters (or no delimiters to assume CSV format).");
            }
            else if ( (xsvInputPaths.size() != xsvDataSourceNames.size()) && (xsvInputPaths.size() != xsvFileTypes.size()) && (xsvInputPaths.size() != xsvKeyColumns.size()) && (xsvInputPaths.size() != xsvPermissiveColumns.size()) ) {
                throw new UserException.BadInput("Must specify the same number of XSV input arguments for all XSV specifications.");
            }

            // Quick existence checks here:
            for ( final String inputPath : xsvInputPaths ) {
                final Path filePath = IOUtils.getPath(inputPath);

                if ( !Files.exists(filePath) ) {
                    throw new UserException.BadInput("Specified XSV file does not exist: " + filePath.toUri().toString());
                }

                if ( !Files.isReadable(filePath) ) {
                    throw new UserException.BadInput("Cannot read specified XSV file: " + filePath.toUri().toString());
                }

                if ( Files.isDirectory(filePath) ) {
                    throw new UserException.BadInput("Given XSV file path is a directory, but should be a file: " + filePath.toUri().toString());
                }
            }
        }
        else if ( xsvInputPaths.isEmpty() && !(xsvInputDelimiters.isEmpty() && xsvDataSourceNames.isEmpty() && xsvFileTypes.isEmpty() && xsvKeyColumns.isEmpty() && xsvPermissiveColumns.isEmpty()) ) {
            throw new UserException.BadInput("Must specify the same number of XSV input arguments for all XSV specifications.");
        }

    }

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

        final List<Feature> featureList = new ArrayList<>();

        // TODO: Need to add all features here for any kind of Data Source:
        featureList.addAll( featureContext.getValues(gtfVariants) );

        // Create a place to keep our funcotations:
        final List<Funcotation> funcotations = new ArrayList<>();

        // Annotate with Gencode first:
        final List<Funcotation> funcotationsFromGencodeFactory = gencodeFuncotationFactory.createFuncotations(variant, referenceContext, featureList);
        funcotations.addAll( funcotationsFromGencodeFactory );

        // Create a list of GencodeFuncotation to use for other Data Sources:
        final List<GencodeFuncotation> gencodeFuncotations =
                funcotationsFromGencodeFactory.stream()
                    .map(f -> (GencodeFuncotation) f).collect(Collectors.toList());

        // Annotate with the rest of the data sources:
        for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {

            // Make sure we don't double up on the Gencodes:
            if ( funcotationFactory.equals(gencodeFuncotationFactory) ) {
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
}
