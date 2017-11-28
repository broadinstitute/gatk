package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.codecs.GENCODE.GencodeGtfFeature;

import java.io.File;
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
    private FeatureInput<GencodeGtfFeature> gtfVariants;

    //-----------------------------------------------------
    // Optional args:

    @Argument(
            shortName = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_LONG_NAME,
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_SHORT_NAME,
            optional = true,
            doc = "Method of detailed transcript selection."
    )
    protected FuncotatorArgumentDefinitions.TranscriptSelectionMode transcriptSelectionMode = FuncotatorArgumentDefinitions.TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE;

    @Argument(
            shortName = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_LONG_NAME,
            fullName  = FuncotatorArgumentDefinitions.TRANSCRIPT_LIST_SHORT_NAME,
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

        // Set up our gencode factory:
        gencodeFuncotationFactory = new GencodeFuncotationFactory(gencodeTranscriptFastaFile,
                                                                 transcriptSelectionMode,
                                                                 transcriptList,
                                                                 annotationOverridesMap);

        // Set up our data source factories:
        // TODO: Set up ALL datasource factories based on config file:
        dataSourceFactories.add(
                gencodeFuncotationFactory
        );

        // Need to determine which annotations are accounted for (by the funcotation factories) and which are not.
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
            for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {
                if ( !funcotationFactory.getSupportedFuncotationFields().contains(field) ) {
                    outAnnotations.put(field, annotationMap.get(field));
                }
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
