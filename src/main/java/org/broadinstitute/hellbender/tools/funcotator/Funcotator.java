package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.codecs.GENCODE.*;

import java.io.File;
import java.util.*;

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

    public static final String GTF_FILE_ARG_LONG_NAME = "gtfFile";
    public static final String GTF_FILE_ARG_SHORT_NAME = "gtf";
    public static final String GENCODE_FASTA_ARG_NAME = "fasta";

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
            shortName = GENCODE_FASTA_ARG_NAME,
            fullName  = GENCODE_FASTA_ARG_NAME,
            doc = "GENCODE Transcript FASTA File.")
    protected File gencodeTranscriptFastaFile;

    @Argument(
            fullName = GTF_FILE_ARG_LONG_NAME,
            shortName = GTF_FILE_ARG_SHORT_NAME,
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

    //==================================================================================================================

    private OutputRenderer outputRenderer;
    private final List<DataSourceFuncotationFactory> dataSourceFactories = new ArrayList<>();

    //==================================================================================================================

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        logger.info("Beginning traversal");

        // Set up our data source factories:
        dataSourceFactories.add(new GencodeFuncotationFactory(gencodeTranscriptFastaFile, transcriptSelectionMode, transcriptList));

        // Set up our output renderer:
        // TODO: in the future this should be encapsulated into a factory for output renderers based on an input argument.
        outputRenderer = new VcfOutputRenderer(getHeaderForVariants(), createVCFWriter(outputFile), dataSourceFactories);
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
     * Creates an annotation on the given {@code variant} or enqueues it to be processed during a later call to this method.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureContext {@link FeatureContext} corresponding to the given {@code variant}.
     */
    private void enqueueAndHandleVariant(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        final List<Feature> featureList = new ArrayList<>();
        featureList.addAll( featureContext.getValues(gtfVariants) );

        final List<Funcotation> funcotations = new ArrayList<>();
        for ( final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories ) {
            funcotations.addAll( funcotationFactory.createFuncotations(variant, referenceContext, featureList) );
        }
        outputRenderer.write(variant, funcotations);
    }
}
