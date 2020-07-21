package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureWalker;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Perform functional annotation on a segment file (tsv).  Outputs two files.  The first is a tsv where each row
 * is a segment and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints.
 * The second file has a gene for each row.  The rest of the information is the segment on which it is covered.
 * The first file will have the name as specified by the output parameter.  The second will have '.gene_list.txt' appended.
 * The functionality here is the same as seen in Oncotator with a SEG file as input, but with both SIMPLE_TSV and GENE_LIST as outputs. </p>
 * <p>Note that FuncotateSegments will support seg files from the GATK, whereas Oncotator will not.</p>
 * <p>FuncotateSegments can support hg38 seg files from the GATK.</p>
 * <p>FuncotateSegments does not support direct reading of cloud-based datasources.</p>
 * <p>For more information about Oncotator, on small variants, see Ramos AH, Lichtenstein L, et al. Oncotator: Cancer Variant Annotation Tool. Human Mutation (2015). http://dx.doi.org/10.1002/humu.22771</p>
 * <h3>Usage example</h3>
 * <pre>
 *   ./gatk FuncotateSegments \
 *    --data-sources-path dataSourcesFolder/ \
 *   --ref-version hg19 \
 *   --output-file-format SEG \
 *   -R reference.fasta \
 *   --segments input.seg \
 *   -O input.seg.funcotated.tsv \
 *   --transcript-list tx_list.txt
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Perform functional annotation on a segment file (tsv).  Outputs two files.  The first is a tsv where each row" +
                " is a segment and the annotations are the covered genes and which genes+exon is overlapped by the segment breakpoints." +
                "  The second file has a gene for each row.  The rest of the information is the segment on which it is covered.\n" +
                "  The first file will have the name as specified by the output parameter.  The second will have '.gene_list.txt' appended.\n" +
                "The functionality here is the same as seen in Oncotator with a SEG file as input, but with both SIMPLE_TSV and GENE_LIST as outputs.  " +
                "\nNote that FuncotateSegments will support seg files from the GATK, whereas Oncotator will not." +
                "\nFuncotateSegments can support hg38 seg files from the GATK." +
                "\nFuncotateSegments does not support direct reading of cloud-based datasources." +
                "\nFor more information about Oncotator, on small variants, see Ramos AH, Lichtenstein L, et al. Oncotator: Cancer Variant Annotation Tool. Human Mutation (2015). http://dx.doi.org/10.1002/humu.22771",
        oneLineSummary = "Functional annotation for segment files.  The output formats are not well-defined and subject to change.",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FuncotateSegments extends FeatureWalker<AnnotatedInterval> {
    private static final Logger logger = LogManager.getLogger(FuncotateSegments.class);

    private static final List<String> MAPPING_DEFAULT = Arrays.asList(
            "MEAN_LOG2_COPY_RATIO:Segment_Mean",
            "CALL:Segment_Call",
            "sample:Sample",
            "sample_id:Sample",
            "NUM_POINTS_COPY_RATIO:Num_Probes"
            );

    private static final String MAPPING_FULL_NAME = "alias-to-key-mapping";

    @Argument(
            doc = "Input segment file (tab-separated values).  Must have a call column.",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME
    )
    private GATKPath segmentFile;

    @ArgumentCollection
    private final FuncotatorSegmentArgumentCollection funcotatorArgs = new FuncotatorSegmentArgumentCollection();

    /** Mapping goes from old name to new name.  Note that this is not necessarily what is output, since the output
        renderer might further transform the name.*/
    @Argument(doc="(Advanced) Mapping between an alias and key values that are recognized by the backend.  Users should not typically have to specify this.",
            fullName = MAPPING_FULL_NAME)
    private List<String> aliasToKeyMappingAsString = new ArrayList<>(MAPPING_DEFAULT);

    private LinkedHashMap<String, String> aliasToKeyMapping;

    private FuncotatorEngine funcotatorEngine;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(AnnotatedInterval.class);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private OutputRenderer outputRenderer;

    @Override
    public void onTraversalStart() {

        Utils.validateArg(funcotatorArgs.transcriptSelectionMode != TranscriptSelectionMode.ALL, "Cannot funcotate segments with the ALL transcript selection mode.  Please select another mode.");

        // Get our overrides for annotations:
        final LinkedHashMap<String, String> annotationDefaultsMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationOverrides);

        // Next set up our transcript list:
        final Set<String> finalUserTranscriptIdSet = FuncotatorEngine.processTranscriptList(funcotatorArgs.userTranscriptIdSet);

        // Initialize the datasources (and make sure to filter to handle only segment-enabled funcotation factories.
        // Initialize all of our data sources:
        // Sort data sources to make them process in the same order each time:
        funcotatorArgs.dataSourceDirectories.sort(Comparator.naturalOrder());
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths(funcotatorArgs.referenceVersion, funcotatorArgs.dataSourceDirectories);

        // Create the data sources from the input:
        // This will also create and register the FeatureInputs (created by the Data Sources)
        // with the GATK Engine, so we do not have to plumb them in after the fact.
        //  Only take datasources that support the annotation of segments.
        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                annotationOverridesMap,
                funcotatorArgs.transcriptSelectionMode,
                finalUserTranscriptIdSet,
                this,
                funcotatorArgs.lookaheadFeatureCachingInBp,
                new FlankSettings(0,0),
                true,
                funcotatorArgs.minNumBasesForValidSegment
        ).stream()
         .filter(DataSourceFuncotationFactory::isSupportingSegmentFuncotation)
         .collect(Collectors.toList());

        // Log the datasources
        logger.info("The following datasources support funcotation on segments: ");
        dataSourceFuncotationFactories.forEach(ff -> logger.info(" " + ff.getInfoString()));

        // Initialize a funcotator engine to handle segments.
        funcotatorEngine = new FuncotatorEngine(funcotatorArgs,
                getBestAvailableSequenceDictionary(),
                VcfFuncotationMetadata.create(
                        Arrays.asList(
                                new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                                        VCFHeaderLineType.Integer, "End coordinate of the variant")
                        )
                ),
                dataSourceFuncotationFactories
        );

        aliasToKeyMapping = FuncotatorEngine.splitAnnotationArgsIntoMap(aliasToKeyMappingAsString);

        // Create a composite output renderer -- happens in the engine as long as the output format is SEG.
        outputRenderer = funcotatorEngine.createOutputRenderer(annotationDefaultsMap, annotationOverridesMap,
                            new VCFHeader(), getDefaultToolVCFHeaderLines(), this);
    }

    @Override
    public void apply(final AnnotatedInterval segment, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // Convert segment into a VariantContext while honoring the funcotation engine's necessary conversions.
        final VariantContext segmentVariantContext = funcotatorEngine.getDefaultVariantTransformer()
                .andThen(getTransformAttributesToStandardNames())
                .apply(AnnotatedIntervalToSegmentVariantContextConverter.convert(segment, referenceContext));

        // Get the correct reference for B37/HG19 compliance:
        // This is necessary because of the variant transformation that gets applied in VariantWalkerBase::apply.
        final ReferenceContext correctReferenceContext = funcotatorEngine.getCorrectReferenceContext(segmentVariantContext, referenceContext);

        // funcotate
        //  The resulting funcotation map should only have one transcript ID (which is the "no transcript" ID).
        final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForSegment(segmentVariantContext, correctReferenceContext, featureContext);

        // This will propagate input variant context attributes to the output
        for (final String txId : funcotationMap.getTranscriptList()) {
            funcotationMap.add(txId, FuncotatorUtils.createFuncotationsFromMetadata(segmentVariantContext, createMetadata(), "FUNCOTATE_SEGMENTS"));
        }

        // Force the final output to have the same contig convention as the input.
        final VariantContext finalVC = new VariantContextBuilder(segmentVariantContext)
                .chr(segment.getContig())
                .make();

        // write the variant context
        outputRenderer.write(finalVC, funcotationMap);
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return segmentFile;
    }

    @Override
    public Object onTraversalSuccess() {
        return true;
    }

    @Override
    public void closeTool() {
        if ( funcotatorEngine != null) {
            funcotatorEngine.close();
        }

        if ( outputRenderer != null ) {
            outputRenderer.close();
        }
    }

    @Override
    protected <T extends Feature> SimpleInterval makeFeatureInterval(final T feature) {
        if (funcotatorArgs.referenceVersion.equals(BaseFuncotatorArgumentCollection.FuncotatorReferenceVersionHg19)) {
            return new SimpleInterval(FuncotatorUtils.convertB37ContigToHg19Contig(feature.getContig()), feature.getStart(), feature.getEnd());
        } else {
            return new SimpleInterval(feature);
        }
    }

    private FuncotationMetadata createMetadata() {
        return VcfFuncotationMetadata.create(
                Arrays.asList(
                        new VCFInfoHeaderLine("Segment_Mean",1, VCFHeaderLineType.Float, "Mean for the segment.  Units will be the same as the input file."),
                        new VCFInfoHeaderLine("Num_Probes",1, VCFHeaderLineType.Integer, "Number of probes/targets/bins overlapping the segment."),
                        new VCFInfoHeaderLine("Segment_Call",1, VCFHeaderLineType.String, "Segment call (whether the segment is amplified, deleted, etc)."),
                        new VCFInfoHeaderLine("Sample",1, VCFHeaderLineType.String, "Sample name for the segment."),
                        new VCFInfoHeaderLine("build",1, VCFHeaderLineType.String, "Genome build (e.g. 'hg19' or 'hg38').")
                )
        );
    }

    private VariantContext transformAttributesToStandardNames(final VariantContext vc) {

        final Map<String,Object> transformedAttributes = vc.getAttributes().entrySet().stream()
                .collect(Collectors.toMap(e-> aliasToKeyMapping.getOrDefault(e.getKey(), e.getKey()), e -> e.getValue()));
        final VariantContextBuilder vcb = new VariantContextBuilder(vc).attributes(transformedAttributes);
        return vcb.make();
    }

    private VariantTransformer getTransformAttributesToStandardNames() {
        return this::transformAttributesToStandardNames;
    }
}
