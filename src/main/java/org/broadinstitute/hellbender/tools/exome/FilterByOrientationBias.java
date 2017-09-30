package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasFilterer;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.PreAdapterOrientationScorer;
import picard.analysis.artifacts.SequencingArtifactMetrics;
import org.broadinstitute.hellbender.utils.artifacts.Transition;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Additionally filter Mutect2 somatic variant calls for sequence-context dependent artifacts.
 *
 * <p>
 *     This tool is complementary to {@link org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls}.
 *     The tool requires the pre-adapter detailed metrics calculated by Picard CollectSequencingArtifactMetrics.
 *     Specify the base substitution to consider for orientation bias. For a given base substitution specified with
 *     the --artifactModes argument, the tool considers both the forward and reverse complement for filtering.
 * </p>
 *
 * <p>
 *     The metrics from CollectSequencingArtifactMetrics provide a global view across a sample's genome that empowers
 *     decision making in ways site-specific analyses cannot. CollectSequencingArtifactMetrics should be run for both
 *     the normal sample and the tumor sample, if the matched normal is available. The detailed metrics measure orientation
 *     bias for all 3-base contexts and help determine whether variant filtering for a sequence context is necessary.
 * </p>
 *
 * <p>
 *     The most common artifact to consider for filtering is the OxoG (8-Oxoguanine) artifact. OxoG artifacts stem from oxidation of
 *     guanine to 8-oxoguanine, which results in G to T transversion mutations. If samples were derived from histological
 *     slides, then consider the FFPE (formalin-fixed paraffin-embedded) artifact. FFPE artifacts stem from formaldehyde
 *     deamination of cytosines, which results in C to T transition mutations.
 * </p>
 *
 * <h3>Example</h3>
 *
 * For OxoG artifacts, specify G to T artifacts.
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" FilterByOrientationBias \
 *   --artifactModes 'G/T' \
 *   -V tumor_unfiltered.vcf.gz \
 *   -P tumor.pre_adapter_detail_metrics \
 *   --output oxog_filtered.vcf.gz
 * </pre>
 *
 */

@CommandLineProgramProperties(
        summary = "Filter Mutect2 somatic variant calls using the Orientation Bias Filter.\n" +
                "Used for the OxoG (G/T) and Deamination (FFPE) (C/T) artifacts that get introduced into our SNV calling.\n" +
                "\n" +
                "Notes:  All variants are held in RAM.\n This tool will only catch artifacts in diploid organisms.\n" +
                " - Triallelic sites may not be considered for filtering -- the behavior of this tool is undefined for triallelic sites.\n" +
                " - If you do not wish to filter, you must skip the filter entirely.  You cannot do a dummy/non-filter operation with this tool.\n" +
                " - Sites-only VCF files are NOT supported.  At least one sample must be present in the VCF file for any filtering to be done.  If no samples are present in the VCF file, no filtering annotation is actually performed.\n" +
                " - This tool was tested only with output for GATK4 Mutect and makes assumptions about the existence of fields produced by GATK4 Mutect.\n" +
                " - ALT_F1R2 and ALT_F2R1 tags must be present in all SNP variants.\n" +
                " - Do NOT specify artifact modes that are reverse complements of each other.  Behavior of this tool is undefined when that happens.  For example, do not specify C/A and G/T in the same run.\n" +
                " - Any variants that are filtered in the input file are not considered for filtering here, nor are these variants used in deciding cutoff.\n" +
                " - The Orientation Bias puts a filter tag in both the genotype (FORMAT) and variant (FILTER) fields.\n" +
                " - In multiallelic sites, only the first alternate allele is used for filtering.\n" +
                " - This filter should be applied last in any M2 toolchain.\n" +
                " - Common artifacts:\n G/T (OxoG)\n C/T (deamination) ",
        oneLineSummary = "(Experimental) Filter Mutect2 somatic variant calls using orientation bias",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FilterByOrientationBias extends VariantWalker {

    public static final String PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME = "P";
    public static final String PRE_ADAPTER_METRICS_DETAIL_FILE_FULL_NAME = "preAdapterDetailFile";
    public static final String ARTIFACT_MODES_SHORT_NAME = "A";
    public static final String ARTIFACT_MODES_FULL_NAME = "artifactModes";

    public static final String DEFAULT_ARTIFACT_MODE = "G/T";
    public static final String SUMMARY_FILE_SUFFIX = ".summary";

    @Argument(
            doc="Output Somatic SNP/Indel VCF file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    protected File outputFile;

    @Argument(
            doc = "PreAdapter Detail metrics file.  Usually, generated by CollectSequencingArtifactMetrics.",
            shortName = PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME,
            fullName = PRE_ADAPTER_METRICS_DETAIL_FILE_FULL_NAME,
            optional = false)
    protected File preAdapterMetricsFile;

    @Argument(
            doc = "PreAdapter Detail artifacts of interest on the forward strand.  'C/A' for a single artifact.  " +
                    "Specify this parameter multiple times to process multiple artifacts at the same time:  'C/A,T/G'  " +
                    "Artifacts must be one base to one base (e.g. 'CC/CA' is illegal).  G>T is OxoG.",
            shortName = ARTIFACT_MODES_SHORT_NAME,
            fullName = ARTIFACT_MODES_FULL_NAME,
            optional = true
    )
    protected List<String> transitions = new ArrayList<>();

    private Map<Transition, Double> transitionToPreAdapterScoreMap;

    private SortedSet<Transition> relevantTransitions = new TreeSet<>();

    private VariantContextWriter vcfWriter;

    /** Each has an OxoQ annotation */
    private List<VariantContext> firstPassVariants = new ArrayList<>();

    @Override
    public void onTraversalStart() {

        // Gets around issue 2274 in gatk public
        if (transitions.size() == 0) {
            transitions.add(DEFAULT_ARTIFACT_MODE);
        }

        // Sort the input artifacts argument
        transitions.sort(null);

        final MetricsFile<SequencingArtifactMetrics.PreAdapterDetailMetrics, Comparable<?>> mf = new MetricsFile<>();

        try {
            mf.read(new FileReader(preAdapterMetricsFile));
        } catch (final FileNotFoundException fnfe) {
            throw new UserException("Could not find file: " + preAdapterMetricsFile.getAbsolutePath());
        }

        // Collect all of the transitions that were specified in the parameters.
        relevantTransitions.addAll(transitions.stream().map(s -> convertParameterToTransition(s)).collect(Collectors.toSet()));

        // Get the PreAdapterQ score from the picard metrics tool, which gives an indication of how badly infested the bam file is.
        transitionToPreAdapterScoreMap = PreAdapterOrientationScorer.scoreOrientationBiasMetricsOverContext(mf.getMetrics());
        logger.info("preAdapter scores:");
        transitionToPreAdapterScoreMap.keySet().stream().forEach(k -> logger.info(k + ": " + transitionToPreAdapterScoreMap.get(k)));

        setupVCFWriter();
    }

    private Transition convertParameterToTransition(final String stringTransition) {
        final String[] splitTransition = stringTransition.split("/");
        if (!isValidTransition(splitTransition)) {
            throw new UserException("Invalid artifact mode: " + String.join("/", splitTransition));
        }

        return Transition.transitionOf(splitTransition[0].charAt(0), splitTransition[1].charAt(0));
    }

    private boolean isValidTransition(final String[] splitTransition) {
        if (splitTransition[0].length() != 1) {
            throw new UserException.BadInput("First base invalid - must be of length 1: " + splitTransition[0]+ ". Artifact modes must be specified as one base-slash-one base.  E.g. G/T");
        } else if (splitTransition[1].length() != 1) {
            throw new UserException.BadInput("Second base invalid - must be of length 1: " + splitTransition[1]+ ". Artifact modes must be specified as one base-slash-one base.  E.g. G/T");
        }
        return true;
    }

    /**
     *  Just adds the Pre Adapter Q annotation to the variant and creates a new variant.
     *
     *  Note: the writing of the VCF is not done here, since we need to aggregate once we have OxoQ scores.
     *  Note:  No variant is dropped, if no additional annotation was needed, then the original format field is
     *   preserved
     * @param variant See {@link VariantWalker}
     * @param readsContext See {@link VariantWalker}
     * @param referenceContext See {@link VariantWalker}
     * @param featureContext See {@link VariantWalker}
     */
    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final VariantContext updatedVariant = OrientationBiasFilterer.annotateVariantContextWithPreprocessingValues(variant, relevantTransitions, transitionToPreAdapterScoreMap);
        firstPassVariants.add(updatedVariant);

        // See onTraversalSuccess for the actual filtering.
    }

    private void setupVCFWriter() {
        vcfWriter = createVCFWriter(outputFile);
        vcfWriter.writeHeader(OrientationBiasFilterer.createVCFHeader(getHeaderForVariants(), getCommandLine(), transitions));
    }

    @Override
    public Object onTraversalSuccess() {

        logger.info("Tagging whether genotypes are in one of the artifact modes.");

        // Calculate how many artifacts need to be cut
        double fdrThreshold = 0.01;

        final List<VariantContext> finalVariants = OrientationBiasFilterer.annotateVariantContextsWithFilterResults(fdrThreshold, relevantTransitions, firstPassVariants, transitionToPreAdapterScoreMap);

        logger.info("Writing variants to VCF...");
        finalVariants.forEach(vcfWriter::add);

        logger.info("Writing a simple summary table...");
        List<String> sampleNames = new ArrayList<>();
        if (finalVariants.size() != 0) {
            sampleNames = finalVariants.get(0).getSampleNamesOrderedByName();
        }
        final List<Pair<String, Transition>> sampleTransitionCombinations =  new ArrayList<>();
        for (Transition relevantTransition : relevantTransitions) {
            for (String sampleName : sampleNames) {
                sampleTransitionCombinations.add(Pair.of(sampleName, relevantTransition));
            }
        }

        OrientationBiasUtils.writeOrientationBiasSummaryTable(sampleTransitionCombinations, finalVariants, transitionToPreAdapterScoreMap, new File(outputFile.getAbsolutePath() + SUMMARY_FILE_SUFFIX));
        return null;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
