package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Calculates copy number posteriors for a given set of structural variants.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 *     <li>
 *         Germline CNV intervals VCF
 *     </li>
 *     <li>
 *         Germline CNV contig ploidy calls
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF with copy number posteriors
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCopyNumberPosteriors
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Collects read counts at specified intervals",
        oneLineSummary = "Collects read counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class SVCopyNumberPosteriors extends VariantWalker {
    public static final String COPY_NUMBER_INTERVALS_LONG_NAME = "cnv-intervals-vcf";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "ploidy-calls-file";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String COPY_NEUTRAL_PRIOR_LONG_NAME = "copy-neutral-prior";
    public static final String CNV_BND_PROB_THRESHOLD_LONG_NAME = "cnv-bnd-max-phred";
    public static final String CNV_BND_SAMPLE_THRESHOLD_LONG_NAME = "cnv-bnd-carrier-fraction";
    public static final String GENOTYPE_DEPTH_CALLS_LONG_NAME = "genotype-depth-calls";

    @Argument(
            doc = "Germline copy number intervals VCF. Can be specified more than once for runs with different bin sizes.",
            fullName = COPY_NUMBER_INTERVALS_LONG_NAME
    )
    private List<File> copyNumberPosteriorsFiles;

    @Argument(
            doc = "Contig ploidy calls file. Can be specified for multiple samples.",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME
    )
    private List<File> contigPloidyCallFiles;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Argument(
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 0;

    @Argument(
            doc = "Prior probability of neutral copy number for " + DepthEvidenceAggregator.COPY_NEUTRAL_PRIOR_BASIS_LENGTH +
                    "bp in regions where copy number calls are not available.",
            fullName = COPY_NEUTRAL_PRIOR_LONG_NAME,
            minValue = 0,
            maxValue = 1
    )
    private double copyNeutralPrior = 0.99;

    @Argument(
            doc = "Max phred-scaled alt copy state probability for CNV carrier genotypes.",
            fullName = CNV_BND_PROB_THRESHOLD_LONG_NAME,
            minValue = 0
    )
    private int cnvBndMaxStatePhred = 20;

    @Argument(
            doc = "Min carrier sample fraction meeting --" + CNV_BND_PROB_THRESHOLD_LONG_NAME
                    + " for conversion of CNVs to BNDs. Set to 0 to bypass this step.",
            fullName = CNV_BND_SAMPLE_THRESHOLD_LONG_NAME,
            minValue = 0.,
            maxValue = 1.
    )
    private double cnvBndMinSampleFraction = 0.5;

    @Argument(
            doc = "Perform genotyping on depth-only calls.",
            fullName = GENOTYPE_DEPTH_CALLS_LONG_NAME
    )
    private boolean genotypeDepthCalls = false;

    @Argument(
            doc = "Depth-only call min included intervals overlap",
            fullName = SVCluster.DEPTH_ONLY_INCLUDE_INTERVAL_OVERLAP_LONG_NAME,
            minValue = 0,
            maxValue = 1,
            optional = true
    )
    private double minDepthOnlyIncludeOverlap = 0.5;

    private final Map<String,IntervalTree<Object>> includedIntervalsTreeMap = new HashMap<>();
    private List<VCFFileReader> posteriorsReaders;
    private List<String> samples;
    private VariantContextWriter outputWriter;
    private String currentContig;
    private SVGenotypeEngineDepthOnly depthOnlyGenotypeEngine;
    private DepthEvidenceAggregator depthEvidenceAggregator;
    private SAMSequenceDictionary dictionary;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        if (genotypeDepthCalls) {
            depthOnlyGenotypeEngine = new SVGenotypeEngineDepthOnly();
        }
        posteriorsReaders = copyNumberPosteriorsFiles.stream().map(VCFFileReader::new).collect(Collectors.toList());
        samples = getHeaderForVariants().getSampleNamesInOrder();
        final Collection<CalledContigPloidyCollection> contigPloidyCollections = contigPloidyCallFiles.stream()
                .map(CalledContigPloidyCollection::new)
                .collect(Collectors.toList());
        validateSampleSets(contigPloidyCollections);
        depthEvidenceAggregator = new DepthEvidenceAggregator(posteriorsReaders, contigPloidyCollections, copyNeutralPrior, samples, dictionary);
        outputWriter = createVCFWriter(outputFile);
        outputWriter.writeHeader(composeHeader());
        currentContig = null;
        loadIntervalTree();
    }

    @Override
    public Object onTraversalSuccess() {
        outputWriter.close();
        return null;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord call = SVCallRecordUtils.create(variant);
        if (!SVCallRecordUtils.isValidSize(call, minEventSize)
                || !SVCallRecordUtils.intervalIsIncluded(call, includedIntervalsTreeMap, minDepthOnlyIncludeOverlap)) {
            return;
        }

        VariantContext finalVariant = depthEvidenceAggregator.apply(call, variant);
        if (SVGenotypeEngine.CNV_TYPES.contains(finalVariant.getStructuralVariantType())
                && !SVGenotypeEngineFromModel.isDepthOnlyVariant(finalVariant)
                && !cnvHasDepthSupport(finalVariant, cnvBndMinSampleFraction, cnvBndMaxStatePhred)) {
            finalVariant = convertCnvToBnd(finalVariant);
        }
        if (genotypeDepthCalls && SVGenotypeEngineFromModel.isDepthOnlyVariant(finalVariant)) {
            finalVariant = depthOnlyGenotypeEngine.genotypeVariant(finalVariant);
        }
        outputWriter.add(finalVariant);
    }

    private void validateSampleSets(final Collection<CalledContigPloidyCollection> contigPloidyCollections) {
        final Set<String> svSamples = new HashSet<>(samples);
        final Set<String> cnvSamples = getCNVSamples();
        final Set<String> contigPloidySamples = getContigPloidySamples(contigPloidyCollections);
        validateSubset(cnvSamples, svSamples, "The following samples from the SV VCF were not in the CNV VCF");
        validateSubset(contigPloidySamples, svSamples, "The following samples from the SV VCF were not in the contig ploidy calls");
    }

    private void validateSubset(final Set<String> set, final Set<String> subset, final String msg) {
        final Set<String> diff = Sets.difference(subset, set);
        Utils.validate(diff.isEmpty(), msg + ":" + String.join(", ", diff));
    }

    private Set<String> getContigPloidySamples(final Collection<CalledContigPloidyCollection> contigPloidyCollections) {
        return contigPloidyCollections.stream()
                .map(p -> p.getMetadata().getSampleName())
                .collect(Collectors.toSet());
    }

    private Set<String> getCNVSamples() {
        final Set<String> cnvSamplesSet = new HashSet<>(posteriorsReaders.get(0).getFileHeader().getSampleNamesInOrder());
        for (int i = 1; i < posteriorsReaders.size(); i++) {
            final Set<String> set = new HashSet<>(posteriorsReaders.get(i).getFileHeader().getSampleNamesInOrder());
            if (!set.equals(cnvSamplesSet)) {
                throw new UserException.BadInput("CNV VCFs do not contain identical samples.");
            }
        }
        return Collections.unmodifiableSet(cnvSamplesSet);
    }

    private void loadIntervalTree() {
        final List<SimpleInterval> intervals = getRequestedIntervals();
        if (intervals == null) {
            throw new UserException.MissingReference("Reference dictionary is required");
        }
        for (final SimpleInterval interval : intervals) {
            includedIntervalsTreeMap.putIfAbsent(interval.getContig(), new IntervalTree<>());
            includedIntervalsTreeMap.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), null);
        }
    }

    private VCFHeader composeHeader() {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.addAll(getDefaultToolVCFHeaderLines());
        headerInfo.add(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, 1,
                VCFHeaderLineType.Integer, "Phred-scaled copy number posterior over the event region"));
        headerInfo.add(new VCFFormatHeaderLine(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY, 1,
                VCFHeaderLineType.Integer, "Neutral copy number"));
        if (genotypeDepthCalls) {
            headerInfo.addAll(depthOnlyGenotypeEngine.getHeaderLines());
        }
        final VCFHeader vcfHeader = new VCFHeader(getHeaderForVariants());
        headerInfo.stream().forEach(line -> vcfHeader.addMetaDataLine(line));
        return vcfHeader;
    }

    private static VariantContext convertCnvToBnd(final VariantContext variant) {
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        builder.alleles(Lists.newArrayList(Allele.REF_N, SVGenotypeEngine.BND_SYMBOLIC_ALLELE));
        builder.attribute(VCFConstants.SVTYPE, StructuralVariantType.BND);
        builder.attribute(GATKSVVCFConstants.SVLEN, -1);
        return builder.make();
    }

    protected static boolean cnvHasDepthSupport(final VariantContext variant,
                                                final double cnvBndMinSampleFraction,
                                                final double cnvBndMaxStatePhred) {
        final StructuralVariantType svType = variant.getStructuralVariantType();
        if (!SVGenotypeEngine.CNV_TYPES.contains(svType)) {
            throw new IllegalArgumentException("Variant was not a CNV");
        }
        if (cnvBndMinSampleFraction == 0) {
            return true;
        }
        int samplesWithDepthSupport = 0;
        for (final Genotype genotype : variant.getGenotypes()) {
            if (genotypeHasDepthSupport(genotype, svType, cnvBndMaxStatePhred)) {
                samplesWithDepthSupport++;
            }
        }
        final int numCarriers = (int) variant.getGenotypes().stream()
                .filter(g -> VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE) == GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE)
                .count();
        final int minSamples = (int) Math.ceil(numCarriers * cnvBndMinSampleFraction);
        return samplesWithDepthSupport >= minSamples;
    }

    private static boolean genotypeHasDepthSupport(final Genotype genotype,
                                                   final StructuralVariantType svType,
                                                   final double cnvBndMaxStatePhred) {
        final CopyNumberPosteriorDistribution copyStatePosterior = SVGenotypeEngineDepthOnly.getCopyNumberStatePosterior(genotype);
        final List<IntegerCopyNumberState> states = copyStatePosterior.getIntegerCopyNumberStateList();
        final int samplePloidy = SVGenotypeEngineFromModel.getNeutralCopyNumber(genotype);
        if (svType.equals(StructuralVariantType.DEL)) {
            for (int i = 0; i < samplePloidy; i++) {
                if (QualityUtils.logProbToPhred(copyStatePosterior.getCopyNumberPosterior(states.get(i))) <= cnvBndMaxStatePhred) {
                    return true;
                }
            }
        } else {
            // DUP or CNV
            for (int i = samplePloidy + 1; i < states.size(); i++) {
                if (QualityUtils.logProbToPhred(copyStatePosterior.getCopyNumberPosterior(states.get(i))) <= cnvBndMaxStatePhred) {
                    return true;
                }
            }
        }
        return false;
    }
}
