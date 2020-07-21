package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceSummaryRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Check variants in a VCF against tumor-normal bams representing the same samples, though not the ones from the actual calls.\n" +
                "Bare-bones implementation heavily inspired by MutationValidator from Broad CGA group.\n" +
                "The algorithm is not the same.\n" +
                "This tool can only handle exactly one validation PAIR at a time and this should not be RNA.\n" +
                "Multiallelics in a VCF are not supported and will be skipped.\n" +
                "This tool will validate germline mutations as true positives.\n",
        oneLineSummary = "Check variants against tumor-normal bams representing the same samples, though not the ones from the actual calls.",
        programGroup = VariantEvaluationProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public class ValidateBasicSomaticShortMutations extends VariantWalker {
    public static final String ANNOTATED_VCF_LONG_NAME = "annotated-vcf";
    public static final String SAMPLE_NAME_DISCOVERY_VCF_LONG_NAME = "discovery-sample-name";
    public static final String SAMPLE_NAME_VALIDATION_CASE = "val-case-sample-name";
    public static final String SAMPLE_NAME_VALIDATION_CONTROL = "val-control-sample-name";
    public final static int DEFAULT_MIN_BQ_CUTOFF = 20;
    public final static String CUTOFF_LONG_NAME = "min-base-quality-cutoff";
    public final static String MIN_POWER_LONG_NAME = "min-power";
    public final static String MAX_VALIDATION_NORMAL_COUNT_LONG_NAME = "max-validation-normal-count";

    private static final double DEFAULT_MIN_POWER = 0.9;
    private static final int DEFAULT_MAX_VALIDATION_NORMAL_COUNT = 1;

    @Argument(fullName = SAMPLE_NAME_DISCOVERY_VCF_LONG_NAME,
            doc = "sample name for discovery in VCF.")
    protected String discoverySampleInVcf;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The output file, which will be a validation table (tsv).")
    protected GATKPath outputFile;

    @Argument(fullName = ANNOTATED_VCF_LONG_NAME,
            doc = "Optional output vcf containing original variants annotated with validation info.",
            optional = true)
    protected GATKPath annotatedVcf;

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)",
            fullName = Concordance.SUMMARY_LONG_NAME,
            optional = true)
    protected GATKPath summary;

    @Argument(doc = "Minimum power for an unvalidated variant to be considered a false positive.",
            fullName = MIN_POWER_LONG_NAME,
            optional = true)
    protected double minPower = DEFAULT_MIN_POWER;

    @Argument(doc = "Maximum read count in the validation normal for a variant to validate.  More counts is considered an artifact.",
            fullName = MAX_VALIDATION_NORMAL_COUNT_LONG_NAME,
            optional = true)
    protected int maxValidationNormalCount = DEFAULT_MAX_VALIDATION_NORMAL_COUNT;

    @Argument(fullName = SAMPLE_NAME_VALIDATION_CASE,
            doc = "validation case sample name (in the bam)")
    protected String validationCaseName;

    @Argument(fullName = SAMPLE_NAME_VALIDATION_CONTROL,
            doc = "validation control sample name (in the bam)")
    protected String validationControlName;

    @Argument(fullName = CUTOFF_LONG_NAME,
            doc = "minimum base quality to count a read toward validation.",
            optional = true)
    public int minBqCutoff = DEFAULT_MIN_BQ_CUTOFF;

    private VariantContextWriter vcfWriter;

    @Override
    public boolean requiresReference() {
        return true;
    }

    // for the optional vcf
    public final static String POWER_INFO_FIELD_KEY = "POWER";
    public final static String VALIDATION_AD_INFO_FIELD_KEY = "VAL_AD";
    public final static String JUDGMENT_INFO_FIELD_KEY = "JUDGMENT";

    public final VCFInfoHeaderLine POWER_HEADER_LINE = new VCFInfoHeaderLine(POWER_INFO_FIELD_KEY, 1, VCFHeaderLineType.Float,
            "Power to validate variant in validation bam.");

    public final VCFInfoHeaderLine VALIDATION_AD_HEADER_LINE = new VCFInfoHeaderLine(VALIDATION_AD_INFO_FIELD_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer,
            "Ref and alt allele count in validation bam.");

    public final VCFInfoHeaderLine JUDGMENT_HEADER_LINE = new VCFInfoHeaderLine(JUDGMENT_INFO_FIELD_KEY, 1, VCFHeaderLineType.String,
            "Validation judgment: validated, unvalidated, or skipped.");


    public enum Judgment {
        VALIDATED, UNVALIDATED, SKIPPED;
    }

    private List<BasicValidationResult> results = new ArrayList<>();

    private final MutableInt snpTruePositiveCount = new MutableInt(0);
    private final MutableInt snpFalsePositiveCount = new MutableInt(0);
    private final MutableInt indelTruePositiveCount = new MutableInt(0);
    private final MutableInt indelFalsePositiveCount = new MutableInt(0);

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> result = new ArrayList<>(super.getDefaultReadFilters());
        result.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        result.add(ReadFilterLibrary.NOT_DUPLICATE);
        result.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        return result;
    }

    @Override
    public boolean requiresReads() {return true;}

    @Override
    public void onTraversalStart() {
        if (annotatedVcf != null) {
            final VCFHeader inputHeader = getHeaderForVariants();
            final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
            headerLines.addAll(Arrays.asList(POWER_HEADER_LINE, VALIDATION_AD_HEADER_LINE, JUDGMENT_HEADER_LINE));
            headerLines.addAll(getDefaultToolVCFHeaderLines());
            final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
            vcfWriter = createVCFWriter(annotatedVcf);
            vcfWriter.writeHeader(vcfHeader);
        }
    }

    @Override
    public void apply(VariantContext discoveryVariantContext, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        final Genotype genotype = discoveryVariantContext.getGenotype(discoverySampleInVcf);

        // Skip any symbolic reference alleles
        final Allele referenceAllele = discoveryVariantContext.getReference();
        if (referenceAllele.isSymbolic()) {
            logger.warn("Skipping variant with symbolic reference allele: " + discoveryVariantContext);
            return;
        }

        // If we cannot validate this genotype, we should simply skip it.
        final boolean isAbleToValidate = BasicSomaticShortMutationValidator.isAbleToValidateGenotype(genotype, referenceAllele);
        if (!isAbleToValidate) {
            if (annotatedVcf != null) {
                vcfWriter.add(new VariantContextBuilder(discoveryVariantContext).attribute(JUDGMENT_INFO_FIELD_KEY, Judgment.SKIPPED).make());
            }
            return;
        }

        // We assume that there is only one alternate allele that we are interested in.
        // Multiallelics are not supported.
        final Allele altAllele = genotype.getAllele(1);

        final SAMFileHeader samFileHeader = getHeaderForReads();

        final ReadPileup readPileup = new ReadPileup(discoveryVariantContext, readsContext);
        final Map<String, ReadPileup> pileupsBySample = readPileup.splitBySample(samFileHeader, "__UNKNOWN__");

        // This could happen when read realignment moves reads such that a particular locus has zero reads
        if (pileupsBySample.isEmpty()){
            return;
        }

        final ReadPileup validationNormalPileup = pileupsBySample.get(validationControlName);
        // TODO: handle this case more carefully
        if (validationNormalPileup == null){
            return;
        }

        final ReadPileup validationTumorPileup = pileupsBySample.get(validationCaseName);

        final Map<Allele, MutableInt> validationTumorAllelicCounts = new AllelePileupCounter(referenceAllele, discoveryVariantContext.getAlternateAlleles(), minBqCutoff, validationTumorPileup)
                .getCountMap();

        final int validationTumorAltCount = validationTumorAllelicCounts.getOrDefault(altAllele, new MutableInt(0)).intValue();
        final int validationTumorRefCount = validationTumorAllelicCounts.get(referenceAllele).intValue();

        final String filterString = discoveryVariantContext.getFilters().stream().sorted().collect(Collectors.joining(";"));
        final BasicValidationResult basicValidationResult = BasicSomaticShortMutationValidator.calculateBasicValidationResult(
                genotype, referenceAllele, validationNormalPileup, validationTumorAltCount,
                validationTumorRefCount + validationTumorAltCount, minBqCutoff,
                new SimpleInterval(discoveryVariantContext.getContig(), discoveryVariantContext.getStart(), discoveryVariantContext.getEnd()),
                filterString);

        if (basicValidationResult != null) {
            results.add(basicValidationResult);
        }

        final boolean normalArtifact = basicValidationResult.getNumAltSupportingReadsInNormal() > maxValidationNormalCount;
        final boolean validated = !normalArtifact && basicValidationResult != null && basicValidationResult.isOutOfNoiseFloor();
        final boolean powered = normalArtifact || (basicValidationResult != null && basicValidationResult.getPower() > minPower);

        if (discoveryVariantContext.isSNP()) {
            if (validated) {
                snpTruePositiveCount.increment();
            } else if (powered) {
                snpFalsePositiveCount.increment();
            }
        } else {
            if (validated) {
                indelTruePositiveCount.increment();
            } else if (powered) {
                indelFalsePositiveCount.increment();
            }
        }

        if (annotatedVcf != null) {
            vcfWriter.add(new VariantContextBuilder(discoveryVariantContext)
                    .attribute(JUDGMENT_INFO_FIELD_KEY, validated ? Judgment.VALIDATED : Judgment.UNVALIDATED)
                    .attribute(POWER_INFO_FIELD_KEY, basicValidationResult == null ? 0 : basicValidationResult.getPower())
                    .attribute(VALIDATION_AD_INFO_FIELD_KEY, new int[] {validationTumorRefCount, validationTumorAltCount}).make());
        }
    }

    @Override
    public Object onTraversalSuccess(){
        BasicValidationResult.write(results, outputFile);

        if (summary != null) {
            try (ConcordanceSummaryRecord.Writer concordanceSummaryWriter = ConcordanceSummaryRecord.getWriter(summary)) {
                concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.SNP,
                        snpTruePositiveCount.getValue(),
                        snpFalsePositiveCount.getValue(),
                      0));
                concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.INDEL,
                        indelTruePositiveCount.getValue(),
                        indelFalsePositiveCount.getValue(),
                        0));
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception writing the concordance summary table", e);
            }
        }
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
