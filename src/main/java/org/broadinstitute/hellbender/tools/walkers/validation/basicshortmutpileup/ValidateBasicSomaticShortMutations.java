package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@CommandLineProgramProperties(

        summary = "(Experimental) Bare-bones implementation heavily inspired by MutationValidator from Broad CGA group.\n" +
                "The algorithm is not the same.\n" +
                "This tool can only handle exactly one validation PAIR at a time and this should not be RNA.\n" +
                "Multiallelics in a VCF are not supported and will be skipped.\n" +
                "This tool will validate germline mutations as true positives.\n",
        oneLineSummary = "(Experimental) Check the variants in a VCF against a tumor-normal pair of bams representing the same samples, though not the ones from the actual calls.",
        programGroup = VariantProgramGroup.class
)
@BetaFeature
public class ValidateBasicSomaticShortMutations extends VariantWalker {
    public static final String SAMPLE_NAME_DISCOVERY_VCF = "dsv";
    public static final String SAMPLE_NAME_VALIDATION_CASE = "valcase";
    public static final String SAMPLE_NAME_VALIDATION_CONTROL = "valcontrol";
    public final static int DEFAULT_MIN_BQ_CUTOFF = 20;
    public final static String CUTOFF_SHORT_NAME = "bqcutoff";
    public final static String CUTOFF_LONG_NAME = "min-base-quality-cutoff";

    @Argument(shortName = SAMPLE_NAME_DISCOVERY_VCF,
            doc = "sample name for discovery in VCF.")
    protected String discoverySampleInVcf;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "The output file, which will be a validation table (tsv).")
    protected String outputFile;

    @Argument(shortName = SAMPLE_NAME_VALIDATION_CASE,
            doc = "validation case sample name (in the bam)")
    protected String validationCaseName;

    @Argument(shortName = SAMPLE_NAME_VALIDATION_CONTROL,
            doc = "validation control sample name (in the bam)")
    protected String validationControlName;

    @Argument(shortName = CUTOFF_SHORT_NAME,
            fullName = CUTOFF_LONG_NAME,
            doc = "minimum base quality to count a read toward validation.",
            optional = true)
    public Integer minBqCutoff = DEFAULT_MIN_BQ_CUTOFF;

    @Override
    public boolean requiresReference() {
        return true;
    }

    public final static String CONTIG = SimpleAnnotatedGenomicRegion.CONTIG_HEADER;
    public final static String START = SimpleAnnotatedGenomicRegion.START_HEADER;
    public final static String END = SimpleAnnotatedGenomicRegion.END_HEADER;
    public final static String REF = "ref_allele";
    public final static String ALT = "alt_allele";
    public final static String DISCOVERY_ALT_COVERAGE = "t_alt_count";
    public final static String DISCOVERY_REF_COVERAGE = "t_ref_count";
    public final static String VALIDATION_ALT_COVERAGE = "tv_alt_count";
    public final static String VALIDATION_REF_COVERAGE = "tv_ref_count";
    public final static String MIN_VAL_COUNT = "min_val_count";
    public final static String POWER = "power";
    public final static String IS_NOT_NOISE = "validated";
    public final static String IS_ENOUGH_VALIDATION_COVERAGE = "sufficient_tv_alt_coverage";
    public final static String DISCOVERY_VCF_FILTER = "discovery_vcf_filter";

    public static String[] headers = {CONTIG, START, END, REF, ALT, DISCOVERY_ALT_COVERAGE, DISCOVERY_REF_COVERAGE,
            VALIDATION_ALT_COVERAGE, VALIDATION_REF_COVERAGE, MIN_VAL_COUNT, POWER, IS_NOT_NOISE, IS_ENOUGH_VALIDATION_COVERAGE, DISCOVERY_VCF_FILTER};

    private List<BasicValidationResult> results = new ArrayList<>();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> result = new ArrayList<>(super.getDefaultReadFilters());
        result.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        result.add(ReadFilterLibrary.NOT_DUPLICATE);
        return result;
    }

    @Override
    public boolean requiresReads() {return true;}

    @Override
    public void apply(VariantContext discoveryVariantContext, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        final Genotype genotype = discoveryVariantContext.getGenotype(discoverySampleInVcf);

        // If we cannot validate this genotype, we should simple skip it.
        final Allele referenceAllele = discoveryVariantContext.getReference();
        final boolean isAbleToValidate = BasicSomaticShortMutationValidator.isAbleToValidateGenotype(genotype, referenceAllele);
        if (!isAbleToValidate) {
            return;
        }

        // We assume that there is only one alternate allele that we are interested in.
        //   Please note that multiallelics are not supported.
        final Allele altAllele = genotype.getAllele(1);

        final SAMFileHeader samFileHeader = getHeaderForReads();

        final ReadPileup readPileup = GATKProtectedVariantContextUtils.getPileup(discoveryVariantContext, readsContext);
        final Map<String, ReadPileup> stringPileupElementMap = readPileup.splitBySample(samFileHeader, "__UNKNOWN__");

        final ReadPileup validationNormalPileup = stringPileupElementMap.get(validationControlName);

        final ReadPileup validationTumorPileup = stringPileupElementMap.get(validationCaseName);

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

        results.add(basicValidationResult);
    }

    @Override
    public boolean requiresIntervals() {return true;}

    @Override
    public Object onTraversalSuccess(){
        final TableColumnCollection tableColumnCollection = new TableColumnCollection(headers);
        try (final TableWriter<BasicValidationResult> writer = new TableWriter<BasicValidationResult>(new File(outputFile), tableColumnCollection) {
            @Override
            protected void composeLine(BasicValidationResult record, DataLine dataLine) {
                dataLine.set(CONTIG, record.getContig());
                dataLine.set(START, record.getStart());
                dataLine.set(END, record.getEnd());
                dataLine.set(REF, record.getReference().getBaseString());
                dataLine.set(ALT, record.getAlternate().getBaseString());
                dataLine.set(DISCOVERY_ALT_COVERAGE, record.getDiscoveryAltCount());
                dataLine.set(DISCOVERY_REF_COVERAGE, record.getDiscoveryRefCount());
                dataLine.set(VALIDATION_ALT_COVERAGE, record.getValidationAltCount());
                dataLine.set(VALIDATION_REF_COVERAGE, record.getValidationRefCount());
                dataLine.set(MIN_VAL_COUNT, record.getMinValidationReadCount());
                dataLine.set(POWER, record.getPower());
                dataLine.set(IS_NOT_NOISE, record.isOutOfNoiseFloor());
                dataLine.set(IS_ENOUGH_VALIDATION_COVERAGE, record.isEnoughValidationReads());
                dataLine.set(DISCOVERY_VCF_FILTER, record.getFilters() == null ? "" : record.getFilters());
            }
        }) {
            writer.writeHeaderIfApplies();
            writer.writeAllRecords(results);
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile(new File(outputFile), "Could not create file: " + new File(outputFile).getAbsolutePath());
        }
        return "SUCCESS";
    }
}
