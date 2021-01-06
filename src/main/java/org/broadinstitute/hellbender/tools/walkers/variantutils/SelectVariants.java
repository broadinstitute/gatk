package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;

import java.nio.file.Path;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBArgumentCollection;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.samples.MendelianViolation;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.*;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Select a subset of variants from a VCF file
 *
 * <p>This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain
 * analyses. Examples of such analyses include comparing and contrasting cases vs. controls, extracting variant or
 * non-variant loci that meet certain requirements, or troubleshooting some unexpected results, to name a few.</p>
 *
 * <p>
 * There are many different options for selecting subsets of variants from a larger callset:
 * <ul>
 *     <li>Extract one or more samples from a callset based on either a complete sample name or a pattern match.</li>
 *     <li>Specify criteria for inclusion that place thresholds on annotation values, e.g. "DP > 1000" (depth of
 *     coverage greater than 1000x), "AF < 0.25" (sites with allele frequency less than 0.25). These criteria are written
 *     as "JEXL expressions", which are documented in the
 *     <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions">article about using JEXL expressions</a>.</li>
 *     <li>Provide concordance or discordance tracks in order to include or exclude variants that are also present
 *     in other given callsets.</li>
 *     <li>Select variants based on criteria like their type (e.g. INDELs only), evidence of mendelian violation,
 *     filtering status, allelicity, etc.</li>
 * </ul>
 * </p>
 *
 * <p>There are also several options for recording the original values of certain annotations which are recalculated
 * when one subsets the new callset, trims alleles, etc.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant call set in VCF format from which a subset can be selected.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A new VCF file containing the selected subset of variants.
 * </p>
 *
 * * <h3>Usage examples</h3>
 * <h4>Select SNPs</h4>
 * <pre>
 *     gatk SelectVariants \
 *     -R Homo_sapiens_assembly38.fasta \
 *     -V input.vcf \
 *     --select-type-to-include SNP \
 *     -O output.vcf
 * </pre>
 *
 * <h4>Query Chromosome 20 Variants from a GenomicsDB</h4>
 * <pre>
 *     gatk SelectVariants \
 *     -R Homo_sapiens_assembly38.fasta \
 *     -V gendb://genomicsDB \
 *     -L 20 \
 *     -O output.chr20.vcf
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain " +
                  "analyses. Examples include comparing and contrasting cases vs. controls, extracting variant or non-variant loci that meet " +
                  "certain requirements, or troubleshooting some unexpected results, to name a few.",
        oneLineSummary = "Select a subset of variants from a VCF file",
        programGroup = VariantManipulationProgramGroup.class
)
@DocumentedFeature
public final class SelectVariants extends VariantWalker {

    private static final int MAX_FILTERED_GENOTYPES_DEFAULT_VALUE  = Integer.MAX_VALUE;
    private static final double MAX_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE = 1.0;

    private static final int MIN_FILTERED_GENOTYPES_DEFAULT_VALUE  = 0;
    private static final double MIN_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE = 0.0;

    private static final int MAX_NOCALL_NUMBER_DEFAULT_VALUE = Integer.MAX_VALUE;
    private static final double MAX_NOCALL_FRACTION_DEFAULT_VALUE = 1.0;

    /**
     * A site is considered discordant if there exists some sample in the variant track that has a non-reference
     * genotype and either the site isn't present in this track, the sample isn't present in this track, or the
     * sample is called reference in this track.
     */
    @Argument(fullName="discordance", shortName = "disc",
                    doc="Output variants not called in this comparison track", optional=true)
    private FeatureInput<VariantContext> discordanceTrack;

    /**
     * A site is considered concordant if (1) we are not looking for specific samples and there is a variant called
     * in both the variant and concordance tracks or (2) every sample present in the variant track is present in the
     * concordance track and they have the sample genotype call.
     */
    @Argument(fullName="concordance", shortName = "conc",
                    doc="Output variants also called in this comparison track", optional=true)
    private FeatureInput<VariantContext> concordanceTrack;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              doc="Path to which variants should be written")
    public GATKPath vcfOutput = null;

    /**
     * This argument can be specified multiple times in order to provide multiple sample names, or to specify
     * the name of one or more files containing sample names. File names must use the extension ".args", and the
     * expected file format is simply plain text with one sample name per line. Note that sample exclusion takes
     * precedence over inclusion, so that if a sample is in both lists it will be excluded.
     */
    @Argument(fullName=StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME, shortName=StandardArgumentDefinitions.SAMPLE_NAME_SHORT_NAME, doc="Include genotypes from this sample", optional=true)
    private Set<String> sampleNames = new LinkedHashSet<>(0);

    /**
     * Using a regular expression allows you to match multiple sample names that have that pattern in common. This
     * argument can be specified multiple times in order to use multiple different matching patterns.
     */
    @Argument(fullName="sample-expressions", shortName="se",
                    doc="Regular expression to select multiple samples", optional=true)
    private Set<String> sampleExpressions = new LinkedHashSet<>(0);

    /**
     * Note that sample exclusion takes precedence over inclusion, so that if a sample is in both lists it will be
     * excluded. This argument can be specified multiple times in order to provide multiple sample names, or to
     * specify the name of one or more files containing sample names. File names must use the extension ".args",
     * and the expected file format is simply plain text with one sample name per line.
     */
    @Argument(fullName="exclude-sample-name", shortName="xl-sn", doc="Exclude genotypes from this sample", optional=true)
    private Set<String> XLsampleNames = new LinkedHashSet<>(0);

    /**
     * Using a regular expression allows you to match multiple sample names that have that pattern in common. Note that
     * sample exclusion takes precedence over inclusion, so that if a sample is in both lists it will be excluded. This
     * argument can be specified multiple times in order to use multiple different matching patterns.
     */
    @Argument(fullName="exclude-sample-expressions", shortName="xl-se",
                    doc="List of sample expressions to exclude", optional=true)
    private Set<String> XLsampleExpressions = new LinkedHashSet<>(0);

    /**
     * See example commands above for detailed usage examples. Note that these expressions are evaluated *after* the
     * specified samples are extracted and the INFO field annotations are updated.
     */
    @Argument(shortName="select", doc="One or more criteria to use when selecting the data", optional=true)
    private ArrayList<String> selectExpressions = new ArrayList<>();

    /**
     * Invert the selection criteria for -select.
     */
    @Argument(shortName="invert-select", doc="Invert the selection criteria for -select", optional=true)
    private boolean invertSelect = false;

    /*
     * If this flag is enabled, sites that are found to be non-variant after the subsetting procedure (i.e. where none
     * of the selected samples display evidence of variation) will be excluded from the output.
     */
    @Argument(fullName="exclude-non-variants", doc="Don't include non-variant sites", optional=true)
    private boolean XLnonVariants = false;

    /**
     * If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS`
     * in the FILTER field) will be excluded from the output.
     */
    @Argument(fullName="exclude-filtered", doc="Don't include filtered sites", optional=true)
    private boolean XLfiltered = false;

    /**
     * The default behavior of this tool is to remove bases common to all remaining alleles after subsetting
     * operations have been completed, leaving only their minimal representation. If this flag is enabled, the original
     * alleles will be preserved as recorded in the input VCF.
     */
    @Argument(fullName="preserve-alleles", doc="Preserve original alleles, do not trim", optional=true)
    private boolean preserveAlleles = false;

    /**
     * When this flag is enabled, all alternate alleles that are not present in the (output) samples will be removed.
     * Note that this even extends to biallelic SNPs - if the alternate allele is not present in any sample, it will be
     * removed and the record will contain a '.' in the ALT column. Note also that sites-only VCFs, by definition, do
     * not include the alternate allele in any genotype calls.
     */
    @Argument(fullName="remove-unused-alternates",
                    doc="Remove alternate alleles not present in any genotypes", optional=true)
    private boolean removeUnusedAlternates = false;

    /**
     * When this argument is used, we can choose to include only multiallelic or biallelic sites, depending on how many
     * alleles are listed in the ALT column of a VCF. For example, a multiallelic record such as:
     *     1    100 .   A   AAA,AAAAA
     * will be excluded if `-restrict-alleles-to BIALLELIC` is used, because there are two alternate alleles, whereas a
     * record such as:
     *     1    100 .   A  T
     * will be included in that case, but would be excluded if `-restrict-alleles-to MULTIALLELIC` is used.
     * Valid options are ALL (default), MULTIALLELIC or BIALLELIC.
     */
    @Argument(fullName="restrict-alleles-to",
                    doc="Select only variants of a particular allelicity", optional=true)
    private NumberAlleleRestriction alleleRestriction = NumberAlleleRestriction.ALL;

    /**
     * When subsetting a callset, this tool recalculates the AC, AF, and AN values corresponding to the contents of the
     * subset. If this flag is enabled, the original values of those annotations will be stored in new annotations called
     * AC_Orig, AF_Orig, and AN_Orig.
     */
    @Argument(fullName="keep-original-ac",
                    doc="Store the original AC, AF, and AN values after subsetting", optional=true)
    private boolean keepOriginalChrCounts = false;

    /**
     * When subsetting a callset, this tool recalculates the site-level (INFO field) DP value corresponding to the
     * contents of the subset. If this flag is enabled, the original value of the DP annotation will be stored in
     * a new annotation called DP_Orig.
     */
    @Argument(fullName="keep-original-dp",
                    doc="Store the original DP value after subsetting", optional=true)
    private boolean keepOriginalDepth = false;

    /**
     * If this flag is enabled, this tool will select only variants that correspond to a mendelian violation as
     * determined on the basis of family structure. Requires passing a pedigree file using the engine-level
     * `-ped` argument.
     */
    @Argument(fullName="mendelian-violation", doc="Output mendelian violation sites only", optional=true)
    private Boolean mendelianViolations = false;

    /**
     * If this flag is enabled, this tool will select only variants that do not correspond to a mendelian violation as
     * determined on the basis of family structure. Requires passing a pedigree file using the engine-level
     * `-ped` argument.
     */
    @Argument(fullName="invert-mendelian-violation",
                    doc="Output non-mendelian violation sites only", optional=true)
    private Boolean invertMendelianViolations = false;

    /**
     * This argument specifies the genotype quality (GQ) threshold that all members of a trio must have in order
     * for a site to be accepted as a mendelian violation. Note that the `-mv` flag must be set for this argument
     * to have an effect.
     */
    @Argument(fullName="mendelian-violation-qual-threshold",
                    doc="Minimum GQ score for each trio member to accept a site as a violation", optional=true)
    private double mendelianViolationQualThreshold = 0;

    @Argument(fullName=StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName=StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file", optional=true)
    private GATKPath pedigreeFile = null;

    /**
     * The value of this argument should be a number between 0 and 1 specifying the fraction of total variants to be
     * randomly selected from the input callset. Note that this is done using a probabilistic function, so the final
     * result is not guaranteed to carry the exact fraction requested. Can be used for large fractions.
     */
    @Argument(fullName="select-random-fraction", shortName="fraction",
                    doc="Select a fraction of variants at random from the input", optional=true)
    private double fractionRandom = 0;

    /**
     * The value of this argument should be a number between 0 and 1 specifying the fraction of total variants to be
     * randomly selected from the input callset and set to no-call (./). Note that this is done using a probabilistic
     * function, so the final result is not guaranteed to carry the exact fraction requested. Can be used for large fractions.
     */
    @Argument(fullName="remove-fraction-genotypes",
                        doc="Select a fraction of genotypes at random from the input and sets them to no-call", optional=true)
    private double fractionGenotypes = 0;

    /**
     * This argument selects particular kinds of variants out of a list. If left empty, there is no type selection
     * and all variant types are considered for other selection criteria. Valid types are INDEL, SNP, MIXED, MNP,
     * SYMBOLIC, NO_VARIATION. Can be specified multiple times.
     */
    @Argument(fullName="select-type-to-include", shortName="select-type",
                    doc="Select only a certain type of variants from the input file", optional=true)
    private List<VariantContext.Type> typesToInclude = new ArrayList<>();

    /**
     * This argument excludes particular kinds of variants out of a list. If left empty, there is no type selection
     * and all variant types are considered for other selection criteria. Valid types are INDEL, SNP, MIXED, MNP,
     * SYMBOLIC, NO_VARIATION. Can be specified multiple times.
     */
    @Argument(fullName="select-type-to-exclude", shortName="xl-select-type",
                    doc="Do not select certain type of variants from the input file", optional=true)
    private List<VariantContext.Type> typesToExclude = new ArrayList<>();

    /**
     * List of IDs (or a .list file containing ids) to select. The tool will only select variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. If a file, the file
     * name must end in ".list", and the expected file format is simply plain text with one ID per line.
     */
    @Argument(fullName="keep-ids", shortName="ids", doc="List of variant rsIDs to select", optional=true)
    private Set<String> rsIDsToKeep = new HashSet<>();

    /**
     * List of IDs (or a .list file containing ids) to exclude. The tool will exclude variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. If a file, the
     * file name must end in ".list", and the expected file format is simply plain text with one ID per line.
     */
    @Argument(fullName="exclude-ids", shortName="xl-ids", doc="List of variant rsIDs to exclude", optional=true)
    private Set<String> rsIDsToRemove = new HashSet<>();

    @Hidden
    @Argument(fullName="fully-decode", doc="If true, the incoming VariantContext will be fully decoded", optional=true)
    private boolean fullyDecode = false;

    /**
     * If this argument is provided, indels that are larger than the specified size will be excluded.
     */
    @Argument(fullName="max-indel-size", optional=true, doc="Maximum size of indels to include")
    private int maxIndelSize = Integer.MAX_VALUE;

    /**
     * If this argument is provided, indels that are smaller than the specified size will be excluded.
     */
    @Argument(fullName="min-indel-size", optional=true, doc="Minimum size of indels to include")
    private int minIndelSize = 0;

    /**
     * If this argument is provided, select sites where at most a maximum number of samples are filtered at the genotype level.
     */
    @Argument(fullName="max-filtered-genotypes", optional=true, doc="Maximum number of samples filtered at the genotype level")
    private int maxFilteredGenotypes = MAX_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where at least a minimum number of samples are filtered at
     * the genotype level.
     */
    @Argument(fullName="min-filtered-genotypes", optional=true,
                    doc="Minimum number of samples filtered at the genotype level")
    private int minFilteredGenotypes = MIN_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where a fraction or less of the samples are filtered at
     * the genotype level.
     */
    @Argument(fullName="max-fraction-filtered-genotypes",
                    optional=true, doc="Maximum fraction of samples filtered at the genotype level")
    private double maxFractionFilteredGenotypes = MAX_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where a fraction or more of the samples are filtered at
     * the genotype level.
     */
    @Argument(fullName="min-fraction-filtered-genotypes", optional=true,
                    doc="Maximum fraction of samples filtered at the genotype level")
    private double minFractionFilteredGenotypes = MIN_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where at most the given number of samples have no-call genotypes.
     */
    @Argument(fullName="max-nocall-number", optional=true,
            doc="Maximum number of samples with no-call genotypes")
    private int maxNOCALLnumber = MAX_NOCALL_NUMBER_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where at most the given fraction of samples have no-call genotypes.
     */
    @Argument(fullName="max-nocall-fraction", optional=true,
            doc="Maximum fraction of samples with no-call genotypes")
    private double maxNOCALLfraction = MAX_NOCALL_FRACTION_DEFAULT_VALUE;

    /**
     * If this argument is provided, set filtered genotypes to no-call (./.).
     */
    @Argument(fullName="set-filtered-gt-to-nocall", optional=true, doc="Set filtered genotypes to no-call")
    private boolean setFilteredGenotypesToNocall = false;

    /**
     * Info annotation fields to be dropped (specified by key)
     */
    @Argument(fullName = "drop-info-annotation", shortName = "DA", optional = true, doc = "Info annotations to drop from output vcf.  Annotations to be dropped are specified by their key.")
    private List<String> infoAnnotationsToDrop = new ArrayList<>();

    /**
     * Genotype annotation fields to be dropped (specified by key)
     */
    @Argument(fullName = "drop-genotype-annotation", shortName = "DGA", optional = true, doc = "Genotype annotations to drop from output vcf.  Annotations to be dropped are specified by their key.")
    private List<String> genotypeAnnotationsToDrop = new ArrayList<>();

    @Hidden
    @Argument(fullName="allow-nonoverlapping-command-line-samples", optional=true,
                    doc="Allow samples other than those in the VCF to be specified on the command line. These samples will be ignored.")
    private boolean allowNonOverlappingCommandLineSamples = false;

    @Hidden
    @Argument(fullName="suppress-reference-path", optional=true,
            doc="Suppress reference path in output for test result differencing")
    private boolean suppressReferencePath = false;

    @ArgumentCollection
    private GenomicsDBArgumentCollection genomicsdbArgs = new GenomicsDBArgumentCollection();

    private VariantContextWriter vcfWriter = null;

    private enum NumberAlleleRestriction {
        ALL,
        BIALLELIC,
        MULTIALLELIC
    }

    private SortedSet<String> samples = new TreeSet<>();
    private boolean noSamplesSpecified = false;

    private Set<VariantContext.Type> selectedTypes = new LinkedHashSet<>();
    private final ArrayList<String> selectNames = new ArrayList<>();
    private List<VariantContextUtils.JexlVCMatchExp> jexls = null;

    private boolean discordanceOnly = false;
    private boolean concordanceOnly = false;

    private MendelianViolation mv = null;
    private SampleDB sampleDB = null;

    /* variables used by the SELECT RANDOM modules */
    private boolean selectRandomFraction = false;

    // Random number generator for the genotypes to remove
    private final Random randomGenotypes = new Random();

    private final List<Allele> diploidNoCallAlleles = GATKVariantContextUtils.noCallAlleles(2);

    private final Map<Integer, Integer> ploidyToNumberOfAlleles = new LinkedHashMap<Integer, Integer>();

    @Override
    protected GenomicsDBOptions getGenomicsDBOptions() {
        if (genomicsDBOptions == null) {
            genomicsDBOptions = new GenomicsDBOptions(referenceArguments.getReferencePath(), genomicsdbArgs);
        }
        return genomicsDBOptions;
    }

    final private PriorityQueue<VariantContext> pendingVariants = new PriorityQueue<>(Comparator.comparingInt(VariantContext::getStart));

    /**
     * Set up the VCF writer, the sample expressions and regexs, filters inputs, and the JEXL matcher
     *
     */
    @Override
    public void onTraversalStart() {
        final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());

        // Initialize VCF header lines
        final Set<VCFHeaderLine> headerLines = createVCFHeaderLineList(vcfHeaders);

        for (int i = 0; i < selectExpressions.size(); i++) {
            // It's not necessary that the user supply select names for the JEXL expressions, since those
            // expressions will only be needed for omitting records.  Make up the select names here.
            selectNames.add(String.format("select-%d", i));
        }

        jexls = VariantContextUtils.initializeMatchExps(selectNames, selectExpressions);

        // Prepare the sample names and types to be used by the corresponding filters
        samples = createSampleNameInclusionList(vcfHeaders);
        selectedTypes = createSampleTypeInclusionList();

        // Look at the parameters to decide which analysis to perform
        discordanceOnly = discordanceTrack != null;
        if (discordanceOnly) {
            logger.info("Selecting only variants discordant with the track: " + discordanceTrack.getName());
        }

        concordanceOnly = concordanceTrack != null;
        if (concordanceOnly) {
            logger.info("Selecting only variants concordant with the track: " + concordanceTrack.getName());
        }

        if (mendelianViolations) {
            sampleDB = SampleDB.createSampleDBFromPedigree(pedigreeFile);
            mv = new MendelianViolation(mendelianViolationQualThreshold, false, true);
        }

        selectRandomFraction = fractionRandom > 0;
        if (selectRandomFraction) {
            logger.info("Selecting approximately " + 100.0*fractionRandom + "% of the variants at random from the variant track");
        }

        //TODO: this should be refactored/consolidated as part of
        // https://github.com/broadinstitute/gatk/issues/121 and
        // https://github.com/broadinstitute/gatk/issues/1116
        Set<VCFHeaderLine> actualLines = null;
        SAMSequenceDictionary sequenceDictionary = null;
        if (hasReference()) {
            Path refPath = referenceArguments.getReferencePath();
            sequenceDictionary= this.getReferenceDictionary();
            actualLines = VcfUtils.updateHeaderContigLines(headerLines, refPath, sequenceDictionary, suppressReferencePath);
        }
        else {
            sequenceDictionary = getHeaderForVariants().getSequenceDictionary();
            if (null != sequenceDictionary) {
                actualLines = VcfUtils.updateHeaderContigLines(headerLines, null, sequenceDictionary, suppressReferencePath);
            }
            else {
                actualLines = headerLines;
            }
        }
        if (!infoAnnotationsToDrop.isEmpty()) {
            for (final String infoField : infoAnnotationsToDrop) {
                logger.info(String.format("Will drop info annotation: %s",infoField));
            }
        }
        if (!genotypeAnnotationsToDrop.isEmpty()) {
            for (final String genotypeAnnotation : genotypeAnnotationsToDrop) {
                logger.info(String.format("Will drop genotype annotation: %s",genotypeAnnotation));
            }
        }

        final Path outPath = vcfOutput.toPath();
        vcfWriter = createVCFWriter(outPath);
        vcfWriter.writeHeader(new VCFHeader(actualLines, samples));
    }

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {

        /*check for pending variants to write out
        since variant starts will only be moved further right, we can write out a pending variant if the current variant start is after the pending variant start
        variant record locations can move to the right due to allele trimming if preserveAlleles is false
         */
        while (!pendingVariants.isEmpty() && (pendingVariants.peek().getStart()<=vc.getStart() || !(pendingVariants.peek().getContig().equals(vc.getContig())))) {
            vcfWriter.add(pendingVariants.poll());
        }

        if (fullyDecode) {
            vc = vc.fullyDecode(getHeaderForVariants(), lenientVCFProcessing);
        }

        if (mendelianViolations && invertLogic((mv.countFamilyViolations(sampleDB, samples, vc) == 0), invertMendelianViolations)) {
            return;
        }

        if (discordanceOnly && !isDiscordant(vc, featureContext.getValues(discordanceTrack))) {
            return;
        }

        if (concordanceOnly && !isConcordant(vc, featureContext.getValues(concordanceTrack))) {
            return;
        }

        if (alleleRestriction.equals(NumberAlleleRestriction.BIALLELIC) && !vc.isBiallelic()) {
            return;
        }

        if (alleleRestriction.equals(NumberAlleleRestriction.MULTIALLELIC) && vc.isBiallelic()) {
            return;
        }

        if (containsIndelLargerOrSmallerThan(vc, maxIndelSize, minIndelSize)) {
            return;
        }

        if (considerFilteredGenotypes()) {
            final int numFilteredSamples = numFilteredGenotypes(vc);
            final double fractionFilteredGenotypes = samples.isEmpty() ? 0.0 : numFilteredSamples / samples.size();
            if (numFilteredSamples > maxFilteredGenotypes || numFilteredSamples < minFilteredGenotypes ||
                    fractionFilteredGenotypes > maxFractionFilteredGenotypes || fractionFilteredGenotypes < minFractionFilteredGenotypes)
                return;
        }

        if (considerNoCallGenotypes()) {
            final int numNoCallSamples = numNoCallGenotypes(vc);
            final double fractionNoCallGenotypes = samples.isEmpty() ? 0.0 : ((double) numNoCallSamples) / samples.size();
            if (numNoCallSamples > maxNOCALLnumber || fractionNoCallGenotypes > maxNOCALLfraction)
                return;
        }

        // Initialize the cache of PL index to a list of alleles for each ploidy.
        initalizeAlleleAnyploidIndicesCache(vc);

        final VariantContext sub = subsetRecord(vc, preserveAlleles, removeUnusedAlternates);
        final VariantContext filteredGenotypeToNocall;

        if ( setFilteredGenotypesToNocall ) {
            final VariantContextBuilder builder = new VariantContextBuilder(sub);
            GATKVariantContextUtils.setFilteredGenotypeToNocall(builder, sub, setFilteredGenotypesToNocall, this::getGenotypeFilters);
            filteredGenotypeToNocall = builder.make();
        } else {
            filteredGenotypeToNocall = sub;
        }

        // Not excluding non-variants OR (subsetted polymorphic variants AND not spanning deletion) AND (including filtered loci OR subsetted variant) is not filtered
        // If exclude non-variants argument is not called, filtering will NOT occur.
        // If exclude non-variants is called, and a spanning deletion exists, the spanning deletion will be filtered
        // If exclude non-variants is called, it is a polymorphic variant, but not a spanning deletion, filtering will not occur
        // True iff exclude-filtered is not called or the filteredGenotypeToNocall is not already filtered

        if ((!XLnonVariants || (filteredGenotypeToNocall.isPolymorphicInSamples() && !GATKVariantContextUtils.isSpanningDeletionOnly(filteredGenotypeToNocall)))
                && (!XLfiltered || !filteredGenotypeToNocall.isFiltered()))
        {

            // Write the subsetted variant if it matches all of the expressions
            boolean failedJexlMatch = false;

            try {
                for (VariantContextUtils.JexlVCMatchExp jexl : jexls) {
                    if (invertLogic(!VariantContextUtils.match(filteredGenotypeToNocall, jexl), invertSelect)){
                        failedJexlMatch = true;
                        break;
                    }
                }
            } catch (IllegalArgumentException e) {
                // The IAE thrown by htsjdk already includes an informative error message ("Invalid JEXL
                //  expression detected...")
                throw new UserException(e.getMessage() +
                        "\nSee https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions for documentation on using JEXL in GATK", e);
            }

            if (!failedJexlMatch &&
                    (!selectRandomFraction || Utils.getRandomGenerator().nextDouble() < fractionRandom)) {
                //remove annotations being dropped and write variantcontext
                final VariantContext variantContextToWrite = buildVariantContextWithDroppedAnnotationsRemoved(filteredGenotypeToNocall);
                pendingVariants.add(variantContextToWrite);
            }
        }
    }

    /**
     * write out all remaining pending variants
     */
    @Override
    public Object onTraversalSuccess() {
        while(!pendingVariants.isEmpty()) {
            vcfWriter.add(pendingVariants.poll());
        }
        return null;
    }

    private VariantContext buildVariantContextWithDroppedAnnotationsRemoved(final VariantContext vc) {
        if (infoAnnotationsToDrop.isEmpty() && genotypeAnnotationsToDrop.isEmpty()) {
            return vc;
        }
        final VariantContextBuilder rmAnnotationsBuilder = new VariantContextBuilder(vc);
        for (String infoField : infoAnnotationsToDrop) {
            rmAnnotationsBuilder.rmAttribute(infoField);
        }
        if (!genotypeAnnotationsToDrop.isEmpty()) {
            final ArrayList<Genotype> genotypesToWrite = new ArrayList<>();
            for (Genotype genotype : vc.getGenotypes()) {
                final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype).noAttributes();
                final Map<String, Object> attributes = new HashMap<>(genotype.getExtendedAttributes());
                for (String genotypeAnnotation : genotypeAnnotationsToDrop) {
                    attributes.remove(genotypeAnnotation);
                }
                genotypeBuilder.attributes(attributes);
                genotypesToWrite.add(genotypeBuilder.make());
            }
            rmAnnotationsBuilder.genotypes(GenotypesContext.create(genotypesToWrite));
        }
        return rmAnnotationsBuilder.make();
    }

    /**
     * Get the genotype filters
     *
     * @param vc the variant context
     * @param g the genotype
     * @return list of genotype filter names
     */
    private List<String> getGenotypeFilters(final VariantContext vc, final Genotype g) {
        final List<String> filters = new ArrayList<>();
        if (g.isFiltered()) {
            filters.add(g.getFilters());
        }

        return filters;
    }

    /**
     * Initialize cache of allele anyploid indices
     *
     * Initialize the cache of PL index to a list of alleles for each ploidy.
     *
     * @param vc    Variant Context
     */
    private void initalizeAlleleAnyploidIndicesCache(final VariantContext vc) {
        if (vc.getType() != VariantContext.Type.NO_VARIATION &&
                vc.getGenotypes().stream().anyMatch(Genotype::hasLikelihoods)) { // Bypass if not a variant or no genotype has Pls
            for (final Genotype g : vc.getGenotypes()) {
                // Make a new entry if the we have not yet cached a PL to allele indices map for this ploidy and allele count
                // skip if there are no PLs -- this avoids hanging on high-allelic somatic samples, for example, where
                // there's no need for the PL indices since they don't exist
                final int genotypePloidy = g.getPloidy();
                if (genotypePloidy != 0 && (!ploidyToNumberOfAlleles.containsKey(genotypePloidy) || ploidyToNumberOfAlleles.get(genotypePloidy) < vc.getNAlleles())) {
                    GenotypeLikelihoods.initializeAnyploidPLIndexToAlleleIndices(vc.getNAlleles() - 1, genotypePloidy);
                    ploidyToNumberOfAlleles.put(genotypePloidy, vc.getNAlleles());
                }
            }
        }

    }

    /**
     * Close out the new variants file.
     */
    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    /**
     * Create filters for variant types, ids, and genomic intervals.
     */
    @Override
    protected CountingVariantFilter makeVariantFilter() {
        CountingVariantFilter compositeFilter = new CountingVariantFilter(VariantFilterLibrary.ALLOW_ALL_VARIANTS);

        if (!selectedTypes.isEmpty()) {
            compositeFilter = compositeFilter.and(new CountingVariantFilter(new VariantTypesVariantFilter(selectedTypes)));
        }

        if (rsIDsToKeep != null && !rsIDsToKeep.isEmpty()) {
            compositeFilter = compositeFilter.and(new CountingVariantFilter(new VariantIDsVariantFilter(rsIDsToKeep)));
        }

        if (rsIDsToRemove != null && !rsIDsToRemove.isEmpty()) {
            compositeFilter = compositeFilter.and(new CountingVariantFilter(new VariantIDsVariantFilter(rsIDsToRemove).negate()));
        }

        return compositeFilter;
    }

    /**
     * Prepare the sample names to be included(/excluded) in the output by the names filter.
     */
    private SortedSet<String> createSampleNameInclusionList(Map<String, VCFHeader> vcfHeaders) {
        final SortedSet<String> vcfSamples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
        final Collection<String> samplesFromExpressions = Utils.filterCollectionByExpressions(vcfSamples, sampleExpressions, false);

        // first, find any samples that were listed on the command line but which don't exist in in the header
        final Set<String> samplesNotInHeader = new LinkedHashSet<>(samplesFromExpressions.size()+sampleNames.size());
        samplesNotInHeader.addAll(samplesFromExpressions);
        samplesNotInHeader.addAll(sampleNames);
        samplesNotInHeader.removeAll(vcfSamples);

        // second, add the requested samples
        samples.addAll(sampleNames);
        samples.addAll(samplesFromExpressions);

        // report any requested samples that don't exist in the header and remove them from the list we're accumulating
        logger.debug(Utils.join(",", samplesNotInHeader));
        if (!samplesNotInHeader.isEmpty()) {
            if (allowNonOverlappingCommandLineSamples) {
                logger.warn("Samples present on command line input that are not present in the VCF. These samples will be ignored.");
                samples.removeAll(samplesNotInHeader);
            }
            else {
                throw new UserException.BadInput(String.format("%s%n%n%s%n%n%s%n%n%s",
                        "Samples entered on command line (through -sf or -sn) that are not present in the VCF.",
                        "A list of these samples:",
                        Utils.join(",", samplesNotInHeader),
                        "To ignore these samples, run with --allow-nonoverlapping-command-line-samples"));
            }
        }

        // if none were requested, we want all of them
        if (samples.isEmpty()) {
            samples.addAll(vcfSamples);
            noSamplesSpecified = true;
        }

        // Exclude samples take precedence over include - remove any excluded samples
        final Collection<String> XLsamplesFromExpressions = Utils.filterCollectionByExpressions(vcfSamples, XLsampleExpressions, false);
        samples.removeAll(XLsampleNames);
        samples.removeAll(XLsamplesFromExpressions);
        noSamplesSpecified = noSamplesSpecified &&
                                XLsampleNames.isEmpty() &&
                                XLsamplesFromExpressions.isEmpty();

        if (samples.isEmpty() && !noSamplesSpecified) {
            throw new UserException("All samples requested to be included were also requested to be excluded.");
        }

        if (!noSamplesSpecified) {
            for (String sample : samples) {
                logger.info("Including sample '" + sample + "'");
            }
        }

        return samples;
    }

    /**
     * Prepare the type inclusion list to be used by the type filter
     */
    private Set<VariantContext.Type> createSampleTypeInclusionList() {

        // if user specified types to include, add these, otherwise, add all possible variant context types to list of vc types to include
        if (typesToInclude.isEmpty()) {
            selectedTypes.addAll(Arrays.asList(VariantContext.Type.values()));
        }
        else {
            selectedTypes.addAll(typesToInclude);
        }

        // Exclude types take precedence over include - remove specified exclude types
        selectedTypes.removeAll(typesToExclude);

        return selectedTypes;
    }

    /**
     * Prepare the VCF header lines
     */
    private Set<VCFHeaderLine> createVCFHeaderLineList(Map<String, VCFHeader> vcfHeaders) {

        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfHeaders.values(), true);
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // need AC, AN and AF since output if set filtered genotypes to no-call
        if (setFilteredGenotypesToNocall) {
            GATKVariantContextUtils.addChromosomeCountsToHeader(headerLines);
        }

        if (keepOriginalChrCounts) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AC_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AF_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AN_KEY));
        }
        if (keepOriginalDepth) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_DP_KEY));
        }

        headerLines.addAll(Arrays.asList(ChromosomeCounts.descriptions));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));

        //remove header lines for info field and genotype annotations being dropped
        headerLines.removeIf(l->l instanceof VCFInfoHeaderLine && infoAnnotationsToDrop.contains(((VCFInfoHeaderLine)l).getID()));
        headerLines.removeIf(l->l instanceof VCFFormatHeaderLine && genotypeAnnotationsToDrop.contains(((VCFFormatHeaderLine)l).getID()));

        return headerLines;
    }

    /**
     * Invert logic if specified
     *
     * @param logic boolean logical operation value
     * @param invert whether to invert logic
     * @return invert logic if invert flag is true, otherwise leave the logic
     */
    private static boolean invertLogic(final boolean logic, final boolean invert){
        return invert ? !logic : logic;
    }

    /*
     * Determines if any of the alternate alleles are greater than the max indel size or less than the min indel size
     *
     * @param vc            the variant context to check
     * @param max-indel-size  the maximum size of allowed indels
     * @param min-indel-size  the minimum size of allowed indels
     * @return true if the VC contains an indel larger than max-indel-size or less than the min-indel-size, false otherwise
     *
     * Protected for unit test access
     */
    protected static boolean containsIndelLargerOrSmallerThan(final VariantContext vc, final int maxIndelSize, final int minIndelSize) {
        final List<Integer> lengths = vc.getIndelLengths();
        if (lengths == null)
            return false;

        for (final Integer indelLength : lengths) {
            if (Math.abs(indelLength) > maxIndelSize || Math.abs(indelLength) < minIndelSize)
                return true;
        }

        return false;
    }

    /**
     * Find the number of no-call genotypes
     *
     * @param vc the variant context
     * @return number of filtered samples
     */
    private int numNoCallGenotypes(final VariantContext vc) {
        return numGenotypes(vc, Genotype::isNoCall);
    }

    /**
     * Find the number of filtered samples
     *
     * @param vc the variant context
     * @return number of filtered samples
     */
    private int numFilteredGenotypes(final VariantContext vc)  {
        return numGenotypes(vc, g -> g.isFiltered() && !g.getFilters().isEmpty());
    }

    /**
     * Find the number of samples passing the given filter.
     *
     * @param vc the variant
     * @param f predicate by which to filter genotypes
     * @return number of filtered samples
     */
    private int numGenotypes(final VariantContext vc, final Predicate<Genotype> f)  {
        return vc == null ? 0 : (int)vc.getGenotypes(samples).stream().filter(f).count();
    }

    /**
     * Checks if vc has a variant call for (at least one of) the samples.
     *
     * @param vc the variant rod VariantContext. Here, the variant is the dataset you're looking for discordances to (e.g. HapMap)
     * @param compVCs the comparison VariantContext (discordance)
     * @return true VariantContexts are discordant, false otherwise
     */
    private boolean isDiscordant (final VariantContext vc, final Collection<VariantContext> compVCs) {
        if (vc == null) {
            return false;
        }

        // if we're not looking at specific samples then the absence of a compVC means discordance
        if (noSamplesSpecified) {
            return (compVCs == null || compVCs.isEmpty());
        }

        // check if we find it in the variant rod
        final GenotypesContext genotypes = vc.getGenotypes(samples);
        for (final Genotype g : genotypes) {
            if (sampleHasVariant(g)) {
                // There is a variant called (or filtered with not exclude filtered option set) that is not HomRef for at least one of the samples.
                if (compVCs == null) {
                    return true;
                }
                // Look for this sample in the all vcs of the comp ROD track.
                boolean foundVariant = false;
                for (final VariantContext compVC : compVCs) {
                    if (haveSameGenotypes(g, compVC.getGenotype(g.getSampleName()))) {
                        foundVariant = true;
                        break;
                    }
                }
                // if (at least one sample) was not found in all VCs of the comp ROD, we have discordance
                if (!foundVariant)
                    return true;
            }
        }
        return false; // we only get here if all samples have a variant in the comp rod.
    }

    /**
     * Checks if the two variants have the same genotypes for the selected samples
     *
     * @param vc the variant rod VariantContext.
     * @param compVCs the comparison VariantContext
     * @return true if VariantContexts are concordant, false otherwise
     */
    private boolean isConcordant (final VariantContext vc, final Collection<VariantContext> compVCs) {
        if (vc == null || compVCs == null || compVCs.isEmpty()) {
            return false;
        }

        // if we're not looking for specific samples then the fact that we have both VCs is enough to call it concordant.
        if (noSamplesSpecified) {
            return true;
        }

        // make a list of all samples contained in this variant VC that are being tracked by the user command line arguments.
        final Set<String> variantSamples = vc.getSampleNames();
        variantSamples.retainAll(samples);

        // check if we can find all samples from the variant rod in the comp rod.
        for (final String sample : variantSamples) {
            boolean foundSample = false;
            for (final VariantContext compVC : compVCs) {
                final Genotype varG = vc.getGenotype(sample);
                final Genotype compG = compVC.getGenotype(sample);
                if (haveSameGenotypes(varG, compG)) {
                    foundSample = true;
                    break;
                }
            }
            // if at least one sample doesn't have the same genotype, we don't have concordance
            if (!foundSample) {
                return false;
            }
        }
        return true;
    }

    private boolean sampleHasVariant(final Genotype g) {
        return (g !=null && !g.isHomRef() && (g.isCalled() || (g.isFiltered() && !XLfiltered)));
    }

    private boolean haveSameGenotypes(final Genotype g1, final Genotype g2) {
        if (g1 == null || g2 == null) {
            return false;
        }

        if ((g1.isCalled() && g2.isFiltered()) ||
                (g2.isCalled() && g1.isFiltered()) ||
                (g1.isFiltered() && g2.isFiltered() && XLfiltered)) {
            return false;
        }

        final List<Allele> a1s = g1.getAlleles();
        final List<Allele> a2s = g2.getAlleles();
        return (a1s.containsAll(a2s) && a2s.containsAll(a1s));
    }

    /**
     * Helper method to subset a VC record, modifying some metadata stored in the INFO field (i.e. AN, AC, AF).
     *
     * @param vc       the VariantContext record to subset
     * @param preserveAlleles should we trim constant sequence from the beginning and/or end of all alleles, or preserve it?
     * @param removeUnusedAlternates removes alternate alleles with AC=0
     * @return the subsetted VariantContext
     */
    private VariantContext subsetRecord(final VariantContext vc, final boolean preserveAlleles, final boolean removeUnusedAlternates) {
        //subContextFromSamples() always decodes the vc, which is a fairly expensive operation.  Avoid if possible
        if (noSamplesSpecified && !removeUnusedAlternates) {
            return vc;
        }

        // strip out the alternate alleles that aren't being used
        final VariantContext sub = vc.subContextFromSamples(samples, removeUnusedAlternates);

        // If no subsetting happened, exit now
        if (sub.getNSamples() == vc.getNSamples() && sub.getNAlleles() == vc.getNAlleles()) {
            return vc;
        }

        // fix the PL and AD values if sub has fewer alleles than original vc and remove a fraction of the genotypes if needed
        final GenotypesContext oldGs = sub.getGenotypes();
        GenotypesContext newGC = sub.getNAlleles() == vc.getNAlleles() ? oldGs :
                AlleleSubsettingUtils.subsetAlleles(oldGs, 0, vc.getAlleles(), sub.getAlleles(), null, GenotypeAssignmentMethod.DO_NOT_ASSIGN_GENOTYPES, vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));

        if (fractionGenotypes > 0) {
            final List<Genotype> genotypes = newGC.stream().map(genotype -> randomGenotypes.nextDouble() > fractionGenotypes ? genotype :
                    new GenotypeBuilder(genotype).alleles(getNoCallAlleles(genotype.getPloidy())).noGQ().make()).collect(Collectors.toList());
            newGC = GenotypesContext.create(new ArrayList<>(genotypes));
        }

        // since the VC has been subset (either by sample or allele), we need to strip out the MLE tags
        final VariantContextBuilder builder = new VariantContextBuilder(sub);
        builder.rmAttributes(Arrays.asList(GATKVCFConstants.MLE_ALLELE_COUNT_KEY,GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        builder.genotypes(newGC);
        addAnnotations(builder, vc, sub.getSampleNames());
        final VariantContext subset = builder.make();

        return preserveAlleles? subset : GATKVariantContextUtils.trimAlleles(subset,true,true);
    }

    /**
     * Get the ploidy number of NO-CALL alleles
     *
     * @param ploidy    number of sets of chromosomes
     * @return  the NO-CALL alleles
     */
    private List<Allele> getNoCallAlleles(final int ploidy) {
        return ploidy == 2 ? diploidNoCallAlleles : GATKVariantContextUtils.noCallAlleles(ploidy);
    }

    /*
     * Add annotations to the new VC
     *
     * @param builder     the new VC to annotate
     * @param originalVC  the original VC
     * @param selectedSampleNames the post-selection list of sample names
     */
    private void addAnnotations(final VariantContextBuilder builder, final VariantContext originalVC, final Set<String> selectedSampleNames) {
        if (fullyDecode) {
            return; // TODO -- annotations are broken with fully decoded data
        }

        if (keepOriginalChrCounts) {
            final int[] indexOfOriginalAlleleForNewAllele;
            final List<Allele> newAlleles = builder.getAlleles();
            final int numOriginalAlleles = originalVC.getNAlleles();

            // if the alleles already match up, we can just copy the previous list of counts
            if (numOriginalAlleles == newAlleles.size()) {
                indexOfOriginalAlleleForNewAllele = null;
            }
            // otherwise we need to parse them and select out the correct ones
            else {
                indexOfOriginalAlleleForNewAllele = new int[newAlleles.size() - 1];
                Arrays.fill(indexOfOriginalAlleleForNewAllele, -1);

                // note that we don't care about the reference allele at position 0
                for (int newI = 1; newI < newAlleles.size(); newI++) {
                    final Allele newAlt = newAlleles.get(newI);
                    for (int oldI = 0; oldI < numOriginalAlleles - 1; oldI++) {
                        if (newAlt.equals(originalVC.getAlternateAllele(oldI), false)) {
                            indexOfOriginalAlleleForNewAllele[newI - 1] = oldI;
                            break;
                        }
                    }
                }
            }

            if (originalVC.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                builder.attribute(GATKVCFConstants.ORIGINAL_AC_KEY,
                        getReorderedAttributes(originalVC.getAttribute(VCFConstants.ALLELE_COUNT_KEY), indexOfOriginalAlleleForNewAllele));
            }
            if (originalVC.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                builder.attribute(GATKVCFConstants.ORIGINAL_AF_KEY,
                        getReorderedAttributes(originalVC.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), indexOfOriginalAlleleForNewAllele));
            }
            if (originalVC.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
                builder.attribute(GATKVCFConstants.ORIGINAL_AN_KEY, originalVC.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
            }
        }

        VariantContextUtils.calculateChromosomeCounts(builder, false);

        if (keepOriginalDepth && originalVC.hasAttribute(VCFConstants.DEPTH_KEY)) {
            builder.attribute(GATKVCFConstants.ORIGINAL_DP_KEY, originalVC.getAttribute(VCFConstants.DEPTH_KEY));
        }

        boolean sawDP = false;
        int depth = 0;
        for (final String sample : selectedSampleNames ) {
            final Genotype g = originalVC.getGenotype(sample);
            if (!g.isFiltered()) {
                if (g.hasDP()) {
                    depth += g.getDP();
                    sawDP = true;
                }
            }
        }

        if (sawDP) {
            builder.attribute(VCFConstants.DEPTH_KEY, depth);
        }
    }

    /**
     * Pulls out the appropriate tokens from the old ordering of an attribute to the new ordering
     *
     * @param attribute               the non-null attribute (from the INFO field)
     * @param oldToNewIndexOrdering   the mapping from new to old ordering
     * @return non-null Object attribute
     */
    private Object getReorderedAttributes(final Object attribute, final int[] oldToNewIndexOrdering) {
        // if the ordering is the same, then just use the original attribute
        if (oldToNewIndexOrdering == null) {
            return attribute;
        }

        // break the original attributes into separate tokens; unfortunately, this means being smart about class types
        final Object[] tokens;
        if (attribute.getClass().isArray()) {
            tokens = (Object[]) attribute;
        } else if (List.class.isAssignableFrom(attribute.getClass())) {
            tokens = ((List) attribute).toArray();
        } else {
            tokens = attribute.toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
        }

        Utils.validateArg(Arrays.stream(oldToNewIndexOrdering).allMatch(index -> index < tokens.length), () ->
                "the old attribute has an incorrect number of elements: " + attribute);
        return Arrays.stream(oldToNewIndexOrdering).mapToObj(index -> tokens[index]).collect(Collectors.toList());
    }

    /**
     * Should the number of filtered genotypes be considered for filtering?
     *
     * @return true if any of the filtered genotype samples arguments is used (not the default value), false otherwise
     */
    private boolean considerFilteredGenotypes(){
        return maxFilteredGenotypes != MAX_FILTERED_GENOTYPES_DEFAULT_VALUE ||
                minFilteredGenotypes != MIN_FILTERED_GENOTYPES_DEFAULT_VALUE ||
                maxFractionFilteredGenotypes != MAX_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE ||
                minFractionFilteredGenotypes != MIN_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE;
    }

    /**
     * Should the number of no-call genotypes be considered for filtering?
     *
     * @return true if any of the filtered genotype samples arguments is used (not the default value), false otherwise
     */
    private boolean considerNoCallGenotypes(){
        return maxNOCALLnumber != MAX_NOCALL_NUMBER_DEFAULT_VALUE ||
                maxNOCALLfraction != MAX_NOCALL_FRACTION_DEFAULT_VALUE;
    }
}
