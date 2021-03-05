package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantLocusWalker;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBArgumentCollection;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
 *
 * <p>
 * This tool is designed to perform joint genotyping on a single input, which may contain one or many samples. In any
 * case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or
 * `-ERC BP_RESOLUTION`.
 *
 *
 * <h3>Input</h3>
 * <p>
 * The GATK4 GenotypeGVCFs tool can take only one input track.  Options are 1) a single single-sample GVCF 2) a single
 * multi-sample GVCF created by CombineGVCFs or 3) a GenomicsDB workspace created by GenomicsDBImport.
 * A sample-level GVCF is produced by HaplotypeCaller with the `-ERC GVCF` setting.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A final VCF in which all samples have been jointly genotyped.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Perform joint genotyping on a singular sample by providing a single-sample GVCF or on a cohort by providing a combined multi-sample GVCF</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" GenotypeGVCFs \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V input.g.vcf.gz \
 *   -O output.vcf.gz
 * </pre>
 *
 * <h4>Perform joint genotyping on GenomicsDB workspace created with GenomicsDBImport</h4>
 * <pre>
 * gatk --java-options "-Xmx4g" GenotypeGVCFs \
 *   -R Homo_sapiens_assembly38.fasta \
 *   -V gendb://my_database \
 *   -O output.vcf.gz \
 *   --tmp-dir /path/to/large/tmp
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 *   <li>Only GVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call GVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</li>
 *   <li>Cannot take multiple GVCF files in one command.</li>
 *   <li>The amount of temporary disk storage required by GenomicsDBImport may exceed what is available in the default location: `/tmp`. The command line argument `--tmp-dir` can be used to specify an alternate temperary storage location with sufficient space.</li>
 * </ul>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
@CommandLineProgramProperties(summary = "Perform joint genotyping on a single-sample GVCF from HaplotypeCaller or a multi-sample GVCF from CombineGVCFs or GenomicsDBImport",
        oneLineSummary = "Perform joint genotyping on one or more samples pre-called with HaplotypeCaller",
        programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
public final class GenotypeGVCFs extends VariantLocusWalker {

    public static final String PHASED_HOM_VAR_STRING = "1|1";
    public static final String ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME = "only-output-calls-starting-in-intervals";
    public static final String ALL_SITES_LONG_NAME = "include-non-variant-sites";
    public static final String ALL_SITES_SHORT_NAME = "all-sites";
    public static final String KEEP_COMBINED_LONG_NAME = "keep-combined-raw-annotations";
    public static final String KEEP_COMBINED_SHORT_NAME = "keep-combined";
    public static final String FORCE_OUTPUT_INTERVALS_NAME = "force-output-intervals";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private GATKPath outputFile;

    @Argument(fullName=ALL_SITES_LONG_NAME, shortName=ALL_SITES_SHORT_NAME, doc="Include loci found to be non-variant after genotyping", optional=true)
    private boolean includeNonVariants = false;

    /**
     * Import all data between specified intervals.   Improves performance using large lists of intervals, as in exome
     * sequencing, especially if GVCF data only exists for specified intervals.  Use with
     * --only-output-calls-starting-in-intervals if input GVCFs contain calls outside the specified intervals.
     */
    @Argument(fullName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            shortName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            doc = "Boolean flag to import all data in between intervals.")
    private boolean mergeInputIntervals = false;

    /**
     * "Genotype" somatic GVCFs, outputting genotypes according to confidently called alt alleles, which may lead to inconsistent ploidy
     * Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.
     */
    @Argument(fullName= CombineGVCFs.SOMATIC_INPUT_LONG_NAME, doc = "Finalize input GVCF according to somatic (i.e. Mutect2) TLODs (BETA feature)")
    protected boolean somaticInput = false;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     */
    @Argument(fullName=M2ArgumentCollection.EMISSION_LOD_LONG_NAME, shortName = M2ArgumentCollection.EMISSION_LOG_SHORT_NAME,
    doc = "LOD threshold to emit variant to VCF.")
    protected double tlodThreshold = 3.5;  //allow for some lower quality variants


    /**
     * Margin of error in allele fraction to consider a somatic variant homoplasmic, i.e. if there is less than a 0.1% reference allele fraction, those reads are likely errors
     */
    @Argument(fullName=CombineGVCFs.ALLELE_FRACTION_DELTA_LONG_NAME, doc = "Margin of error in allele fraction to consider a somatic variant homoplasmic")
    protected double afTolerance = 1e-3;  //based on Q30 as a "good" base quality score

    /**
     * If specified, keep the combined raw annotations (e.g. AS_SB_TABLE) after genotyping.  This is applicable to Allele-Specific annotations
     */
    @Argument(fullName=KEEP_COMBINED_LONG_NAME, shortName = KEEP_COMBINED_SHORT_NAME, doc = "If specified, keep the combined raw annotations")
    protected boolean keepCombined = false;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    @ArgumentCollection
    private GenomicsDBArgumentCollection genomicsdbArgs = new GenomicsDBArgumentCollection();

    /**
     * This option can only be activated if intervals are specified.
     */
    @Advanced
    @Argument(fullName= ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME,
            doc="Restrict variant output to sites that start within provided intervals",
            optional=true)
    private boolean onlyOutputCallsStartingInIntervals = false;


    @Argument(fullName = FORCE_OUTPUT_INTERVALS_NAME,
            suppressFileExpansion = true, doc = "sites at which to output genotypes even if non-variant in samples", optional = true)
    protected final List<String> forceOutputIntervalStrings = new ArrayList<>();

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set
     * when appropriate. Note that dbSNP is not used in any way for the genotyping calculations themselves.
     */
    @ArgumentCollection
    private final DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private ReferenceConfidenceVariantContextMerger merger;

    private VariantContextWriter vcfWriter;

    /** these are used when {@link #onlyOutputCallsStartingInIntervals) is true */
    private List<SimpleInterval> intervals;

    private OverlapDetector<GenomeLoc> forceOutputIntervals;

    private boolean forceOutputIntervalsPresent;

    private GenotypeGVCFsEngine gvcfEngine;

    /**
     * Get the largest interval per contig that contains the intervals specified on the command line.
     * @param getIntervals intervals to be transformed
     * @param sequenceDictionary used to validate intervals
     * @return a list of one interval per contig spanning the input intervals after processing and validation
     */
    @Override
    protected List<SimpleInterval> transformTraversalIntervals(final List<SimpleInterval> getIntervals, final SAMSequenceDictionary sequenceDictionary) {
        if (mergeInputIntervals) {
            return IntervalUtils.getSpanningIntervals(getIntervals, sequenceDictionary);
        } else {
            return getIntervals;
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected GenomicsDBOptions getGenomicsDBOptions() {
        //extract called genotypes so hom refs with no PLs aren't ambiguous
        genomicsdbArgs.callGenotypes = true;
        if (genomicsDBOptions == null) {
            genomicsDBOptions = new GenomicsDBOptions(referenceArguments.getReferencePath(), genomicsdbArgs, genotypeArgs);
        }
        return genomicsDBOptions;
    }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class);
    }

    @Override
    public void onTraversalStart() {

        if (somaticInput) {
            logger.warn("Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.");
        }

        forceOutputIntervalsPresent = !forceOutputIntervalStrings.isEmpty();

        if (includeNonVariants && forceOutputIntervalsPresent ) {
            throw new CommandLineException.BadArgumentValue(String.format("Force output (--%s) is incompatible with including " +
                    "non-variants (--%s and --%s).  Use the latter to force genotyping at all sites and the former to force genotyping only at given sites." +
                    "In both cases, variant sites are genotyped as usual.", FORCE_OUTPUT_INTERVALS_NAME, ALL_SITES_LONG_NAME, ALL_SITES_SHORT_NAME));
        }

        final GenomeLocSortedSet forceOutputLocs = IntervalUtils.loadIntervals(forceOutputIntervalStrings, IntervalSetRule.UNION,
                 IntervalMergingRule.ALL, 0, new GenomeLocParser(getSequenceDictionaryForDrivingVariants()));

        forceOutputIntervals = OverlapDetector.create(forceOutputLocs.toList());

        if (!(includeNonVariants || forceOutputIntervalsPresent)) {
            changeTraversalModeToByVariant();
        }

        final VCFHeader inputVCFHeader = getHeaderForVariants();

        if(onlyOutputCallsStartingInIntervals) {
            if( !hasUserSuppliedIntervals()) {
                throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required if --" + ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME + " was specified.");
            }
        }

        intervals = hasUserSuppliedIntervals() ? intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()) :
                Collections.emptyList();

        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, keepCombined);

        merger = new ReferenceConfidenceVariantContextMerger(annotationEngine, getHeaderForVariants(), somaticInput, false, true);

        //methods that cannot be called in engine bc its protected
        Set<VCFHeaderLine> defaultToolVCFHeaderLines = getDefaultToolVCFHeaderLines();
        vcfWriter = createVCFWriter(outputFile);

        //create engine object
        gvcfEngine = new GenotypeGVCFsEngine(annotationEngine, genotypeArgs, includeNonVariants, inputVCFHeader);

        //call initialize method in engine class that creates VCFWriter object and writes a header to it
        vcfWriter = gvcfEngine.setupVCFWriter(defaultToolVCFHeaderLines, keepCombined, dbsnp, vcfWriter);

    }

    @Override
    public void apply(final Locatable loc, List<VariantContext> variants, ReadsContext reads, ReferenceContext ref, FeatureContext features) {

        final boolean inForceOutputIntervals = forceOutputIntervalsPresent && forceOutputIntervals.overlapsAny(loc);
        final boolean forceOutput = includeNonVariants || inForceOutputIntervals;
        final VariantContext regenotypedVC = gvcfEngine.callRegion(loc, variants, ref, features, merger, somaticInput, tlodThreshold, afTolerance, forceOutput);

        if (regenotypedVC != null) {
            final SimpleInterval variantStart = new SimpleInterval(regenotypedVC.getContig(), regenotypedVC.getStart(), regenotypedVC.getStart());
            if ((forceOutput || !GATKVariantContextUtils.isSpanningDeletionOnly(regenotypedVC)) &&
                    (!onlyOutputCallsStartingInIntervals || intervals.stream().anyMatch(interval -> interval.contains (variantStart)))) {
                vcfWriter.add(regenotypedVC);
            }
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
