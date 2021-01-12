package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Set of arguments for the {@link HaplotypeCallerEngine}
 */
public class HaplotypeCallerArgumentCollection extends AssemblyBasedCallerArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    public static final String GQ_BAND_LONG_NAME = "gvcf-gq-bands";
    public static final String GQ_BAND_SHORT_NAME = "GQB";
    public static final String DO_NOT_CORRECT_OVERLAPPING_BASE_QUALITIES_LONG_NAME = "do-not-correct-overlapping-quality";
    public static final String OUTPUT_BLOCK_LOWER_BOUNDS = "floor-blocks";
    public static final String DRAGEN_GATK_MODE_LONG_NAME = "dragen-mode";


    @ArgumentCollection
    public StandardCallerArgumentCollection standardArgs = new StandardCallerArgumentCollection();

    @Override
    protected int getDefaultMaxMnpDistance() { return 0; }

    @Override
    protected ReadThreadingAssemblerArgumentCollection getReadThreadingAssemblerArgumentCollection() {
        return new HaplotypeCallerReadThreadingAssemblerArgumentCollection();
    }

    /**
     * You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This
     * is especially useful if your samples are all in the same file but you need to run them individually through HC
     * in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     */
    @Argument(fullName = StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME, shortName = StandardArgumentDefinitions.SAMPLE_ALIAS_SHORT_NAME, doc = "Name of single sample to use from a multi-sample bam", optional = true)
    public String sampleNameToUse = null;

    /**
     * rsIDs from this file are used to populate the ID column of the output. Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    public DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Argument(fullName = StandardArgumentDefinitions.COMPARISON_LONG_NAME, shortName = StandardArgumentDefinitions.COMPARISON_SHORT_NAME, doc = "Comparison VCF file(s)", optional = true)
    public List<FeatureInput<VariantContext>> comps = new ArrayList<>();

    /**
     * When HC is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar genotype quality (GQ) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     *
     * This argument allows you to set the GQ bands. HC expects a list of strictly increasing GQ values
     * that will act as exclusive upper bounds for the GQ bands. To pass multiple values,
     * you provide them one by one with the argument, as in `-GQB 10 -GQB 20 -GQB 30` and so on
     * (this would set the GQ bands to be `[0, 10), [10, 20), [20, 30)` and so on, for example).
     * Note that GQ values are capped at 99 in the GATK, so values must be integers in [1, 100].
     * If the last value is strictly less than 100, the last GQ band will start at that value (inclusive)
     * and end at 100 (exclusive).
     */
    @Advanced
    @Argument(fullName = GQ_BAND_LONG_NAME, shortName = GQ_BAND_SHORT_NAME, doc= "Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)", optional = true)
    public List<Integer> GVCFGQBands = new ArrayList<>(70);
    {
            for (int i=1; i<=60; ++i) {
                GVCFGQBands.add(i);
            }
            GVCFGQBands.add(70); GVCFGQBands.add(80); GVCFGQBands.add(90); GVCFGQBands.add(99);
    };

    /**
     * Output the band lower bound for each GQ block instead of the min GQ -- for better compression
     */
    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS, doc = "Output the band lower bound for each GQ block regardless of the data it represents", optional = true)
    public boolean floorBlocks = false;

    /**
     * This parameter determines the maximum size of an indel considered as potentially segregating in the
     * reference model.  It is used to eliminate reads from being indel informative at a site, and determines
     * by that mechanism the certainty in the reference base.  Conceptually, setting this parameter to
     * X means that each informative read is consistent with any indel of size < X being present at a specific
     * position in the genome, given its alignment to the reference.
     */
    @Advanced
    @Argument(fullName = "indel-size-to-eliminate-in-ref-model", doc = "The size of an indel to check for in the reference model", optional = true)
    public int indelSizeToEliminateInRefModel = 10;

    /**
     * If set, certain "early exit" optimizations in HaplotypeCaller, which aim to save compute and time by skipping
     * calculations if an ActiveRegion is determined to contain no variants, will be disabled. This is most likely to be useful if
     * you're using the -bamout argument to examine the placement of reads following reassembly and are interested in seeing the mapping of
     * reads in regions with no variations. Setting the --force-active flag may also be necessary.
     */
    @Advanced
    @Argument(fullName = "disable-optimizations", doc="Don't skip calculations in ActiveRegions with no variants",
            optional = true)
    public boolean disableOptimizations = false;

    /**
     * These arguments are associated with DRAGEN-GATK
     */

    /**
     * DRAGEN-GATK mode changes a long list of arguments to support running DRAGEN-GATK with FRD + BQD + STRE (with or without
     * a provided STRE table provided):
     *
     *
     */
    @Argument(fullName = DRAGEN_GATK_MODE_LONG_NAME, optional = true, doc="Single argument for enabling the bulk of DRAGEN-GATK features. NOTE: THIS WILL OVERWRITE PROVIDED ARGUMENT CHECK TOOL INFO TO SEE WHICH ARGUMENTS ARE SET).")
    public Boolean dragenMode = false;
    @Advanced
    @Argument(fullName = "apply-bqd", doc = "If enabled this argument will apply the DRAGEN-GATK BaseQualityDropout model to the genotyping model for filtering sites due to Linked Error mode.", optional = true)
    public boolean applyBQD = false;
    @Advanced
    @Argument(fullName = "apply-frd", doc = "If enabled this argument will apply the DRAGEN-GATK ForeignReadDetection model to the genotyping model for filtering sites.", optional = true)
    public boolean applyFRD = false;
    @Advanced
    @Argument(fullName = "disable-spanning-event-genotyping", doc = "If enabled this argument will disable inclusion of the '*' spanning event when genotyping events that overlap deletions", optional = true)
    public boolean disableSpanningEventGenotyping = false;
    @Advanced
    @Argument(fullName = "transform-dragen-mapping-quality", doc = "If enabled this argument will map DRAGEN aligner aligned reads with mapping quality <=250 to scale up to MQ 50", optional = true)
    public boolean transformDRAGENMapQ = false;
    //TODO NOTE TO THE REVIEWER, THIS ARGUMENT IS INSUFFICIENT BOTH THIS AND --minimum-mapping-quality must be set, unfortunatley
    //TODO they can't be unified since the toolDefaultReadFilters get instantiated before this field gets populated, and we can't
    //TODO pull the threshold from that filter since it might or might not exist by the time we go to filter for threading, really
    //TODO we should unify on the readFilter version of this check i think but perhaps they are seperate for athropological historical reasons and it is thus culturally protected?
    @Advanced
    @Argument(fullName = "mapping-quality-threshold-for-genotyping", doc = "Control the threshold for discounting reads from the genotyper due to mapping quality after the active region detection and assembly steps but before genotyping. NOTE: this is in contrast to the --"+ ReadFilterArgumentDefinitions.MINIMUM_MAPPING_QUALITY_NAME+" argument which filters reads from all parts of the HaplotypeCaller. If you would like to call genotypes with a different threshold both arguments must be set.", optional = true)
    public int mappingQualityThreshold = HaplotypeCallerEngine.DEFAULT_READ_QUALITY_FILTER_THRESHOLD;
    @Advanced
    @Argument(fullName = "max-effective-depth-adjustment-for-frd", doc = "Set the maximum depth to modify FRD adjustment to in the event of high depth sites (0 to disable)", optional = false)
    public int maxEffectiveDepthAdjustment = 0;

    @Hidden
    @Argument(fullName = "keep-rg", doc = "Only use reads from this read group when making calls (but use all reads to build the assembly)", optional = true)
    public String keepRG = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName = "just-determine-active-regions", doc = "Just determine ActiveRegions, don't perform assembly or calling", optional = true)
    public boolean justDetermineActiveRegions = false;


    @Hidden
    @Advanced
    @Argument(fullName="debug-assembly-region-state", doc="Write output files for assembled regions with read summaries and called haplotypes to the specified path", optional = true)
    public GATKPath assemblyStateOutput = null;

    @Hidden
    @Advanced
    @Argument(fullName="debug-genotyper-output", doc ="Location to write genotyper debug stream that contains detailed information about the internal state of the genotyepr", optional = true)
    public String genotyperDebugOutStream = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName = "dont-genotype", doc = "Perform assembly but do not genotype variants", optional = true)
    public boolean dontGenotype = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName = DO_NOT_RUN_PHYSICAL_PHASING_LONG_NAME,  doc = "Disable physical phasing", optional = true)
    public boolean doNotRunPhysicalPhasing = false;

    /**
     * Base quality is capped at half of PCR error rate for bases where read and mate overlap, to account for full correlation of PCR errors at these bases.  This argument disables that correction.
     */
    @Advanced
    @Argument(fullName = DO_NOT_CORRECT_OVERLAPPING_BASE_QUALITIES_LONG_NAME, doc = "Disable overlapping base quality correction")
    public boolean doNotCorrectOverlappingBaseQualities = false;

    @Advanced
    @Argument(fullName= USE_FILTERED_READS_FOR_ANNOTATIONS_LONG_NAME, doc = "Use the contamination-filtered read maps for the purposes of annotating variants", optional=true)
    public boolean useFilteredReadMapForAnnotations = false;
}
