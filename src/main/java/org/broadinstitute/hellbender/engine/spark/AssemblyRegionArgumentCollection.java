package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;

public class AssemblyRegionArgumentCollection implements Serializable {
    public static final String ASSEMBLY_REGION_OUT_LONG_NAME = "assembly-region-out";
    public static final String FORCE_ACTIVE_REGIONS_LONG_NAME = "force-active";
    private static final long serialVersionUID = 1L;

    public static final String MIN_ASSEMBLY_LONG_NAME = "min-assembly-region-size";
    public static final String MAX_ASSEMBLY_LONG_NAME = "max-assembly-region-size";
    public static final String ASSEMBLY_PADDING_LONG_NAME = "assembly-region-padding";
    public static final String MAX_STARTS_LONG_NAME = "max-reads-per-alignment-start";
    public static final String THRESHOLD_LONG_NAME = "active-probability-threshold";
    public static final String PROPAGATION_LONG_NAME = "max-prob-propagation-distance";

    public static final int DEFAULT_MIN_ASSEMBLY_REGION_SIZE = 50;
    public static final int DEFAULT_MAX_ASSEMBLY_REGION_SIZE = 300;
    public static final int DEFAULT_ASSEMBLY_REGION_PADDING = 100;
    public static final int DEFAULT_MAX_READS_PER_ALIGNMENT = 50;
    public static final double DEFAULT_ACTIVE_PROB_THRESHOLD = 0.002;
    public static final int DEFAULT_MAX_PROB_PROPAGATION_DISTANCE = 50;
    public static final String INDEL_PADDING_LONG_NAME = "padding-around-indels";
    public static final String SNP_PADDING_LONG_NAME = "padding-around-snps";
    public static final String STR_PADDING_LONG_NAME = "padding-around-strs";

    /**
     * The following parameters can be confusing due to the overlap between active regions, assembly regions, and genotyping regions,
     * all of which are implemented as the {@link org.broadinstitute.hellbender.engine.AssemblyRegion} class.
     *
     * An active region is a genomic interval in which we call variants.  In order to output phased calls we try to put nearby variants
     * in the same active region.  The size of active regions are determined by the variants present, {@link AssemblyRegionArgumentCollection.activeProbThreshold},
     * and {@link AssemblyRegionArgumentCollection.maxProbPropagationDistance}, and are bounded by the parameters {@link AssemblyRegionArgumentCollection.minAssemblyRegionSize}
     * and {@link AssemblyRegionArgumentCollection.maxAssemblyRegionSize}.  Importantly, active regions define responsibility for variant calling
     * and have nothing to do with the size of the span over which we perform local assembly and pair-HMM.  That is, a variant may belong to two
     * overlapping assembly or genotyping regions but it only belongs to a single active region.
     *
     * An assembly region is the padded interval surrounding an active region over which we perform local assembly.  Padding is useful to phase with
     * any variation that may be just outside the active region, to avoid dangling ends in the assembly region, and to resolve indels.  The
     * {@code assemblyRegionPadding} parameter determines the number of extra bases an assembly region contains on either side of the active region
     * it surrounds.
     *
     * A genotyping region is an interval surrounding an active region in which we perform pair-HMM and Smith-Waterman alignment.  Even though
     * we ultimately only genotype and call variants within the active region, we call entire haplotypes as an intermediate step, hence the
     * needed for an expanded genotyping region.  The parameters {@link AssemblyRegionArgumentCollection.snpPaddingForGenotyping},
     * {@link AssemblyRegionArgumentCollection.snpPaddingForGenotyping} and {@link AssemblyRegionArgumentCollection.indelPaddingForGenotyping} determine
     * the size of genotyping regions.
     *
     * The overall flow is as follows:
     *
     * 1) Active regions are found by triaging sites for activity and combining nearby potential variants into active regions.
     * 2) Active regions are padded for assembly and overlapping reads are hard-clipped to this padded assembly region.
     * 3) Reads are assembled into haplotypes.
     * 4) The active region is re-set to the span of all variants found in assembly that overlap the original active region.
     *    This shrinking does not lose variants but it does affect the size of the genotyping region.
     * 5) The genotyping region is determined by padding the shrunk active region.
     * 6) Assembled haplotypes are trimmed and reads are hard-clipped (again) to fit the padded genotyping region.
     */


    /**
     * Parameters that control active regions
     */

    @Argument(fullName = MIN_ASSEMBLY_LONG_NAME, doc = "Minimum size of an assembly region", optional = true)
    public int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Argument(fullName = MAX_ASSEMBLY_LONG_NAME, doc = "Maximum size of an assembly region", optional = true)
    public int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Advanced
    @Argument(fullName = THRESHOLD_LONG_NAME, doc="Minimum probability for a locus to be considered active.", optional = true)
    public double activeProbThreshold = defaultActiveProbThreshold();

    @Advanced
    @Argument(fullName = PROPAGATION_LONG_NAME, doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
    public int maxProbPropagationDistance = defaultMaxProbPropagationDistance();

    @Advanced
    @Argument(fullName = FORCE_ACTIVE_REGIONS_LONG_NAME, doc = "If provided, all regions will be marked as active", optional = true)
    public boolean forceActive = false;

    /**
     * Parameters that control assembly regions
     */

    @Argument(fullName = ASSEMBLY_PADDING_LONG_NAME, doc = "Number of additional bases of context to include around each assembly region", optional = true)
    public int assemblyRegionPadding = defaultAssemblyRegionPadding();

    /**
     * Parameters that control genotyping regions
     */

    @Hidden
    @Argument(fullName= INDEL_PADDING_LONG_NAME, doc = "Include at least this many bases around an event for calling indels", optional = true)
    public int indelPaddingForGenotyping = 75;

    @Hidden
    @Argument(fullName= SNP_PADDING_LONG_NAME, doc = "Include at least this many bases around an event for calling snps", optional = true)
    public int snpPaddingForGenotyping = 20;

    @Hidden
    @Argument(fullName= STR_PADDING_LONG_NAME, doc = "Include at least this many bases around an event for calling STR indels", optional = true)
    public int strPaddingForGenotyping = 75;

    /**
     * The maximum extent into the full active region extension that we're willing to go in genotyping our events
     * NOTE: this is only applicable to the legacy assembly region trimming currently
     */
    @Hidden
    @Argument(fullName="max-extension-into-assembly-region-padding-legacy", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping (-1 to disable). NOTE this only applies for --enable-legacy-assembly-region-trimming mode", optional = true)
    public int maxExtensionIntoRegionPadding = 25;

    /**
     * Other parameters
     */

    @Argument(fullName = MAX_STARTS_LONG_NAME, doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    public int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

    @Hidden
    @Argument(fullName = "enable-legacy-assembly-region-trimming", doc = "Revert changes to the assembly region windows, this will result in less consistent results for assembly window boundaries", optional = true)
    public boolean enableLegacyAssemblyRegionTrimming = false;

    /**
     * @return Default value for the {@link #minAssemblyRegionSize} parameter, if none is provided on the command line
     */
    protected int defaultMinAssemblyRegionSize() { return DEFAULT_MIN_ASSEMBLY_REGION_SIZE; }

    /**
     * @return Default value for the {@link #maxAssemblyRegionSize} parameter, if none is provided on the command line
     */
    protected int defaultMaxAssemblyRegionSize() { return DEFAULT_MAX_ASSEMBLY_REGION_SIZE; }

    /**
     * @return Default value for the {@link #assemblyRegionPadding} parameter, if none is provided on the command line
     */
    protected int defaultAssemblyRegionPadding() { return DEFAULT_ASSEMBLY_REGION_PADDING; }

    /**
     * @return Default value for the {@link #maxReadsPerAlignmentStart} parameter, if none is provided on the command line
     */
    protected int defaultMaxReadsPerAlignmentStart() { return DEFAULT_MAX_READS_PER_ALIGNMENT; }

    /**
     * @return Default value for the {@link #activeProbThreshold} parameter, if none is provided on the command line
     */
    protected double defaultActiveProbThreshold() { return DEFAULT_ACTIVE_PROB_THRESHOLD; }

    /**
     * @return Default value for the {@link #maxProbPropagationDistance} parameter, if none is provided on the command line
     */
    protected int defaultMaxProbPropagationDistance() { return DEFAULT_MAX_PROB_PROPAGATION_DISTANCE; }

    public void validate() {
        if ( minAssemblyRegionSize <= 0 || maxAssemblyRegionSize <= 0 ) {
            throw new CommandLineException.BadArgumentValue("min/max assembly region size must be > 0");
        }

        if ( minAssemblyRegionSize > maxAssemblyRegionSize ) {
            throw new CommandLineException.BadArgumentValue("minAssemblyRegionSize must be <= maxAssemblyRegionSize");
        }

        if ( assemblyRegionPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("assemblyRegionPadding must be >= 0");
        }

        if ( maxReadsPerAlignmentStart < 0 ) {
            throw new CommandLineException.BadArgumentValue("maxReadsPerAlignmentStart must be >= 0");
        }

        if ( snpPaddingForGenotyping < 0 ) {
            throw new CommandLineException.BadArgumentValue("paddingAroundSNPs", "" + snpPaddingForGenotyping + "< 0");
        }

        if ( indelPaddingForGenotyping < 0 ) {
            throw new CommandLineException.BadArgumentValue("paddingAroundIndels", "" + indelPaddingForGenotyping + "< 0");
        }
    }


}
