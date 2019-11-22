package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;

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

    @Argument(fullName = MIN_ASSEMBLY_LONG_NAME, doc = "Minimum size of an assembly region", optional = true)
    public int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Argument(fullName = MAX_ASSEMBLY_LONG_NAME, doc = "Maximum size of an assembly region", optional = true)
    public int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Argument(fullName = ASSEMBLY_PADDING_LONG_NAME, doc = "Number of additional bases of context to include around each assembly region", optional = true)
    public int assemblyRegionPadding = defaultAssemblyRegionPadding();

    @Argument(fullName = MAX_STARTS_LONG_NAME, doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    public int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

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
    }


}
