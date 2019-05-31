package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public abstract class AssemblyRegionArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "min-assembly-region-size", doc = "Minimum size of an assembly region", optional = true)
    public int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Argument(fullName = "max-assembly-region-size", doc = "Maximum size of an assembly region", optional = true)
    public int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Argument(fullName = "assembly-region-padding", doc = "Number of additional bases of context to include around each assembly region", optional = true)
    public int assemblyRegionPadding = defaultAssemblyRegionPadding();

    @Argument(fullName = "max-reads-per-alignment-start", doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    public int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

    @Advanced
    @Argument(fullName = "active-probability-threshold", doc="Minimum probability for a locus to be considered active.", optional = true)
    public double activeProbThreshold = defaultActiveProbThreshold();

    @Advanced
    @Argument(fullName = "max-prob-propagation-distance", doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
    public int maxProbPropagationDistance = defaultMaxProbPropagationDistance();

    /**
     * @return Default value for the {@link #minAssemblyRegionSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultMinAssemblyRegionSize();

    /**
     * @return Default value for the {@link #maxAssemblyRegionSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxAssemblyRegionSize();

    /**
     * @return Default value for the {@link #assemblyRegionPadding} parameter, if none is provided on the command line
     */
    protected abstract int defaultAssemblyRegionPadding();

    /**
     * @return Default value for the {@link #maxReadsPerAlignmentStart} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxReadsPerAlignmentStart();

    /**
     * @return Default value for the {@link #activeProbThreshold} parameter, if none is provided on the command line
     */
    protected abstract double defaultActiveProbThreshold();

    /**
     * @return Default value for the {@link #maxProbPropagationDistance} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxProbPropagationDistance();
}
