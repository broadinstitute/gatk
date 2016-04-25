package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.Hidden;

public class AssemblyRegionTrimmerArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    @Advanced
    @Argument(fullName="dontTrimActiveRegions", shortName="dontTrimActiveRegions", doc="If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping", optional = true)
    protected boolean dontTrimActiveRegions = false;

    /**
     * the maximum extent into the full active region extension that we're willing to go in genotyping our events
     */
    @Hidden
    @Argument(fullName="maxDiscARExtension", shortName="maxDiscARExtension", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping our events for discovery", optional = true)
    protected int discoverExtension = 25;

    @Hidden
    @Argument(fullName="maxGGAARExtension", shortName="maxGGAARExtension", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping our events for GGA mode", optional = true)
    protected int ggaExtension = 300;

    /**
     * Include at least this many bases around an event for calling it
     */
    @Hidden
    @Argument(fullName="paddingAroundIndels", shortName="paddingAroundIndels", doc = "Include at least this many bases around an event for calling indels", optional = true)
    public int indelPadding = 150;

    @Hidden
    @Argument(fullName="paddingAroundSNPs", shortName="paddingAroundSNPs", doc = "Include at least this many bases around an event for calling snps", optional = true)
    public int snpPadding = 20;
}
