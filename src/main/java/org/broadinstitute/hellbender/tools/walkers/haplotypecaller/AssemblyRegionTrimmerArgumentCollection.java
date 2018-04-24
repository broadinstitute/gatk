package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;

public class AssemblyRegionTrimmerArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Advanced
    @Argument(fullName="dont-trim-active-regions", doc="If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping", optional = true)
    protected boolean dontTrimActiveRegions = false;

    /**
     * the maximum extent into the full active region extension that we're willing to go in genotyping our events
     */
    @Hidden
    @Argument(fullName="max-disc-ar-extension", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping our events for discovery", optional = true)
    protected int discoverExtension = 25;

    @Hidden
    @Argument(fullName="max-gga-ar-extension", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping our events for GGA mode", optional = true)
    protected int ggaExtension = 300;

    /**
     * Include at least this many bases around an event for calling it
     */
    @Hidden
    @Argument(fullName="padding-around-indels", doc = "Include at least this many bases around an event for calling indels", optional = true)
    public int indelPadding = 150;

    @Hidden
    @Argument(fullName="padding-around-snps", doc = "Include at least this many bases around an event for calling snps", optional = true)
    public int snpPadding = 20;
}
