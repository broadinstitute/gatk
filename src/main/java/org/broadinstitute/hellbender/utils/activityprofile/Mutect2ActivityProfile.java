package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;

public class Mutect2ActivityProfile extends ActivityProfile {
    public Mutect2ActivityProfile(final AssemblyRegionArgumentCollection args, final SAMFileHeader header) {
        super(args.maxProbPropagationDistance, args.activeProbThreshold, header);
    }
}