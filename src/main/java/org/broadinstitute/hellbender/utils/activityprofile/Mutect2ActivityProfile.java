package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;

public class Mutect2ActivityProfile extends ActivityProfile {
    public Mutect2ActivityProfile(final int maxProbPropagationDistance,
                                   final double activeProbThreshold,
                                   final SAMFileHeader samHeader) {
        super(maxProbPropagationDistance, activeProbThreshold, samHeader);
    }
}
