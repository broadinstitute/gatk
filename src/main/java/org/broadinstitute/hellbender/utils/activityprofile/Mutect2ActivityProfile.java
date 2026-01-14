package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;

import java.util.Collection;
import java.util.Collections;

public class Mutect2ActivityProfile extends ActivityProfile {
    public Mutect2ActivityProfile(final int maxProbPropagationDistance,
                                   final double activeProbThreshold,
                                   final SAMFileHeader samHeader) {
        super(maxProbPropagationDistance, activeProbThreshold, samHeader);
    }

    /**
     * Don't do anything to the state! Just add it!
     * @param justAddedState the most recent ActivityProfileState in traversal
     */
    @Override
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        return Collections.singletonList(justAddedState);
    }
}
