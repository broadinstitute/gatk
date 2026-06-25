package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;

import java.util.Collection;
import java.util.Collections;

public class Mutect2ActivityProfile extends ActivityProfile {
    public Mutect2ActivityProfile(final AssemblyRegionArgumentCollection args, final SAMFileHeader header) {
        super(args.maxProbPropagationDistance, args.activeProbThreshold, header);
    }

    @Override
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        return Collections.singletonList(justAddedState);
    }
}