package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;

import java.util.Collection;
import java.util.Collections;

/**
 * Somatic-specific alternative to the default BandPassActivityProfile.  It is based on the following principles:
 *
 *  1) Two events separated by the maximum phasing distance or less should be put in the same active region in
 *      order to be assembled together.
 *  2) Germline variants are not inherently events of interest BUT they should be included in
 *     the active region of nearby somatic variants.
 *
 *  For example, suppose the phasing distance is 100, there are somatic variants at loci 200, 350, and 600; and
 *  there are germline variants at loci 280, 650, and 800.  The somatic variants at 200 and 350 will be in the same
 *  active region along with the germline variant at 280; the somatic variant at 600 will be grouped with the
 *  germline variant at 650; and the germline variant at 800 will not be in an active region.
 */
public class Mutect2ActivityProfile extends ActivityProfile {
    public Mutect2ActivityProfile(final AssemblyRegionArgumentCollection args, final SAMFileHeader header) {
        // In this class maxProbPropagationDistance is interpreted as a maximum phasing distance.
        // It doesn't "propagate" probability like the BandPassActivityFilter
        super(args.maxProbPropagationDistance, args.activeProbThreshold, header);
    }

    @Override
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        return Collections.singletonList(justAddedState);
    }

    // if two active loci are separated by this distance or less, they go in the same AssemblyRegion
    // eg if maxPhasingDistance is 10, sites 1 and 11 are assembled together but 1 and 12 are not.
    private int maxPhasingDistance() { return getMaxProbPropagationDistance(); }
}