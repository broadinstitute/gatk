package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;

public class Mutect2AssemblyRegionArgumentCollection extends AssemblyRegionArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Override
    protected ActivityProfile.ProfileType defaultActivityProfileType() { return ActivityProfile.ProfileType.MUTECT2; }

    @Override
    protected double defaultGenotypeGermlineSitesFraction() { return 0.0; }
}
