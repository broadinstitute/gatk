package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;

/**
 * Mutect2's assembly region arguments, whose default activity profile is MUTECT2 and which include the arguments that only
 * the MUTECT2 activity profile uses
 */
public class Mutect2AssemblyRegionArgumentCollection extends AssemblyRegionArgumentCollection {
    private static final long serialVersionUID = 1L;

    public static final String GENOTYPE_GERMLINE_SITES_FRACTION_LONG_NAME = "genotype-germline-sites-fraction";
    public static final double DEFAULT_GENOTYPE_GERMLINE_SITES_FRACTION = 1.0;

    /**
     * For the purposes of learning some parameters in Permutect it may be useful to genotype *some* germline variants
     * and emit them as Permutect tensors.  When --genotype-germline-sites is true, setting this parameter to a value less
     * than the default of 1.0 causes only that fraction of clusters of germline variants without nearby somatic variants
     * to be genotyped.  When --genotype-germline-sites is false this argument has no effect and no germline sites are
     * genotyped.  Values other than 1.0 require the MUTECT2 activity profile.
     */
    @Argument(fullName= GENOTYPE_GERMLINE_SITES_FRACTION_LONG_NAME, doc="Fraction of germline sites to be genotyped randomly.", optional = true)
    public double genotypeGermlineSitesFraction = DEFAULT_GENOTYPE_GERMLINE_SITES_FRACTION;

    /**
     * @return Default value for the {@link #activityProfileType} parameter, if none is provided on the command line
     */
    @Override
    protected ActivityProfile.ProfileType defaultActivityProfileType() { return ActivityProfile.ProfileType.MUTECT2; }

    @Override
    protected boolean supportsMutect2ActivityProfile() { return true; }
}
