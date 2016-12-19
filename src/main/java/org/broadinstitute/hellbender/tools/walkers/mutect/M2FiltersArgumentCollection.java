package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

public class M2FiltersArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9345L;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal
     * as an artifact i.e. not as a germline event.
     */
    @Argument(fullName = "normal_artifact_lod", optional = true, doc = "LOD threshold for calling normal artifacts")
    public double NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0;

    /**
     * The LOD threshold for the normal is typically made more strict if the variant has been seen in dbSNP (i.e. another
     * normal sample). We thus require MORE evidence that a variant is NOT seen in this tumor's normal if it has been observed as a germline variant before.
     */
    @Argument(fullName = "dbsnp_normal_lod", optional = true, doc = "LOD threshold for calling normal non-variant at dbsnp sites")
    public double NORMAL_DBSNP_LOD_THRESHOLD = 5.5;

    @Hidden
    @Argument(fullName = "strand_artifact_lod", optional = true, doc = "LOD threshold for calling strand bias")
    public float STRAND_ARTIFACT_LOD_THRESHOLD = 2.0f;

    @Hidden
    @Argument(fullName = "strand_artifact_power_threshold", optional = true, doc = "power threshold for calling strand bias")
    public float STRAND_ARTIFACT_POWER_THRESHOLD = 0.9f;

    @Argument(fullName = "enable_strand_artifact_filter", optional = true, doc = "turn on strand artifact filter")
    public boolean ENABLE_STRAND_ARTIFACT_FILTER = false;

    @Argument(fullName = "enable_clustered_read_position_filter", optional = true, doc = "turn on clustered read position filter")
    public boolean ENABLE_CLUSTERED_READ_POSITION_FILTER = false;

    /**
     * This argument is used for the M1-style read position filter
     */
    @Argument(fullName = "pir_median_threshold", optional = true, doc="threshold for clustered read position artifact median")
    public double PIR_MEDIAN_THRESHOLD = 10;

    /**
     * This argument is used for the M1-style read position filter
     */
    @Argument(fullName = "pir_mad_threshold", optional = true, doc="threshold for clustered read position artifact MAD")
    public double PIR_MAD_THRESHOLD = 3;
}
