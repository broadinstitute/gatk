package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;

import java.io.File;

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

    @Argument(fullName = "maxAltAllelesThreshold", optional = true, doc="filter variants with too many alt alleles")
    public int numAltAllelesThreshold = 1;

    @Argument(fullName = "maxMedianMappingQualityDifference", optional = true, doc="filter variants for which alt reads' median mapping quality is much lower than ref reads' median mapping quality.")
    public int maxMedianMappingQualityDifference = 15;

    @Argument(fullName = "maxMedianBaseQualityDifference", optional = true, doc="filter variants for which alt reads' median base quality is much lower than ref reads' median base quality.")
    public int maxMedianBaseQualityDifference = 10;

    @Argument(fullName = "maxMedianClippingDifference", optional = true, doc="filter variants for which alt reads' median number of clipped bases is too high compared to the median for ref reads.")
    public int maxMedianClippingDifference = 1;

    @Argument(fullName = "maxMedianFragmentLengthDifference", optional = true, doc="filter variants for which alt reads' median fragment length is very different from the median for ref reads.")
    public int maxMedianFragmentLengthDifference = 200;

    @Argument(fullName = "minMedianReadPosition", optional = true, doc="filter variants for which the median position of alt alleles within reads is too near the end of reads.")
    public int minMedianReadPosition = 5;


    @Argument(fullName = "maxEventsInHaplotype", optional = true, doc="Variants coming from a haplotype with more than this many events are filtered")
    public int maxEventsInHaplotype = 2;

    @Argument(shortName = "contaminationTable", fullName = "contaminationTable", optional = true, doc="Table containing contamination information.")
    public File contaminationTable = null;

}
