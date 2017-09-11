package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;

import java.io.File;

public class M2FiltersArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9345L;

    /**
     * Only variants with tumor LODs exceeding this threshold can pass filtering.
     */
    @Argument(fullName = "tumor-lod", optional = true, doc = "LOD threshold for calling tumor variant")
    public double TUMOR_LOD_THRESHOLD = 5.3;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal
     * as an artifact i.e. not as a germline event.
     */
    @Argument(fullName = "normal-artifact-lod", optional = true, doc = "LOD threshold for calling normal artifacts")
    public double NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0;

    /**
     * The default value of this argument was chosen to achieve roughly one false positive per covered megabase according
     * to a back-of-the-envelope calculation assuming the neutral model of allele selection, with parameters estimated from ExAC.
     */
    @Argument(fullName = "max-germline-posterior", optional = true, doc = "Maximum posterior probability that an allele is a germline variant")
    public double maxGermlinePosterior = 0.025;

    @Argument(fullName = "max-alt-allele-count", optional = true, doc = "filter variants with too many alt alleles")
    public int numAltAllelesThreshold = 1;

    @Argument(fullName = "min-median-mapping-quality", optional = true, doc="filter variants for which alt reads' median mapping quality is too low.")
    public int minMedianMappingQuality = 30;

    @Argument(fullName = "min-median-base-quality", optional = true, doc="filter variants for which alt reads' median base quality is too low.")
    public int minMedianBaseQuality = 20;

    @Argument(fullName = "max-median-fragment-length-difference", optional = true, doc="filter variants for which alt reads' median fragment length is very different from the median for ref reads.")
    public int maxMedianFragmentLengthDifference = 10000;

    @Argument(fullName = "min-median-read-position", optional = true, doc = "filter variants for which the median position of alt alleles within reads is too near the end of reads.")
    public int minMedianReadPosition = 5;

    @Argument(fullName = "max-events-in-region", optional = true, doc = "Variants coming from an assembly region with more than this many events are filtered")
    public int maxEventsInRegion = 2;

    @Argument(fullName = "max-strand-artifact-probability", shortName = "strand-prob", optional = true, doc = "Filter a variant if the probability of strand artifact exceeds this number")
    public double strandArtifactPosteriorProbThreshold = 0.99;

    @Argument(fullName = "min-strand-artifact-allele-fraction", shortName = "strand-af", optional = true, doc = "Only filter a variant if the MAP estimate of allele fraction given artifact is below this number")
    public double strandArtifactAlleleFractionThreshold = 0.01;

    @Argument(fullName = "contamination-table", optional = true, doc = "Table containing contamination information.")
    public File contaminationTable = null;

    @Argument(fullName = "unique-alt-read-count", shortName = "unique", optional = true, doc = "Filter a variant if a site contains fewer than this many unique (i.e. deduplicated) reads supporting the alternate allele")
    public int uniqueAltReadCount = 0;

}
