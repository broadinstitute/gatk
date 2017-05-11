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
     * The default value of this argument was chosen to achieve roughly one false positive per covered megabase according
     * to a back-of-the-envelope calculation assuming the neutral model of allele selection, with parameters estimated from ExAC.
     */
    @Argument(fullName = "max_germline_posterior", optional = true, doc = "Maximum posterior probability that an allele is a germline variant")
    public double maxGermlinePosterior = 0.025;

    @Argument(fullName = "maxAltAllelesThreshold", optional = true, doc="filter variants with too many alt alleles")
    public int numAltAllelesThreshold = 1;

    @Argument(fullName = "maxMedianMappingQualityDifference", optional = true, doc="filter variants for which alt reads' median mapping quality is much lower than ref reads' median mapping quality.")
    public int maxMedianMappingQualityDifference = 15;

    @Argument(fullName = "maxMedianBaseQualityDifference", optional = true, doc="filter variants for which alt reads' median base quality is much lower than ref reads' median base quality.")
    public int maxMedianBaseQualityDifference = 10;

    @Argument(fullName = "maxMedianClippingDifference", optional = true, doc="filter variants for which alt reads' median number of clipped bases is too high compared to the median for ref reads.")
    public int maxMedianClippingDifference = 1;

    @Argument(fullName = "maxMedianFragmentLengthDifference", optional = true, doc="filter variants for which alt reads' median fragment length is very different from the median for ref reads.")
    public int maxMedianFragmentLengthDifference = 10000;

    @Argument(fullName = "minMedianReadPosition", optional = true, doc="filter variants for which the median position of alt alleles within reads is too near the end of reads.")
    public int minMedianReadPosition = 5;

    @Argument(fullName = "maxEventsInHaplotype", optional = true, doc="Variants coming from a haplotype with more than this many events are filtered")
    public int maxEventsInHaplotype = 2;

    @Argument(shortName = "strand_prob", fullName = "strandArtifactPosteriorProbability", optional = true, doc = "Filter a variant if the probability of strand artifact exceeds this number")
    public double STRAND_ARTIFACT_POSTERIOR_PROB_THRESHOLD = 0.99;

    @Argument(shortName = "strand_af", fullName = "strandArtifactAlleleFraction", optional = true, doc = "Only filter a variant if the MAP estimate of allele fraction given artifact is below this number")
    public double STRAND_ARTIFACT_ALLELE_FRACTION_THRESHOLD = 0.01;

    @Argument(shortName = "contaminationTable", fullName = "contaminationTable", optional = true, doc="Table containing contamination information.")
    public File contaminationTable = null;

}
