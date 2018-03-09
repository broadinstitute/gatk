package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

import java.io.File;

public class M2FiltersArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9345L;
    public static final String LOG_SOMATIC_PRIOR_LONG_NAME = "log-somatic-prior";
    public static final String TUMOR_LOD_LONG_NAME = "tumor-lod";
    public static final String NORMAL_ARTIFACT_LOD_LONG_NAME = "normal-artifact-lod";
    public static final String MAX_GERMLINE_POSTERIOR_LONG_NAME = "max-germline-posterior";
    public static final String MAX_ALT_ALLELE_COUNT_LONG_NAME = "max-alt-allele-count";
    public static final String MIN_MEDIAN_MAPPING_QUALITY_LONG_NAME = "min-median-mapping-quality";
    public static final String MIN_MEDIAN_BASE_QUALITY_LONG_NAME = "min-median-base-quality";
    public static final String MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_LONG_NAME = "max-median-fragment-length-difference";
    public static final String MIN_MEDIAN_READ_POSITION_LONG_NAME = "min-median-read-position";
    public static final String MAX_EVENTS_IN_REGION_LONG_NAME = "max-events-in-region";
    public static final String MAX_STRAND_ARTIFACT_PROBABILITY_LONG_NAME = "max-strand-artifact-probability";
    public static final String MIN_STRAND_ARTIFACT_ALLELE_FRACTION_LONG_NAME = "min-strand-artifact-allele-fraction";
    public static final String CONTAMINATION_TABLE_LONG_NAME = "contamination-table";
    public static final String MAX_CONTAMINATION_PROBABILITY_LONG_NAME = "max-contamination-probability";
    public static final String UNIQUE_ALT_READ_COUNT_LONG_NAME = "unique-alt-read-count";

    public static final String TUMOR_SEGMENTATION_LONG_NAME = "tumor-segmentation";

    /**
     * A table containing tumor segments and the minor allele fraction of germline hets within each segment.
     * This allows us to refine the germline event filter by, for example, not filtering an allele
     * with an allele fraction significantly different from 0.5 in a segment where the minor allele fraction is 0.5.
     */
    @Argument(fullName = TUMOR_SEGMENTATION_LONG_NAME,
            doc="Pileup summaries for the tumor sample as output by CalculateContamination", optional = true)
    public File tumorSegmentationTable = null;

    /**
     * Prior log-10 probability that any given site has a somatic allele. Impacts germline probability calculation.
     * The workflow uses this parameter only towards the germline event filter. It does NOT relate to the LOD threshold.
     * For example, -6 translates to one in a million or ~3000 somatic mutations per human genome.
     * Depending on tumor type, mutation rate ranges vary (Lawrence et al. Nature 2013), and so adjust parameter accordingly.
     * For higher expected rate of mutation, adjust number up, e.g. -5. For lower expected rate of mutation, adjust number down, e.g. -7.
     */
    @Argument(fullName= LOG_SOMATIC_PRIOR_LONG_NAME,
            doc="Prior probability that a given site has a somatic allele.", optional = true)
    public double log10PriorProbOfSomaticEvent = -6.0;

    /**
     * Only variants with tumor LODs exceeding this threshold can pass filtering.
     */
    @Argument(fullName = TUMOR_LOD_LONG_NAME, optional = true, doc = "LOD threshold for calling tumor variant")
    public double TUMOR_LOD_THRESHOLD = 5.3;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal
     * as an artifact i.e. not as a germline event.
     */
    @Argument(fullName = NORMAL_ARTIFACT_LOD_LONG_NAME, optional = true, doc = "LOD threshold for calling normal artifacts")
    public double NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0;

    /**
     * The default value of this argument was chosen to achieve roughly one false positive per covered megabase according
     * to a back-of-the-envelope calculation assuming the neutral model of allele selection, with parameters estimated from ExAC.
     */
    @Argument(fullName = MAX_GERMLINE_POSTERIOR_LONG_NAME, optional = true, doc = "Maximum posterior probability that an allele is a germline variant")
    public double maxGermlinePosterior = 0.025;

    @Argument(fullName = MAX_ALT_ALLELE_COUNT_LONG_NAME, optional = true, doc = "filter variants with too many alt alleles")
    public int numAltAllelesThreshold = 1;

    @Argument(fullName = MIN_MEDIAN_MAPPING_QUALITY_LONG_NAME, optional = true, doc="filter variants for which alt reads' median mapping quality is too low.")
    public int minMedianMappingQuality = 30;

    @Argument(fullName = MIN_MEDIAN_BASE_QUALITY_LONG_NAME, optional = true, doc="filter variants for which alt reads' median base quality is too low.")
    public int minMedianBaseQuality = 20;

    @Argument(fullName = MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_LONG_NAME, optional = true, doc="filter variants for which alt reads' median fragment length is very different from the median for ref reads.")
    public int maxMedianFragmentLengthDifference = 10000;

    @Argument(fullName = MIN_MEDIAN_READ_POSITION_LONG_NAME, optional = true, doc = "filter variants for which the median position of alt alleles within reads is too near the end of reads.")
    public int minMedianReadPosition = 5;

    @Argument(fullName = MAX_EVENTS_IN_REGION_LONG_NAME, optional = true, doc = "Variants coming from an assembly region with more than this many events are filtered")
    public int maxEventsInRegion = 2;

    @Argument(fullName = MAX_STRAND_ARTIFACT_PROBABILITY_LONG_NAME, shortName = "strand-prob", optional = true, doc = "Filter a variant if the probability of strand artifact exceeds this number")
    public double strandArtifactPosteriorProbThreshold = 0.99;

    @Argument(fullName = MIN_STRAND_ARTIFACT_ALLELE_FRACTION_LONG_NAME, shortName = "strand-af", optional = true, doc = "Only filter a variant if the MAP estimate of allele fraction given artifact is below this number")
    public double strandArtifactAlleleFractionThreshold = 0.01;

    @Argument(fullName = CONTAMINATION_TABLE_LONG_NAME, optional = true, doc = "Table containing contamination information.")
    public File contaminationTable = null;

    @Argument(fullName = MAX_CONTAMINATION_PROBABILITY_LONG_NAME, optional = true, doc = "Filter variants with posterior probability to be due to contamination greater than this.")
    public double maxContaminationProbability = 0.1;

    @Argument(fullName = UNIQUE_ALT_READ_COUNT_LONG_NAME, shortName = "unique", optional = true, doc = "Filter a variant if a site contains fewer than this many unique (i.e. deduplicated) reads supporting the alternate allele")
    public int uniqueAltReadCount = 0;

}
