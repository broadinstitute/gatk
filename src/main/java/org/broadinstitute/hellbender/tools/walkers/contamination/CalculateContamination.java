package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * <p>
 *     Calculates the fraction of reads coming from cross-sample contamination, given results from {@link GetPileupSummaries}.
 *     The resulting contamination table is used with {@link FilterMutectCalls}.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 *
 * <p>This tool borrows from <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/">ContEst</a> by Cibulskis et al the idea of estimating contamination
 * from ref reads at hom alt sites.  However, ContEst uses a probabilistic model that assumes a diploid genotype with no copy number
 * variation and independent contaminating reads.  That is, ContEst assumes that each contaminating read is drawn randomly and
 * independently from a different human.  This tool uses a simpler estimate of contamination that relaxes these assumptions.  In particular,
 * it works in the presence of copy number variations and with an arbitrary number of contaminating samples.  In addition, this tool
 * is designed to work well with no matched normal data.  However, one can run {@link GetPileupSummaries} on a matched normal bam file
 * and input the result to this tool.</p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Tumor-only mode</h4>
 *
 * <pre>
 * gatk CalculateContamination \
 *   -I pileups.table \
 *   -O contamination.table
 * </pre>
 *
 * <h4>Matched normal mode</h4>
 *
 * <pre>
 * gatk CalculateContamination \
 *   -I tumor-pileups.table \
 *   -matched normal-pileups.table \
 *   -O contamination.table
 * </pre>
 * <p>
 *     The resulting table provides the fraction contamination, one line per sample, e.g. SampleID--TAB--Contamination.
 *     The file has no header.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Calculate the fraction of reads coming from cross-sample contamination",
        oneLineSummary = "Calculate the fraction of reads coming from cross-sample contamination",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CalculateContamination extends CommandLineProgram {

    public static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    private static final int MIN_COVERAGE = 10;
    private static final double DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD = 1.0/2;
    private static final double DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD = 3.0;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table")
    private File inputPileupSummariesTable;

    public static final String MATCHED_NORMAL_LONG_NAME = "matched-normal";
    public static final String MATCHED_NORMAL_SHORT_NAME = "matched";
    @Argument(fullName = MATCHED_NORMAL_LONG_NAME,
            shortName = MATCHED_NORMAL_SHORT_NAME,
            doc="The matched normal input table", optional = true)
    private File matchedPileupSummariesTable = null;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table")
    private final File outputTable = null;

    public static final String TUMOR_SEGMENTATION_LONG_NAME = "tumor-segmentation";
    public static final String TUMOR_SEGMENTATION_SHORT_NAME = "segments";
    @Argument(fullName= TUMOR_SEGMENTATION_LONG_NAME,
            shortName= TUMOR_SEGMENTATION_SHORT_NAME,
            doc="The output table containing segmentation of the tumor by minor allele fraction", optional = true)
    private final File outputTumorSegmentation = null;

    public static final String LOW_COVERAGE_RATIO_THRESHOLD_NAME = "low-coverage-ratio-threshold";
    @Argument(fullName = LOW_COVERAGE_RATIO_THRESHOLD_NAME,
            doc="The minimum coverage relative to the median.", optional = true)
    private final double lowCoverageRatioThreshold = DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD;

    public static final String HIGH_COVERAGE_RATIO_THRESHOLD_NAME = "high-coverage-ratio-threshold";
    @Argument(fullName= HIGH_COVERAGE_RATIO_THRESHOLD_NAME,
            doc="The maximum coverage relative to the mean.", optional = true)
    private final double highCoverageRatioThreshold = DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD;

    public static final String HOM_SITES_FOR_CALCULATION_NAME = "hom-sites";
    @Argument(fullName = HOM_SITES_FOR_CALCULATION_NAME, doc = "homozygous sites (alt or ref) used for calculating contamination", optional = true)
    private final File homSitesFile = null;

    public static final String HIGH_COVERAGE_SITES_NAME = "high-coverage-sites";
    @Argument(fullName = HIGH_COVERAGE_SITES_NAME, optional = true)
    private final File highCoverageSitesFile = null;

    public static final String AUXILIARY_INFO = "aux";
    @Argument(fullName = AUXILIARY_INFO, optional = true)
    private final File auxiliaryInfoFile = null;


    @Override
    public Object doWork() {
        final Pair<String, List<PileupSummary>> sampleAndsites = PileupSummary.readFromFile(inputPileupSummariesTable);
        final String sample = sampleAndsites.getLeft();
        final List<PileupSummary> sites = filterSitesByCoverage(sampleAndsites.getRight());

        // used the matched normal to genotype (i.e. find hom alt sites) if available
        final List<PileupSummary> genotypingSites = matchedPileupSummariesTable == null ? sites :
                filterSitesByCoverage(PileupSummary.readFromFile(matchedPileupSummariesTable).getRight());

        final ContaminationModel genotypingModel = new ContaminationModel(genotypingSites, Optional.ofNullable(homSitesFile)); // sato: genotyping sites vs sites?

        if (outputTumorSegmentation != null) {
            final ContaminationModel tumorModel = matchedPileupSummariesTable == null ? genotypingModel :
                    new ContaminationModel(sites, Optional.ofNullable(homSitesFile));
            MinorAlleleFractionRecord.writeToFile(sample, tumorModel.segmentationRecords(), outputTumorSegmentation);
        }

        if (highCoverageSitesFile != null){
            PileupSummary.writeToFile("sample", sites, highCoverageSitesFile);
        }

        final Pair<Double, Double> contaminationAndError = genotypingModel.calculateContaminationFromHoms(sites);
        ContaminationRecord.writeToFile(Arrays.asList(
                new ContaminationRecord(sample, contaminationAndError.getLeft(), contaminationAndError.getRight())), outputTable);

        if (auxiliaryInfoFile != null){
            genotypingModel.writeMessages(auxiliaryInfoFile, true);
        }
        return "SUCCESS";
    }

    private List<PileupSummary> filterSitesByCoverage(final List<PileupSummary> allSites) {
        // Just in case the intervals given to GetPileupSummaries contained un-covered sites, we remove them
        // so that a bunch of zeroes don't throw off the median coverage
        final List<PileupSummary> coveredSites = allSites.stream().filter(s -> s.getTotalCount() > MIN_COVERAGE).collect(Collectors.toList());
        final double[] coverage = coveredSites.stream().mapToDouble(PileupSummary::getTotalCount).toArray();
        final double medianCoverage = new Median().evaluate(coverage);
        final double meanCoverage = new Mean().evaluate(coverage);
        final double lowCoverageThreshold = medianCoverage * lowCoverageRatioThreshold;
        final double highCoverageThreshold = meanCoverage * highCoverageRatioThreshold;
        return coveredSites.stream()
                .filter(ps -> ps.getTotalCount() > lowCoverageThreshold && ps.getTotalCount() < highCoverageThreshold)
                .collect(Collectors.toList());
    }

}
