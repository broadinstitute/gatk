package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang3.Range;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
 * <p>This tool and GetPileupSummaries together replace GATK3's ContEst. Like ContEst, this tool estimates contamination based on the signal
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

    private static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    private static final int MAX_CHANGEPOINTS_PER_CHROMOSOME = 10;
    private static final int MIN_SITES_PER_SEGMENT = 5;

    // our analysis only cares about hom alt and het sites, so we throw away hom refs with a very conservative heuristic
    private static final double ALT_FRACTION_OF_DEFINITE_HOM_REF = 0.05;

    private static final double STRICT_LOH_MAF_THRESHOLD = 0.4;

    private static final double INITIAL_CONTAMINATION_GUESS = 0.05;
    private static final int MAX_ITERATIONS = 10;
    private static final double CONTAMINATION_CONVERGENCE_THRESHOLD = 0.001;

    private static final Range<Double> ALT_FRACTIONS_FOR_SEGMENTATION = Range.between(0.1, 0.9);

    private static final double KERNEL_SEGMENTER_LINEAR_COST = 1.0;
    private static final double KERNEL_SEGMENTER_LOG_LINEAR_COST = 1.0;
    private static final int KERNEL_SEGMENTER_DIMENSION = 100;
    private static final int POINTS_PER_SEGMENTATION_WINDOW = 50;

    private static final int MIN_COVERAGE = 10;
    private static final double DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD = 1.0/2;
    private static final double DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD = 3.0;
    public static final int DESIRED_MINIMUM_HOM_ALT_COUNT = 50;
    public static final double MINOR_ALLELE_FRACTION_STEP_SIZE = 0.05;

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

    private static final double SEGMENTATION_KERNEL_VARIANCE = 0.025;

    private static final BiFunction<PileupSummary, PileupSummary, Double> SEGMENTATION_KERNEL = (ps1, ps2) -> {
        final double maf1 = FastMath.min(ps1.getAltFraction(), 1 - ps1.getAltFraction());
        final double maf2 = FastMath.min(ps2.getAltFraction(), 1 - ps2.getAltFraction());
        return FastMath.exp(-MathUtils.square(maf1 - maf2)/(2 * SEGMENTATION_KERNEL_VARIANCE));
    };

    @Override
    public Object doWork() {
        final List<PileupSummary> sites = filterSites(PileupSummary.readFromFile(inputPileupSummariesTable));

        // used the matched normal to genotype (i.e. find hom alt sites) if available
        final List<PileupSummary> genotypingSites = matchedPileupSummariesTable == null ? sites :
                filterSites(PileupSummary.readFromFile(matchedPileupSummariesTable));

        // we partition the genome into contiguous allelic copy-number segments in order to infer the local minor
        // allele fraction at each site.  This is important because a minor allele fraction close to 1/2 (neutral)
        // allows hets and hom alts to be distinguished easily, while a low minor allele fraction makes it harder
        // to discriminate.  It is crucial to know which site are true hom alts and which sites are hets with
        // loss of heterozygosity.  We do this for the genotyping sample because that is the sample from which
        // the hom alts are deduced.
        final List<List<PileupSummary>> genotypingSegments = findSegments(genotypingSites);


        List<PileupSummary> homAltGenotypingSites = new ArrayList<>();
        final MutableDouble genotypingContamination = new MutableDouble(INITIAL_CONTAMINATION_GUESS);

        for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
            List<List<PileupSummary>> homAltSitesBySegment = Arrays.asList(new ArrayList<>());
            final MutableDouble minorAlleleFractionThreshold = new MutableDouble(STRICT_LOH_MAF_THRESHOLD);
            while (homAltSitesBySegment.stream().mapToInt(List::size).sum() < DESIRED_MINIMUM_HOM_ALT_COUNT && minorAlleleFractionThreshold.doubleValue() > 0) {
                homAltSitesBySegment = genotypingSegments.stream()
                        .map(segment -> segmentHomAlts(segment, genotypingContamination.doubleValue(), minorAlleleFractionThreshold.doubleValue()))
                        .collect(Collectors.toList());
                minorAlleleFractionThreshold.subtract(MINOR_ALLELE_FRACTION_STEP_SIZE);
            }
            homAltGenotypingSites = homAltSitesBySegment.stream().flatMap(List::stream).collect(Collectors.toList());
            final double newGenotypingContamination = calculateContamination(homAltGenotypingSites, errorRate(genotypingSites)).getLeft();
            if (Math.abs(newGenotypingContamination - genotypingContamination.doubleValue()) < CONTAMINATION_CONVERGENCE_THRESHOLD) {
                break;
            }
            genotypingContamination.setValue(newGenotypingContamination);
        }

        if (outputTumorSegmentation != null) {
            final List<List<PileupSummary>> tumorSegments = matchedPileupSummariesTable == null ?
                    genotypingSegments : findSegments(sites);
            List<MinorAlleleFractionRecord> tumorMinorAlleleFractions = tumorSegments.stream()
                    .map(this::makeMinorAlleleFractionRecord).collect(Collectors.toList());
            MinorAlleleFractionRecord.writeToFile(tumorMinorAlleleFractions, outputTumorSegmentation);

        }

        final List<PileupSummary> homAltSites = subsetSites(sites, homAltGenotypingSites);
        final Pair<Double, Double> contaminationAndError = calculateContamination(homAltSites, errorRate(sites));
        final double contamination = contaminationAndError.getLeft();
        final double error = contaminationAndError.getRight();
        ContaminationRecord.writeToFile(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), contamination, error)), outputTable);

        return "SUCCESS";
    }

    private List<List<PileupSummary>> findSegments(final List<PileupSummary> sites) {
        final Map<String, List<PileupSummary>> sitesByContig = sites.stream().collect(Collectors.groupingBy(PileupSummary::getContig));

        return sitesByContig.values().stream()
                .flatMap(contig -> findContigSegments(contig).stream())
                .filter(segment -> segment.size() >= MIN_SITES_PER_SEGMENT)
                .collect(Collectors.toList());
    }

    // in a biallelic site, essentially every non-ref, non-primary alt base is an error, since there are 2 such possible
    // errors out of 3 total, we multiply by 3/2 to get the total base error rate
    private double errorRate(List<PileupSummary> sites) {
        final long totalBases = sites.stream().mapToInt(PileupSummary::getTotalCount).sum();
        final long otherAltBases = sites.stream().mapToInt(PileupSummary::getOtherAltCount).sum();
        return 1.5 * ((double) otherAltBases / totalBases);
    }

    // subset sites in the contaminated sample to hom alt site found in the genotyping sample
    private static List<PileupSummary> subsetSites(final List<PileupSummary> sites, final List<PileupSummary> subsetLoci) {
        final OverlapDetector<PileupSummary> homAltsInMatchedNormalOverlapDetector = OverlapDetector.create(subsetLoci);
        return sites.stream().filter(homAltsInMatchedNormalOverlapDetector::overlapsAny).collect(Collectors.toList());
    }

    private List<PileupSummary> segmentHomAlts(final List<PileupSummary> segment, final double contamination, double minimiumMinorAlleleFraction) {
        final double minorAlleleFraction = calculateMinorAlleleFraction(segment);
        return minorAlleleFraction < minimiumMinorAlleleFraction ? Collections.emptyList() :
                segment.stream().filter(site -> homAltProbability(site, minorAlleleFraction, contamination) > 0.5).collect(Collectors.toList());
    }

    private double calculateMinorAlleleFraction(final List<PileupSummary> segment) {
        final List<PileupSummary> hets = getLikelyHetsBasedOnAlleleFraction(segment);
        final Function<Double, Double> objective = maf -> logLikelihoodOfHetsInSegment(hets, maf);
        return OptimizationUtils.argmax(objective, ALT_FRACTIONS_FOR_SEGMENTATION.getMinimum(), 0.5, 0.4, 0.01, 0.01, 20);
    }

    private MinorAlleleFractionRecord makeMinorAlleleFractionRecord(final List<PileupSummary> segment) {
        final String contig = segment.get(0).getContig();
        final int start = segment.get(0).getStart();
        final int end = segment.get(segment.size() - 1).getEnd();
        final double minorAlleleFraction = calculateMinorAlleleFraction(segment);
        return new MinorAlleleFractionRecord(new SimpleInterval(contig, start, end), minorAlleleFraction);
    }


    // we want log(1/2 (likelihood of alt minor + likelihood of alt major))
    //         =  logSumLog(log likelihood of alt minor, log likelihood of alt major) - log(2)
    private final double logLikelihoodOfHetsInSegment(final List<PileupSummary> hets, final double minorAlleleFraction) {
        return hets.stream().mapToDouble(het -> {
            final int n = het.getTotalCount();
            final int a = het.getAltCount();
            final double altMinorLogLikelihood = new BinomialDistribution(null, n, minorAlleleFraction).logProbability(a);
            final double altMajorLogLikelihood = new BinomialDistribution(null, n, 1 - minorAlleleFraction).logProbability(a);

            return MathUtils.logSumLog(altMinorLogLikelihood, altMajorLogLikelihood) + MathUtils.LOG_ONE_HALF;
        }).sum();
    }

    private List<List<PileupSummary>> findContigSegments(List<PileupSummary> sites) {
        // segment based on obvious hets
        final List<PileupSummary> hetSites = getLikelyHetsBasedOnAlleleFraction(sites);

        if (hetSites.isEmpty()) {
            return Collections.emptyList();
        }

        final List<Integer> changepoints = new ArrayList<>();
        // when the kernel segmenter finds a changepoint at index n, that means index n belongs to the left segment, which goes
        // against the usual end-exclusive intervals of IndexRange etc.  This explains adding in the first changepoint of -1
        // instead of 0 and all the "changepoint + 1" constructions below
        changepoints.add(-1);
        final KernelSegmenter<PileupSummary> segmenter = new KernelSegmenter<>(hetSites);
        changepoints.addAll(segmenter.findChangepoints(MAX_CHANGEPOINTS_PER_CHROMOSOME, SEGMENTATION_KERNEL, KERNEL_SEGMENTER_DIMENSION,
                Arrays.asList(POINTS_PER_SEGMENTATION_WINDOW), KERNEL_SEGMENTER_LINEAR_COST, KERNEL_SEGMENTER_LOG_LINEAR_COST, KernelSegmenter.ChangepointSortOrder.INDEX));
        changepoints.add(hetSites.size()-1);

        final List<SimpleInterval> segments = IntStream.range(0, changepoints.size() - 1)
                .mapToObj(n -> {
                    final PileupSummary firstSiteInSegment = hetSites.get(changepoints.get(n) + 1);
                    final PileupSummary lastSiteInSegment = hetSites.get(changepoints.get(n+1));
                    return new SimpleInterval(firstSiteInSegment.getContig(), firstSiteInSegment.getStart(), lastSiteInSegment.getEnd());
                }).collect(Collectors.toList());

        final OverlapDetector<PileupSummary> od = OverlapDetector.create(sites);

        // for each segment, find overlapping sites and sort by coordinate
        return segments.stream()
                .map(segment -> od.getOverlaps(segment).stream().sorted(Comparator.comparingInt(PileupSummary::getStart)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }

    private List<PileupSummary> getLikelyHetsBasedOnAlleleFraction(List<PileupSummary> sites) {
        return sites.stream()
                    .filter(ps -> ALT_FRACTIONS_FOR_SEGMENTATION.contains(ps.getAltFraction()))
                    .collect(Collectors.toList());
    }

    private static Pair<Double, Double> calculateContamination(List<PileupSummary> homAltSites, final double errorRate) {
        if (homAltSites.isEmpty()) {
            logger.warn("No hom alt sites found!  Perhaps GetPileupSummaries was run on too small of an interval, or perhaps the sample was extremely inbred or haploid.");
            return Pair.of(0.0, 1.0);
        }

        final long totalReadCount = homAltSites.stream().mapToLong(PileupSummary::getTotalCount).sum();
        final long totalRefCount = homAltSites.stream().mapToLong(PileupSummary::getRefCount).sum();

        // if eg ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final long errorRefCount = Math.round(totalReadCount * errorRate / 3);
        final long contaminationRefCount = Math.max(totalRefCount - errorRefCount, 0);
        final double totalDepthWeightedByRefFrequency = homAltSites.stream()
                .mapToDouble(ps -> ps.getTotalCount() * (1 - ps.getAlleleFrequency()))
                .sum();
        final double contamination = contaminationRefCount / totalDepthWeightedByRefFrequency;
        final double standardError = Math.sqrt(contamination / totalDepthWeightedByRefFrequency);

        logger.info(String.format("In %d homozygous variant sites we find %d reference reads due to contamination and %d" +
                    " due to to sequencing error out of a total %d reads.", homAltSites.size(), contaminationRefCount, errorRefCount, totalReadCount));
        logger.info(String.format("Based on population data, we would expect %d reference reads in a contaminant with equal depths at these sites.", (long) totalDepthWeightedByRefFrequency));
        logger.info(String.format("Therefore, we estimate a contamination of %.3f.", contamination));
        logger.info(String.format("The error bars on this estimate are %.5f.", standardError));

        return Pair.of(Math.min(contamination, 1.0), standardError);
    }

    private List<PileupSummary> filterSites(final List<PileupSummary> allSites) {
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
                .filter(ps -> ps.getAltFraction() > ALT_FRACTION_OF_DEFINITE_HOM_REF)
                .collect(Collectors.toList());
    }

    private double homAltProbability(final PileupSummary site, final double minorAlleleFraction, final double contamination) {
        final double alleleFrequency = site.getAlleleFrequency();
        final double homAltPrior = MathUtils.square(alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);

        final int altCount = site.getAltCount();
        final int totalCount = altCount + site.getRefCount();

        if (altCount < totalCount / 2) {
            return 0;
        }

        final double homAltLikelihood = new BinomialDistribution(null, totalCount, 1 - contamination).probability(altCount);
        final double hetLikelihood = new BinomialDistribution(null, totalCount, 1 - minorAlleleFraction).probability(altCount);

        final double unnormalizedHomAltProbability = homAltPrior * homAltLikelihood;
        final double unnormalizedHetProbability = hetPrior * hetLikelihood;

        final double result = unnormalizedHomAltProbability / (unnormalizedHetProbability + unnormalizedHomAltProbability);

        return result;

    }
}
