package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;

/**
 * Detects copy-number events using allelic-count data and ReCapSeg output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Detect copy-number events in a tumor sample using allelic-count data and ReCapSeg output. " +
                "Allelic-count data (reference/alternate counts from the GetHetCoverage tool) is segmented using " +
                "circular binary segmentation; the result is combined with the target coverages " +
                "and segments found by the ReCapSeg tool. Bayesian parameter estimation of models for the " +
                "copy ratios and minor allele fractions in each segment is performed using Markov chain Monte Carlo.",
        oneLineSummary = "Detect copy-number events using allelic-count data and ReCapSeg output.",
        programGroup = CopyNumberProgramGroup.class
)
public class AllelicCapSegger extends CommandLineProgram {
    protected static final String SNP_MAF_SEG_FILE_TAG = "MAF";
    protected static final String UNION_SEG_FILE_TAG = "union";
    protected static final String SMALL_MERGED_SEG_FILE_TAG = "no-small";

    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_LONG_NAME = "smallSegmentThreshold";
    protected static final String SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_SHORT_NAME = "smallTh";

    protected static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    protected static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    protected static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    protected static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    protected static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    protected static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    protected static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    protected static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    protected static final String SIGMA_THRESHOLD_COPY_RATIO_LONG_NAME = "sigmaThresholdCopyRatio";
    protected static final String SIGMA_THRESHOLD_COPY_RATIO_SHORT_NAME = "simThCR";

    protected static final String SIGMA_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "sigmaThresholdAlleleFraction";
    protected static final String SIGMA_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "simThAF";

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Input file for tumor-sample target coverages (output of ReCapSeg tool).",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_SHORT_NAME,
            optional = false
    )
    protected File targetCoveragesFile;

    @Argument(
            doc = "Input file for tumor-sample target-coverage segments (output of ReCapSeg tool).",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File targetSegmentsFile;

    @Argument(
            doc = "Prefix for output files. Will also be used as sample name if that is not provided.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Sample name. If not provided, prefix for output files will be used.",
            fullName = ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME,
            optional = true
    )
    protected String sampleName;

    @Argument(
            doc = "Threshold for small-segment merging. If a segment has strictly less than this number of targets, " +
                    "it is considered small and will be merged with an adjacent segment.",
            fullName = SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_LONG_NAME,
            shortName = SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_SHORT_NAME,
            optional = true
    )
    protected int smallSegmentTargetNumberThreshold = 3;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numSamplesAlleleFraction = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numBurnInAlleleFraction = 50;

    @Argument(
            doc = "Threshold number of standard deviations for copy-ratio similar-segment merging.",
            fullName = SIGMA_THRESHOLD_COPY_RATIO_LONG_NAME,
            shortName = SIGMA_THRESHOLD_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected double sigmaThresholdCopyRatio = 4;

    @Argument(
            doc = "Threshold number of standard deviations for allele-fraction similar-segment merging.",
            fullName = SIGMA_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = SIGMA_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double sigmaThresholdAlleleFraction = 2;

    @Override
    protected Object doWork() {
        if (sampleName == null) {
            sampleName = outputPrefix;
        }

        logger.info("Starting workflow for sample " + sampleName + "...");

        logger.info("Loading input files...");
        //make Genome from input target coverages and SNP counts
        final Genome genome = new Genome(targetCoveragesFile, snpCountsFile, sampleName);
        //load target-coverage segments from input file and fix up start breakpoints
        final List<SimpleInterval> targetSegmentsUnfixed =
                SegmentUtils.readIntervalsFromSegmentFile(targetSegmentsFile);
        final List<SimpleInterval> targetSegments =
                SegmentUtils.fixTargetSegmentStarts(targetSegmentsUnfixed, genome.getTargets());

        //segment SNPs on observed log_2 minor allele fraction
        logger.info("Performing SNP segmentation...");
        final File snpSegmentFile = new File(outputPrefix + "-" + SNP_MAF_SEG_FILE_TAG + ".seg");
        SNPSegmenter.writeSegmentFile(genome.getSNPs(), sampleName, snpSegmentFile);
        final List<SimpleInterval> snpSegments =
                SegmentUtils.readIntervalsFromSegmentFile(snpSegmentFile);

        //combine SNP and target-coverage segments
        logger.info("Combining SNP and target-coverage segments...");
        final List<SimpleInterval> unionedSegments =
                SegmentUtils.unionSegments(targetSegments, snpSegments, genome);
        final File unionedSegmentsFile = new File(outputPrefix + "-" + UNION_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeSegmentFileWithNumTargetsAndNumSNPs(unionedSegmentsFile, unionedSegments, genome);

        //small-segment merging (note that X and Y are always small segments and dropped, since GATK CNV drops them)
        logger.info("Merging small segments...");
        final SegmentedModel segmentedModelWithSmallSegments = new SegmentedModel(unionedSegments, genome);
        final SegmentedModel segmentedModel =
                segmentedModelWithSmallSegments.mergeSmallSegments(smallSegmentTargetNumberThreshold);
        final File segmentedModelFile = new File(outputPrefix + "-" + SMALL_MERGED_SEG_FILE_TAG + ".seg");
        segmentedModel.writeSegmentFileWithNumTargetsAndNumSNPs(segmentedModelFile);

        //initial model fitting
        logger.info("Fitting initial model...");
        final ACSModeller modeller = new ACSModeller(segmentedModel, outputPrefix,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction);

        //similar-segment merging (segment files are output for each merge iteration)
        logger.info("Merging similar segments...");
        modeller.mergeSimilarSegments(sigmaThresholdCopyRatio, sigmaThresholdAlleleFraction,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction);

        //TODO model-parameter output and plotting

        logger.info("AllelicCapSeg run complete for sample " + sampleName + ".");
        return "SUCCESS";
    }
}
