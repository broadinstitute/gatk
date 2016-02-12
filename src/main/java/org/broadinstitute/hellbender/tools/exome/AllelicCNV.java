package org.broadinstitute.hellbender.tools.exome;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Detects copy-number events using allelic-count data and GATK CNV output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Detect copy-number events in a tumor sample using allelic-count data and GATK CNV output. " +
                "Allelic-count data (reference/alternate counts from the GetHetCoverage tool) is segmented using " +
                "circular binary segmentation; the result is combined with the target coverages " +
                "and called segments found by the GATK CNV tool. Bayesian parameter estimation of models for the " +
                "copy ratios and minor allele fractions in each segment is performed using Markov chain Monte Carlo.",
        oneLineSummary = "Detect copy-number events using allelic-count data and GATK CNV output.",
        programGroup = CopyNumberProgramGroup.class
)
public class AllelicCNV extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1l;

    //filename tags for output
    protected static final String SNP_MAF_SEG_FILE_TAG = "MAF";
    protected static final String UNION_SEG_FILE_TAG = "union";
    protected static final String SMALL_MERGED_SEG_FILE_TAG = "no-small";
    protected static final String INITIAL_SEG_FILE_TAG = "sim-0";
    protected static final String INTERMEDIATE_SEG_FILE_TAG = "sim";
    protected static final String FINAL_SEG_FILE_TAG = "sim-final";
    protected static final String GATK_SEG_FILE_TAG = "cnv";
    protected static final String CGA_ACS_SEG_FILE_TAG = "acs";

    private static final int MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS = 25;

    //CLI arguments
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

    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "intervalThresholdCopyRatio";
    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME = "simThCR";

    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "intervalThresholdAlleleFraction";
    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "simThAF";

    protected static final String USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME = "useAllCopyRatioSegments";
    protected static final String USE_ALL_COPY_RATIO_SEGMENTS_SHORT_NAME = "useAllCRSeg";

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Input file for tumor-sample tangent-normalized target coverages (.tn.tsv output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_SHORT_NAME,
            optional = false
    )
    protected File targetCoveragesFile;

    @Argument(
            doc = "Input file for tumor-sample target-coverage segments with calls (output of GATK CNV tool).",
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
    protected int numSamplesAlleleFraction = 200;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numBurnInAlleleFraction = 100;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for copy-ratio similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdCopyRatio = 2;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for allele-fraction similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdAlleleFraction = 2;

    @Argument(
            doc = "Enable use of all copy-ratio--segment breakpoints. " +
                    "(Default behavior uses only breakpoints from segments not called copy neutral, " +
                    "if calls are available in output of GATK CNV provided, and none otherwise.)",
            fullName = USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME,
            shortName = USE_ALL_COPY_RATIO_SEGMENTS_SHORT_NAME,
            optional = true
    )
    protected boolean useAllCopyRatioSegments = false;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        if (sampleName == null) {
            sampleName = outputPrefix;
        }

        logger.info("Starting workflow for sample " + sampleName + "...");

        //make Genome from input target coverages and SNP counts
        logger.info("Loading input files...");
        final Genome genome = new Genome(targetCoveragesFile, snpCountsFile, sampleName);

        //load target-coverage segments from input file
        final List<ModeledSegment> targetSegmentsWithCalls =
                SegmentUtils.readModeledSegmentsFromSegmentFile(targetSegmentsFile);
        logger.info("Number of copy-ratio segments from CNV output: " + targetSegmentsWithCalls.size());

        //merge copy-neutral and uncalled segments (unless disabled)
        final List<SimpleInterval> targetSegmentsUnfixed;
        if (!useAllCopyRatioSegments) {
            logger.info("Merging copy-neutral and uncalled segments...");
            targetSegmentsUnfixed = SegmentMergeUtils.mergeNeutralSegments(targetSegmentsWithCalls);
            logger.info("Number of segments after copy-neutral merging: " + targetSegmentsUnfixed.size());
        } else {
            targetSegmentsUnfixed = targetSegmentsWithCalls.stream().map(ModeledSegment::getSimpleInterval)
                    .collect(Collectors.toList());
        }

        //fix up target-segment start breakpoints (convert from target-end--target-end to target-start--target-end)
        final List<SimpleInterval> targetSegments =
                SegmentUtils.fixTargetSegmentStarts(targetSegmentsUnfixed, genome.getTargets());

        //segment SNPs on observed log_2 minor allele fraction
        logger.info("Performing SNP segmentation...");
        final File snpSegmentFile = new File(outputPrefix + "-" + SNP_MAF_SEG_FILE_TAG + ".seg");
        SNPSegmenter.writeSegmentFile(genome.getSNPs(), sampleName, snpSegmentFile);
        final List<SimpleInterval> snpSegments =
                SegmentUtils.readIntervalsFromSegmentFile(snpSegmentFile);
        logger.info("Number of SNP segments: " + snpSegments.size());

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

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller = new ACNVModeller(segmentedModel,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);
        final File initialModeledSegmentsFile = new File(outputPrefix + "-" + INITIAL_SEG_FILE_TAG + ".seg");
        modeller.writeACNVModeledSegmentFile(initialModeledSegmentsFile);

        //similar-segment merging (segment files are output for each merge iteration)
        logger.info("Merging similar segments...");
        logger.info("Initial number of segments before similar-segment merging: " + modeller.getACNVModeledSegments().size());
        List<ACNVModeledSegment> mergedSegments = new ArrayList<>(modeller.getACNVModeledSegments());
        //perform iterations of similar-segment merging until all similar segments are merged
        int prevNumSegments;
        for (int numIterations = 1; numIterations <= MAX_SIMILAR_SEGMENT_MERGE_ITERATIONS; numIterations++) {
            logger.info("Similar-segment merging iteration: " + numIterations);
            prevNumSegments = modeller.getACNVModeledSegments().size();
            modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction);
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
            final File modeledSegmentsFile = new File(outputPrefix + "-" + INTERMEDIATE_SEG_FILE_TAG + "-" + numIterations + ".seg");
            modeller.writeACNVModeledSegmentFile(modeledSegmentsFile);
        }
        logger.info("Final number of segments after similar-segment merging: " + modeller.getACNVModeledSegments().size());

        //write final model fit to file
        final File finalModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + ".seg");
        modeller.writeACNVModeledSegmentFile(finalModeledSegmentsFile);

        //write file for GATK CNV formatted seg file
        final File finalModeledSegmentsFileAsGatkCNV = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + "." + GATK_SEG_FILE_TAG + ".seg");
        modeller.writeModeledSegmentFile(finalModeledSegmentsFileAsGatkCNV);

        // Write file for ACS- compatible output to help Broad CGA
        final File finalACSModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + "." + CGA_ACS_SEG_FILE_TAG + ".seg");
        modeller.writeACNVModeledSegmentFileAsAllelicCapSegFile(finalACSModeledSegmentsFile);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: Allelic CNV run complete for sample " + sampleName + ".");
    }
}
