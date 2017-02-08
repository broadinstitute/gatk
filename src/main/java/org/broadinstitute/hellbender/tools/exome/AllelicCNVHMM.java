package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.segmentation.AFCRHiddenState;
import org.broadinstitute.hellbender.tools.exome.segmentation.JointAFCRSegmenter;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Detects copy-number events using allelic-count data and GATK CNV output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "(EXPERIMENTAL) Detect copy-number events in a tumor sample using allelic-count data and GATK CNV output.",
        oneLineSummary = "(EXPERIMENTAL) Detect copy-number events using allelic-count data and GATK CNV output",
        programGroup = CopyNumberProgramGroup.class
)
public class AllelicCNVHMM extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1l;

    //filename tags for output
    protected static final String INITIAL_FIT_FILE_TAG = "sim-begin";
    protected static final String FINAL_FIT_FILE_TAG = "sim-final";
    protected static final String SEGMENT_FILE_SUFFIX = ".seg";
    protected static final String CR_PARAMETER_FILE_SUFFIX = ".cr.param";
    protected static final String AF_PARAMETER_FILE_SUFFIX = ".af.param";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String INITIAL_NUM_COPY_RATIO_STATES_LONG_NAME = "initialNumberOfCopyRatioStates";
    protected static final String INITIAL_NUM_COPY_RATIO_STATES_SHORT_NAME = "initialNumCRStates";

    protected static final String INITIAL_NUM_ALLELE_FRACTION_STATES_LONG_NAME = "initialNumberOfAlleleFractionStates";
    protected static final String INITIAL_NUM_ALLELE_FRACTION_STATES_SHORT_NAME = "initialNumAFStates";

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

    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME = "maxNumIterationsSimSeg";
    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME = "maxIterSim";

    protected static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME = "numIterationsSimSegPerFit";
    protected static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME = "numIterSimPerFit";

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage or GetBayesianHetCoverage tools).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Input file for tumor-sample tangent-normalized target log_2 coverages (.tn.tsv output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File tangentNormalizedCoverageFile;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPoNFile;

    @Argument(
            doc = "Prefix for output files. Will also be used as the sample name in downstream plots." +
                    "(Note: if this is a file path or contains slashes (/), " +
                    "the string after the final slash will be used as the sample name in downstream plots.)",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Initial number of hidden copy-ratio states",
            fullName = INITIAL_NUM_COPY_RATIO_STATES_LONG_NAME,
            shortName = INITIAL_NUM_COPY_RATIO_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumCRStates = 10;

    @Argument(
            doc = "Initial number of hidden allele-fraction states",
            fullName = INITIAL_NUM_ALLELE_FRACTION_STATES_LONG_NAME,
            shortName = INITIAL_NUM_ALLELE_FRACTION_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumAFStates = 10;

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
    protected double intervalThresholdCopyRatio = 2.;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for allele-fraction similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for similar-segment merging.",
            fullName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumSimilarSegmentMergingIterations = 25;

    @Argument(
            doc = "Number of similar-segment--merging iterations per MCMC model refit. " +
                    "(Increasing this will decrease runtime, but the final number of segments may be higher. " +
                    "Setting this to 0 will completely disable model refitting between iterations.)",
            fullName = NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME,
            shortName = NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME,
            optional = true
    )
    protected int numSimilarSegmentMergingIterationsPerFit = 5;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        //the string after the final slash in the output prefix (which may be an absolute file path) will be used as the sample name
        final String sampleName = outputPrefix.substring(outputPrefix.lastIndexOf("/") + 1);

        logger.info("Starting workflow for " + sampleName + "...");

        //make Genome from input target coverages and SNP counts
        logger.info("Loading input files...");
        final Genome genome = new Genome(tangentNormalizedCoverageFile, snpCountsFile);

        //load allelic-bias panel of normals if provided
        final AllelicPanelOfNormals allelicPoN =
                allelicPoNFile != null ? AllelicPanelOfNormals.read(allelicPoNFile) : AllelicPanelOfNormals.EMPTY_PON;

        //load target-coverage segments from input file
        final AllelicCountCollection acc = new AllelicCountCollection(snpCountsFile);
        final ReadCountCollection rcc;
        try {
            rcc = ReadCountCollectionUtils.parse(tangentNormalizedCoverageFile);
        } catch (final IOException ex) {
            throw new UserException.BadInput("could not read input file");
        }
        final JointAFCRSegmenter jointSegmenter = JointAFCRSegmenter.createJointSegmenter(initialNumCRStates, rcc, initialNumAFStates, acc);
        final List<Pair<SimpleInterval, AFCRHiddenState>> segmentation = jointSegmenter.findSegments();
        final List<SimpleInterval> segments = segmentation.stream().map(Pair::getLeft).collect(Collectors.toList());
        final SegmentedGenome segmentedGenome = new SegmentedGenome(segments, genome);
        logger.info(String.format("Joint segmentation resulted in %d segments.", segments.size()));

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller = new ACNVModeller(segmentedGenome, allelicPoN,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);

        //write initial segments and parameters to file
        writeACNVModeledSegmentAndParameterFiles(modeller, INITIAL_FIT_FILE_TAG);

        //similar-segment merging (segment files are output for each merge iteration)
        logger.info("Merging similar segments...");
        performSimilarSegmentMergingStep(modeller);

        //write final segments and parameters to file
        writeACNVModeledSegmentAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: Allelic CNV run complete for sample " + sampleName + ".");
    }

    //validate CLI arguments
    private void validateArguments() {
        Utils.validateArg(numSamplesCopyRatio > 0, NUM_SAMPLES_COPY_RATIO_LONG_NAME + " must be positive.");
        Utils.validateArg(numSamplesCopyRatio > numBurnInCopyRatio, NUM_SAMPLES_COPY_RATIO_LONG_NAME + " must be greater than " + NUM_BURN_IN_COPY_RATIO_LONG_NAME);
        Utils.validateArg(numSamplesAlleleFraction > 0, NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME + " must be positive.");
        Utils.validateArg(numSamplesAlleleFraction > numBurnInAlleleFraction, NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME + " must be greater than " + NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME);
        Utils.validateArg(intervalThresholdCopyRatio > 0, INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME + " must be positive.");
        Utils.validateArg(intervalThresholdAlleleFraction > 0, INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME + " must be positive.");
        Utils.validateArg(maxNumSimilarSegmentMergingIterations >= 0, MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME + " must be non-negative.");
        Utils.validateArg(numSimilarSegmentMergingIterationsPerFit >= 0, NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME + " must be non-negative.");
    }

    //similar-segment merging
    private void performSimilarSegmentMergingStep(final ACNVModeller modeller) {
        logger.info("Initial number of segments before similar-segment merging: " + modeller.getACNVModeledSegments().size());
        //perform iterations of similar-segment merging until all similar segments are merged
        for (int numIterations = 1; numIterations <= maxNumSimilarSegmentMergingIterations; numIterations++) {
            logger.info("Similar-segment merging iteration: " + numIterations);
            final int prevNumSegments = modeller.getACNVModeledSegments().size();
            if (numSimilarSegmentMergingIterationsPerFit > 0 && numIterations % numSimilarSegmentMergingIterationsPerFit == 0) {
                //refit model after this merge iteration
                modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction, true);
            } else {
                //do not refit model after this merge iteration (deciles will be unspecified)
                modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction, false);
            }
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
        }
        if (!modeller.isModelFit()) {
            //make sure final model is completely fit (i.e., deciles are specified)
            modeller.fitModel();
        }
        logger.info("Final number of segments after similar-segment merging: " + modeller.getACNVModeledSegments().size());
    }

    //write modeled segments and global parameters to file
    private void writeACNVModeledSegmentAndParameterFiles(final ACNVModeller modeller, final String fileTag) {
        final File modeledSegmentsFile = new File(outputPrefix + "-" + fileTag + SEGMENT_FILE_SUFFIX);
        modeller.writeACNVModeledSegmentFile(modeledSegmentsFile);
        logSegmentFileWrittenMessage(modeledSegmentsFile);
        final File copyRatioParameterFile = new File(outputPrefix + "-" + fileTag + CR_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(outputPrefix + "-" + fileTag + AF_PARAMETER_FILE_SUFFIX);
        modeller.writeACNVModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
    }

    private void logSegmentFileWrittenMessage(final File file) {
        logger.info("Segments written to file " + file);
    }
}