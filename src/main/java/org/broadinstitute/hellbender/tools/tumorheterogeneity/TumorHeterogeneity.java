package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.AllelicCNV;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Model tumor heterogeneity as a Dirichlet mixture of subclones with copy-number variation,
 * starting from {@link AllelicCNV} output.
 * (Alpha version: only a single tumor population is assumed.)
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model tumor heterogeneity as a Dirichlet mixture of subclones with copy-number variation, " +
                "starting from AllelicCNV output. (Alpha version: only a single tumor population is assumed.)",
        oneLineSummary = "Model tumor heterogeneity as a Dirichlet mixture of subclones with copy-number variation, " +
                "starting from AllelicCNV output",
        programGroup = CopyNumberProgramGroup.class
)
public final class TumorHeterogeneity extends SparkCommandLineProgram {
    private static final long serialVersionUID = 19738246L;

    private static final long RANDOM_SEED = 13;

    private static final int NUM_POPULATIONS_CLONAL = TumorHeterogeneityUtils.NUM_POPULATIONS_CLONAL;
    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);

    private static final double INITIAL_WALKER_BALL_SIZE_CLONAL = TumorHeterogeneityUtils.INITIAL_WALKER_BALL_SIZE_CLONAL;
    private static final double INITIAL_WALKER_BALL_SIZE = TumorHeterogeneityUtils.INITIAL_WALKER_BALL_SIZE;

    private static final double EPSILON = TumorHeterogeneityUtils.EPSILON;

    private static final double CONCENTRATION_PRIOR_ALPHA_CLONAL = TumorHeterogeneityUtils.CONCENTRATION_PRIOR_ALPHA_CLONAL;
    private static final double CONCENTRATION_PRIOR_BETA_CLONAL = TumorHeterogeneityUtils.CONCENTRATION_PRIOR_BETA_CLONAL;

    private static final double CONCENTRATION_MIN = TumorHeterogeneityUtils.CONCENTRATION_MIN;
    private static final double CONCENTRATION_MAX = TumorHeterogeneityUtils.CONCENTRATION_MAX;

    private static final double COPY_RATIO_NORMALIZATION_MIN = TumorHeterogeneityUtils.COPY_RATIO_NORMALIZATION_MIN;
    private static final double COPY_RATIO_NORMALIZATION_MAX = TumorHeterogeneityUtils.COPY_RATIO_NORMALIZATION_MAX;

    private static final double COPY_RATIO_NOISE_CONSTANT_MIN = TumorHeterogeneityUtils.COPY_RATIO_NOISE_CONSTANT_MIN;
    private static final double COPY_RATIO_NOISE_CONSTANT_MAX = TumorHeterogeneityUtils.COPY_RATIO_NOISE_CONSTANT_MAX;

    //filename tags for output
    protected static final String SAMPLES_FILE_SUFFIX_CLONAL = ".th.clonal.samples.tsv";
    protected static final String PROFILES_FILE_SUFFIX_CLONAL = ".th.clonal.profiles.tsv";
    protected static final String SUMMARY_FILE_SUFFIX_CLONAL = ".th.clonal.summary.tsv";
    protected static final String SAMPLES_FILE_SUFFIX = ".th.full.samples.tsv";
    protected static final String PROFILES_FILE_SUFFIX = ".th.full.profiles.tsv";
    protected static final String SUMMARY_FILE_SUFFIX = ".th.full.summary.tsv";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME = "maxAllelicCopyNumberClonal";
    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME = "maxACNClonal";

    protected static final String MAX_ALLELIC_COPY_NUMBER_LONG_NAME = "maxAllelicCopyNumber";
    protected static final String MAX_ALLELIC_COPY_NUMBER_SHORT_NAME = "maxACN";

    protected static final String MAX_NUM_POPULATIONS_LONG_NAME = "maxNumPopulations";
    protected static final String MAX_NUM_POPULATIONS_SHORT_NAME = "maxNumPop";

    protected static final String NUM_WALKERS_CLONAL_LONG_NAME = "numWalkersClonal";
    protected static final String NUM_WALKERS_CLONAL_SHORT_NAME = "numWalkClonal";

    protected static final String NUM_WALKERS_LONG_NAME = "numWalkers";
    protected static final String NUM_WALKERS_SHORT_NAME = "numWalk";

    protected static final String NUM_SAMPLES_CLONAL_LONG_NAME = "numSamplesClonal";
    protected static final String NUM_SAMPLES_CLONAL_SHORT_NAME = "numSampClonal";

    protected static final String NUM_SAMPLES_LONG_NAME = "numSamples";
    protected static final String NUM_SAMPLES_SHORT_NAME = "numSamp";

    protected static final String NUM_BURN_IN_CLONAL_LONG_NAME = "numBurnInClonal";
    protected static final String NUM_BURN_IN_CLONAL_SHORT_NAME = "numBurnClonal";

    protected static final String NUM_BURN_IN_LONG_NAME = "numBurnIn";
    protected static final String NUM_BURN_IN_SHORT_NAME = "numBurn";

    protected static final String CONCENTRATION_PRIOR_ALPHA_LONG_NAME = "concentrationPriorAlpha";
    protected static final String CONCENTRATION_PRIOR_ALPHA_SHORT_NAME = "concAlpha";

    protected static final String CONCENTRATION_PRIOR_BETA_LONG_NAME = "concentrationPriorBeta";
    protected static final String CONCENTRATION_PRIOR_BETA_SHORT_NAME = "concBeta";

    protected static final String COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME = "copyRatioNormalizationPriorAlpha";
    protected static final String COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_SHORT_NAME = "crNormAlpha";
    protected static final String COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME = "copyRatioNormalizationPriorBeta";
    protected static final String COPY_RATIO_NORMALIZATION_PRIOR_BETA_SHORT_NAME = "crNormBeta";

    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME = "copyRatioNoiseConstantPriorAlpha";
    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_SHORT_NAME = "crConstAlpha";
    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME = "copyRatioNoiseConstantPriorBeta";
    protected static final String COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_SHORT_NAME = "crConstBeta";

    protected static final String PLOIDY_MISMATCH_PENALTY_LONG_NAME = "ploidyMismatchPenalty";
    protected static final String PLOIDY_MISMATCH_PENALTY_SHORT_NAME = "ploidyPen";

    protected static final String SUBCLONE_VARIANCE_PENALTY_LONG_NAME = "subcloneVariancePenalty";
    protected static final String SUBCLONE_VARIANCE_PENALTY_SHORT_NAME = "subVarPen";

    protected static final String MODE_PURITY_BIN_SIZE_LONG_NAME = "purityBinSize";
    protected static final String MODE_PURITY_BIN_SIZE_SHORT_NAME = "purityBin";

    protected static final String MODE_PLOIDY_BIN_SIZE_LONG_NAME = "ploidyBinSize";
    protected static final String MODE_PLOIDY_BIN_SIZE_SHORT_NAME = "ploidyBin";

    @Argument(
            doc = "Input file for AllelicCNV result.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    protected File allelicCNVFile;

    @Argument(
            doc = "Prefix for output files.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME
    )
    protected String outputPrefix;

    @Argument(
            doc = "Maximum allelic copy number for clonal model.",
            fullName = MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME,
            shortName = MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int maxAllelicCopyNumberClonal = 5;

    @Argument(
            doc = "Maximum allelic copy number for full model.",
            fullName = MAX_ALLELIC_COPY_NUMBER_LONG_NAME,
            shortName = MAX_ALLELIC_COPY_NUMBER_SHORT_NAME,
            optional = true
    )
    protected int maxAllelicCopyNumber = 5;

    @Argument(
            doc = "Maximum number of populations for full model.",
            fullName = MAX_NUM_POPULATIONS_LONG_NAME,
            shortName = MAX_NUM_POPULATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumPopulations = 3;

    @Argument(
            doc = "Number of walkers in MCMC ensemble for clonal model.",
            fullName = NUM_WALKERS_CLONAL_LONG_NAME,
            shortName = NUM_WALKERS_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numWalkersClonal = 50;

    @Argument(
            doc = "Number of walkers in MCMC ensemble for full model.",
            fullName = NUM_WALKERS_LONG_NAME,
            shortName = NUM_WALKERS_SHORT_NAME,
            optional = true
    )
    protected int numWalkers = 100;

    @Argument(
            doc = "Total number of MCMC ensemble samples for clonal model. " +
                    "(Total number of samples will be number of walkers in ensemble multiplied by this number.)",
            fullName = NUM_SAMPLES_CLONAL_LONG_NAME,
            shortName = NUM_SAMPLES_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numSamplesClonal = 100;

    @Argument(
            doc = "Total number of MCMC ensemble samples for full model. " +
                    "(Total number of samples will be number of walkers in ensemble multiplied by this number.)",
            fullName = NUM_SAMPLES_LONG_NAME,
            shortName = NUM_SAMPLES_SHORT_NAME,
            optional = true
    )
    protected int numSamples = 200;

    @Argument(
            doc = "Number of burn-in ensemble samples to discard for clonal model. " +
                    "(Total number of samples will be number of walkers in ensemble multiplied by this number.)",
            fullName = NUM_BURN_IN_CLONAL_LONG_NAME,
            shortName = NUM_BURN_IN_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numBurnInClonal = 50;

    @Argument(
            doc = "Number of burn-in ensemble samples to discard for full model. " +
                    "(Total number of samples will be number of walkers in ensemble multiplied by this number.)",
            fullName = NUM_BURN_IN_LONG_NAME,
            shortName = NUM_BURN_IN_SHORT_NAME,
            optional = true
    )
    protected int numBurnIn = 100;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on concentration parameter.",
            fullName = CONCENTRATION_PRIOR_ALPHA_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorAlpha = 1.;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on concentration parameter.",
            fullName = CONCENTRATION_PRIOR_BETA_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorBeta = 1E1;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on copy-ratio normalization parameter.",
            fullName = COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME,
            shortName = COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNormalizationPriorAlpha = 1.;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on copy-ratio normalization parameter.",
            fullName = COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME,
            shortName = COPY_RATIO_NORMALIZATION_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNormalizationPriorBeta = 1.;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on copy-ratio noise-constant parameter.",
            fullName = COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseConstantPriorAlpha = 1;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on copy-ratio noise-constant parameter.",
            fullName = COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseConstantPriorBeta = 1E2;

    @Argument(
            doc = "Penalty factor for ploidy mismatch in proposal of variant profiles. " +
                    "(A strong (i.e., large) penalty factor prefers correct solutions at the expense of increased mixing time.)",
            fullName = PLOIDY_MISMATCH_PENALTY_LONG_NAME,
            shortName = PLOIDY_MISMATCH_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double ploidyMismatchPenalty = 1E3;

    @Argument(
            doc = "Penalty factor for number of subclones with unique ploidy states per base. " +
                    "(A strong (i.e., large) penalty factor prefers solutions with fewer unique ploidy states per segment.)",
            fullName = SUBCLONE_VARIANCE_PENALTY_LONG_NAME,
            shortName = SUBCLONE_VARIANCE_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double subcloneVariancePenalty = 1E3;

    @Argument(
            doc = "Purity bin size for identifying samples at posterior mode.",
            fullName = MODE_PURITY_BIN_SIZE_LONG_NAME,
            shortName = MODE_PURITY_BIN_SIZE_SHORT_NAME,
            optional = true
    )
    protected double purityModeBinSize = 0.025;

    @Argument(
            doc = "Ploidy bin size for identifying samples at posterior mode.",
            fullName = MODE_PLOIDY_BIN_SIZE_LONG_NAME,
            shortName = MODE_PLOIDY_BIN_SIZE_SHORT_NAME,
            optional = true
    )
    protected double ploidyModeBinSize = 0.025;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        //load ACNV segments from input file
        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(allelicCNVFile);

        //initialize output files
        final File samplesFileClonal = new File(outputPrefix + SAMPLES_FILE_SUFFIX_CLONAL);
        final File profilesFileClonal = new File(outputPrefix + PROFILES_FILE_SUFFIX_CLONAL);
        final File summaryFileClonal = new File(outputPrefix + SUMMARY_FILE_SUFFIX_CLONAL);

        //construct priors using input parameters
        final PloidyStatePrior ploidyStatePriorClonal = calculatePloidyStatePrior(maxAllelicCopyNumberClonal);
        final TumorHeterogeneityPriorCollection priorsClonal = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, ploidyStatePriorClonal,
                CONCENTRATION_PRIOR_ALPHA_CLONAL, CONCENTRATION_PRIOR_BETA_CLONAL,
                copyRatioNormalizationPriorAlpha, copyRatioNormalizationPriorBeta,
                copyRatioNoiseConstantPriorAlpha, copyRatioNoiseConstantPriorBeta,
                ploidyMismatchPenalty, subcloneVariancePenalty);
        
        //initialize data collection from ACNV input and priors
        final TumorHeterogeneityData dataClonal = new TumorHeterogeneityData(segments, priorsClonal);

        //initialize modeller and run MCMC
        final TumorHeterogeneityModeller modellerClonal =
                new TumorHeterogeneityModeller(dataClonal, NUM_POPULATIONS_CLONAL, numWalkersClonal, INITIAL_WALKER_BALL_SIZE_CLONAL, rng);
        modellerClonal.fitMCMC(numSamplesClonal, numBurnInClonal);

        //initialize writer
        final TumorHeterogeneityModellerWriter writerClonal = new TumorHeterogeneityModellerWriter(modellerClonal);

        //write all MCMC samples to file
        writerClonal.writePopulationFractionAndPloidySamples(samplesFileClonal);
        logger.info("Clonal run: MCMC samples output to " + samplesFileClonal + ".");

        //identify samples in purity-ploidy bin centered on posterior mode
        final TumorHeterogeneityState posteriorModeClonal = modellerClonal.getPosteriorMode();
        posteriorModeClonal.values().forEach(p -> logger.info("Clonal run: Posterior mode " + p.getName().name() + ": " + p.getValue()));
        final List<Integer> indicesOfSamplesAtModeClonal = modellerClonal.collectIndicesOfSamplesInBin(posteriorModeClonal, purityModeBinSize, ploidyModeBinSize);

        //average variant profiles of identified samples and write to file
        logger.info("Clonal run: Calculating averaged variant profiles at posterior mode...");
        writerClonal.writeAveragedProfiles(profilesFileClonal, indicesOfSamplesAtModeClonal);
        logger.info("Clonal run: Averaged variant profiles at posterior mode output to " + profilesFileClonal + ".");

        //calculate posterior summaries for global parameters from identified samples and write to file
        logger.info("Clonal run: Calculating summary for posterior mode...");
        writerClonal.writePosteriorSummaries(summaryFileClonal, indicesOfSamplesAtModeClonal, ctx);
        logger.info("Clonal run: Calculating summary for posterior-mode variant profiles output to " + summaryFileClonal + ".");
        
        if (numSamples > 0) {
            //initialize output files
            final File samplesFile = new File(outputPrefix + SAMPLES_FILE_SUFFIX);
            final File profilesFile = new File(outputPrefix + PROFILES_FILE_SUFFIX);
            final File summaryFile = new File(outputPrefix + SUMMARY_FILE_SUFFIX);

            //construct priors using input parameters
            final PloidyStatePrior ploidyStatePrior = calculatePloidyStatePrior(maxAllelicCopyNumber);
            final TumorHeterogeneityPriorCollection priors = new TumorHeterogeneityPriorCollection(
                    NORMAL_PLOIDY_STATE, ploidyStatePrior,
                    concentrationPriorAlpha, concentrationPriorBeta,
                    copyRatioNormalizationPriorAlpha, copyRatioNormalizationPriorBeta,
                    copyRatioNoiseConstantPriorAlpha, copyRatioNoiseConstantPriorBeta,
                    ploidyMismatchPenalty, subcloneVariancePenalty);

            //initialize data collection from using clonal data and full priors
            final TumorHeterogeneityData data = new TumorHeterogeneityData(dataClonal, priors);
            
            //initialize modeller and run MCMC
            logger.info("Full run: Initializing MCMC from clonal result...");
            final TumorHeterogeneityState initialState = TumorHeterogeneityState.initializeFromClonalState(priors, posteriorModeClonal, maxNumPopulations);
            final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(data, initialState, numWalkers, INITIAL_WALKER_BALL_SIZE, rng);
            modeller.fitMCMC(numSamples, numBurnIn);

            //initialize writer
            final TumorHeterogeneityModellerWriter writer = new TumorHeterogeneityModellerWriter(modeller);

            //write all MCMC samples to file
            writer.writePopulationFractionAndPloidySamples(samplesFile);
            logger.info("Full run: MCMC samples output to " + samplesFile + ".");

            //identify samples in purity-ploidy bin centered on posterior mode
            final TumorHeterogeneityState posteriorMode = modeller.getPosteriorMode();
            posteriorMode.values().forEach(p -> logger.info("Full run: Posterior mode " + p.getName().name() + ": " + p.getValue()));
            final List<Integer> indicesOfSamplesAtMode = modeller.collectIndicesOfSamplesInBin(posteriorMode, purityModeBinSize, ploidyModeBinSize);

            //average variant profiles of identified samples and write to file
            logger.info("Full run: Calculating averaged variant profiles at posterior mode...");
            writer.writeAveragedProfiles(profilesFile, indicesOfSamplesAtMode);
            logger.info("Full run: Averaged variant profiles at posterior mode output to " + profilesFile + ".");

            //calculate posterior summaries for global parameters from identified samples and write to file
            logger.info("Full run: Calculating summary for posterior mode...");
            writer.writePosteriorSummaries(summaryFile, indicesOfSamplesAtMode, ctx);
            logger.info("Full run: Calculating summary for posterior-mode variant profiles output to " + summaryFile + ".");
        }

        logger.info("SUCCESS: TumorHeterogeneity run complete.");
    }

    //we implement a simple prior on ploidy states that penalizes copy-number changes, with an additional penalty for homozygous deletions
    private static PloidyStatePrior calculatePloidyStatePrior(final int maxAllelicCopyNumber) {
        final Function<PloidyState, Double> ploidyLogPDF = ps -> {
            if (ps.equals(NORMAL_PLOIDY_STATE)) {
                return Math.log(50.);
            }
//            if (ps.total() % 2 == 0) {
//                return Math.log(2.);
//            }
            return Math.log(1.);
        };
//                Math.log(Math.max(EPSILON, 1. - ploidyStatePriorCompleteDeletionPenalty)) * (ps.m() == 0 && ps.n() == 0 ? 1 : 0)
//                        + Math.log(Math.max(EPSILON, 1. - ploidyStatePriorChangePenalty)) * (Math.abs(NORMAL_PLOIDY_STATE.m() - ps.m()) + Math.abs(NORMAL_PLOIDY_STATE.n() - ps.n()));
        final Map<PloidyState, Double> ploidyStateToUnnormalizedLogProbabilityMap = new LinkedHashMap<>();
        for (int m = 0; m <= maxAllelicCopyNumber; m++) {
            for (int n = 0; n <= maxAllelicCopyNumber; n++) {
                ploidyStateToUnnormalizedLogProbabilityMap.put(new PloidyState(m, n), ploidyLogPDF.apply(new PloidyState(m, n)));
            }
        }
        return new PloidyStatePrior(ploidyStateToUnnormalizedLogProbabilityMap);
    }

    //validate CLI arguments
    private void validateArguments() {
        Utils.regularReadableUserFile(allelicCNVFile);
        Utils.validateArg(maxAllelicCopyNumberClonal > 0, MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(maxAllelicCopyNumber > 0, MAX_ALLELIC_COPY_NUMBER_LONG_NAME + " must be positive.");
        Utils.validateArg(maxNumPopulations > NUM_POPULATIONS_CLONAL, MAX_NUM_POPULATIONS_LONG_NAME + " should be strictly greater than " + NUM_POPULATIONS_CLONAL + ".");
        Utils.validateArg(numSamplesClonal > 0, NUM_SAMPLES_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(numSamples >= 0, NUM_SAMPLES_LONG_NAME + " must be non-negative.");
        Utils.validateArg(numBurnInClonal >= 0 && numBurnInClonal < numSamplesClonal, NUM_BURN_IN_CLONAL_LONG_NAME + " must be non-negative and strictly less than " + NUM_SAMPLES_CLONAL_LONG_NAME + ".");
        Utils.validateArg(numSamples == 0 ? numBurnIn == 0 : numBurnIn >= 0 && numBurnIn < numSamples, NUM_BURN_IN_LONG_NAME + " must be non-negative and less than or equal to " + NUM_SAMPLES_LONG_NAME + ".");
        validatePriorHyperparameters(
                concentrationPriorAlpha, CONCENTRATION_PRIOR_ALPHA_LONG_NAME,
                concentrationPriorBeta, CONCENTRATION_PRIOR_BETA_LONG_NAME,
                CONCENTRATION_MIN, CONCENTRATION_MAX,
                (alpha, beta) -> alpha / beta);
        validatePriorHyperparameters(
                copyRatioNormalizationPriorAlpha, COPY_RATIO_NORMALIZATION_PRIOR_ALPHA_LONG_NAME,
                copyRatioNormalizationPriorBeta, COPY_RATIO_NORMALIZATION_PRIOR_BETA_LONG_NAME,
                COPY_RATIO_NORMALIZATION_MIN, COPY_RATIO_NORMALIZATION_MAX,
                (alpha, beta) -> alpha / beta);
        validatePriorHyperparameters(
                copyRatioNoiseConstantPriorAlpha, COPY_RATIO_NOISE_CONSTANT_PRIOR_ALPHA_LONG_NAME,
                copyRatioNoiseConstantPriorBeta, COPY_RATIO_NOISE_CONSTANT_PRIOR_BETA_LONG_NAME,
                COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX,
                (alpha, beta) -> alpha / beta);
        Utils.validateArg(ploidyMismatchPenalty >= 0., PLOIDY_MISMATCH_PENALTY_LONG_NAME + " must be non-negative.");
        Utils.validateArg(subcloneVariancePenalty >= 0., SUBCLONE_VARIANCE_PENALTY_LONG_NAME + " must be non-negative.");
        Utils.validateArg(0. <= purityModeBinSize && purityModeBinSize <= 1., "Invalid purity bin size for determining mode.");
        Utils.validateArg(0. <= ploidyModeBinSize && ploidyModeBinSize <= 2 * maxAllelicCopyNumberClonal, "Invalid ploidy bin size for determining mode.");
    }

    //validate hyperparameters for Beta or Gamma distributions
    private static void validatePriorHyperparameters(final double alpha, final String alphaName,
                                                     final double beta, final String betaName,
                                                     final double min, final double max,
                                                     final BiFunction<Double, Double, Double> calculateMeanFromHyperparameters) {
        Utils.validateArg(alpha > 0, alphaName + " must be positive.");
        Utils.validateArg(beta > 0, betaName + " must be positive.");
        final double boundedQuantity = calculateMeanFromHyperparameters.apply(alpha, beta);
        Utils.validateArg(min < boundedQuantity && boundedQuantity < max, "Mean calculated from hyperparameters " + alphaName + " and " + betaName + " must be in (" + min + ", " + max + ").");
    }
}