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
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.Function;

/**
 * Detects copy-number events using allelic-count data and GATK CNV output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model tumor heterogeneity using a Dirichlet mixture of subclones with copy number variation.",
        oneLineSummary = "Model tumor heterogeneity using a Dirichlet mixture of subclones with copy number variation",
        programGroup = CopyNumberProgramGroup.class
)
public class TumorHeterogeneity extends SparkCommandLineProgram {
    private static final long serialVersionUID = 19738246L;

    private static final long RANDOM_SEED = 13;
    private static final int NUM_POPULATIONS_CLONAL = 2;
    private static final PloidyState NORMAL_PLOIDY_STATE = new PloidyState(1, 1);

    //filename tags for output
    protected static final String CLONAL_SAMPLES_FILE_SUFFIX = ".th.clonal.samples.tsv";
    protected static final String CLONAL_SUMMARY_FILE_SUFFIX = ".th.clonal.summary.tsv";
    protected static final String RESULT_FILE_SUFFIX = ".th.final.tsv";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME = "maxAllelicCopyNumberClonal";
    protected static final String MAX_ALLELIC_COPY_NUMBER_CLONAL_SHORT_NAME = "maxACNClonal";

    protected static final String MAX_ALLELIC_COPY_NUMBER_LONG_NAME = "maxAllelicCopyNumber";
    protected static final String MAX_ALLELIC_COPY_NUMBER_SHORT_NAME = "maxACN";

    protected static final String MAX_NUM_POPULATIONS_LONG_NAME = "maxNumPopulations";
    protected static final String MAX_NUM_POPULATIONS_SHORT_NAME = "maxNumPop";

    protected static final String NUM_SAMPLES_CLONAL_LONG_NAME = "numSamplesClonal";
    protected static final String NUM_SAMPLES_CLONAL_SHORT_NAME = "numSampClonal";

    protected static final String NUM_BURN_IN_CLONAL_LONG_NAME = "numBurnInClonal";
    protected static final String NUM_BURN_IN_CLONAL_SHORT_NAME = "numBurnClonal";

    protected static final String NUM_SAMPLES_LONG_NAME = "numSamples";
    protected static final String NUM_SAMPLES_SHORT_NAME = "numSamp";

    protected static final String NUM_BURN_IN_LONG_NAME = "numBurnIn";
    protected static final String NUM_BURN_IN_SHORT_NAME = "numBurn";

    protected static final String CONCENTRATION_PRIOR_ALPHA_CLONAL_LONG_NAME = "concentrationPriorAlphaClonal";
    protected static final String CONCENTRATION_PRIOR_ALPHA_CLONAL_SHORT_NAME = "concAlphaClonal";

    protected static final String CONCENTRATION_PRIOR_BETA_CLONAL_LONG_NAME = "concentrationPriorBetaClonal";
    protected static final String CONCENTRATION_PRIOR_BETA_CLONAL_SHORT_NAME = "concBetaClonal";

    protected static final String CONCENTRATION_PRIOR_ALPHA_LONG_NAME = "concentrationPriorAlpha";
    protected static final String CONCENTRATION_PRIOR_ALPHA_SHORT_NAME = "concAlpha";

    protected static final String CONCENTRATION_PRIOR_BETA_LONG_NAME = "concentrationPriorBeta";
    protected static final String CONCENTRATION_PRIOR_BETA_SHORT_NAME = "concBeta";

    protected static final String COPY_RATIO_NOISE_FLOOR_PRIOR_ALPHA_LONG_NAME = "copyRatioNoiseFloorPriorAlpha";
    protected static final String COPY_RATIO_NOISE_FLOOR_PRIOR_ALPHA_SHORT_NAME = "crFloorAlpha";

    protected static final String COPY_RATIO_NOISE_FLOOR_PRIOR_BETA_LONG_NAME = "copyRatioNoiseFloorPriorBeta";
    protected static final String COPY_RATIO_NOISE_FLOOR_PRIOR_BETA_SHORT_NAME = "crFloorBeta";

    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME = "copyRatioNoiseFactorPriorAlpha";
    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME = "crNoiseAlpha";

    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_LONG_NAME = "copyRatioNoiseFactorPriorBeta";
    protected static final String COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME = "crNoiseBeta";

    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME = "minorAlleleFractionNoiseFactorPriorAlpha";
    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME = "mafNoiseAlpha";

    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_LONG_NAME = "minorAlleleFractionNoiseFactorPriorBeta";
    protected static final String MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME = "mafNoiseBeta";

    protected static final String PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME = "completeDeletionPenalty";
    protected static final String PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_SHORT_NAME = "compDelPen";

    protected static final String PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME = "changePenalty";
    protected static final String PLOIDY_STATE_PRIOR_CHANGE_PENALTY_SHORT_NAME = "changePen";

    @Argument(
            doc = "Input file for AllelicCNV result.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File allelicCNVFile;

    @Argument(
            doc = "Prefix for output files.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
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
            doc = "Total number of MCMC samples for clonal model.",
            fullName = NUM_SAMPLES_CLONAL_LONG_NAME,
            shortName = NUM_SAMPLES_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numSamplesClonal = 2000;

    @Argument(
            doc = "Number of burn-in samples to discard for clonal model.",
            fullName = NUM_BURN_IN_CLONAL_LONG_NAME,
            shortName = NUM_BURN_IN_CLONAL_SHORT_NAME,
            optional = true
    )
    protected int numBurnInClonal = 500;

    @Argument(
            doc = "Total number of MCMC samples for full model.",
            fullName = NUM_SAMPLES_LONG_NAME,
            shortName = NUM_SAMPLES_SHORT_NAME,
            optional = true
    )
    protected int numSamples = 5000;

    @Argument(
            doc = "Number of burn-in samples to discard for full model.",
            fullName = NUM_BURN_IN_LONG_NAME,
            shortName = NUM_BURN_IN_SHORT_NAME,
            optional = true
    )
    protected int numBurnIn = 1000;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on concentration parameter for clonal model.",
            fullName = CONCENTRATION_PRIOR_ALPHA_CLONAL_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_ALPHA_CLONAL_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorAlphaClonal = 1E1;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on concentration parameter for clonal model.",
            fullName = CONCENTRATION_PRIOR_BETA_CLONAL_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_BETA_CLONAL_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorBetaClonal = 1E1;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on concentration parameter for full model.",
            fullName = CONCENTRATION_PRIOR_ALPHA_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorAlpha = 1.;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on concentration parameter for full model.",
            fullName = CONCENTRATION_PRIOR_BETA_LONG_NAME,
            shortName = CONCENTRATION_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double concentrationPriorBeta = 1E2;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on copy-ratio noise-floor parameter.",
            fullName = COPY_RATIO_NOISE_FLOOR_PRIOR_ALPHA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_FLOOR_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseFloorPriorAlpha = 1;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on copy-ratio noise-floor parameter.",
            fullName = COPY_RATIO_NOISE_FLOOR_PRIOR_BETA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_FLOOR_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseFloorPriorBeta = 1E2;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on copy-ratio noise-factor parameter.",
            fullName = COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseFactorPriorAlpha = 1;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on copy-ratio noise-factor parameter.",
            fullName = COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_LONG_NAME,
            shortName = COPY_RATIO_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double copyRatioNoiseFactorPriorBeta = 1E3;

    @Argument(
            doc = "Alpha hyperparameter for Gamma-distribution prior on minor-allele-fraction noise-factor parameter.",
            fullName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_LONG_NAME,
            shortName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_ALPHA_SHORT_NAME,
            optional = true
    )
    protected double minorAlleleFractionNoiseFactorPriorAlpha = 1.;

    @Argument(
            doc = "Beta hyperparameter for Gamma-distribution prior on minor-allele-fraction noise-factor parameter.",
            fullName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_LONG_NAME,
            shortName = MINOR_ALLELE_FRACTION_NOISE_FACTOR_PRIOR_BETA_SHORT_NAME,
            optional = true
    )
    protected double minorAlleleFractionNoiseFactorPriorBeta = 1E3;

    @Argument(
        doc = "Penalty for complete allele deletion in ploidy-state prior.",
        fullName = PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME,
        shortName = PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_SHORT_NAME,
        optional = true
    )
    protected double ploidyStatePriorCompleteDeletionPenalty = 1E-1;

    @Argument(
            doc = "Penalty for copy change in ploidy-state prior.",
            fullName = PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME,
            shortName = PLOIDY_STATE_PRIOR_CHANGE_PENALTY_SHORT_NAME,
            optional = true
    )
    protected double ploidyStatePriorChangePenalty = 1E-1;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(allelicCNVFile);

        final TumorHeterogeneityData data = new TumorHeterogeneityData(segments);

        final PloidyStatePrior ploidyStatePriorClonal = calculatePloidyStatePrior(ploidyStatePriorCompleteDeletionPenalty, ploidyStatePriorChangePenalty, maxAllelicCopyNumberClonal);
        final TumorHeterogeneityPriorCollection priorsClonal = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, ploidyStatePriorClonal,
                concentrationPriorAlphaClonal, concentrationPriorBetaClonal,
                copyRatioNoiseFloorPriorAlpha, copyRatioNoiseFloorPriorBeta,
                copyRatioNoiseFactorPriorAlpha, copyRatioNoiseFactorPriorBeta,
                minorAlleleFractionNoiseFactorPriorAlpha, minorAlleleFractionNoiseFactorPriorBeta);

        final PloidyStatePrior ploidyStatePrior = calculatePloidyStatePrior(ploidyStatePriorCompleteDeletionPenalty, ploidyStatePriorChangePenalty, maxAllelicCopyNumber);
        final TumorHeterogeneityPriorCollection priors = new TumorHeterogeneityPriorCollection(
                NORMAL_PLOIDY_STATE, ploidyStatePrior,
                concentrationPriorAlpha, concentrationPriorBeta,
                copyRatioNoiseFloorPriorAlpha, copyRatioNoiseFloorPriorBeta,
                copyRatioNoiseFactorPriorAlpha, copyRatioNoiseFactorPriorBeta,
                minorAlleleFractionNoiseFactorPriorAlpha, minorAlleleFractionNoiseFactorPriorBeta);

        final File samplesFileClonal = new File(outputPrefix + CLONAL_SAMPLES_FILE_SUFFIX);
        final File summaryFileClonal = new File(outputPrefix + CLONAL_SUMMARY_FILE_SUFFIX);
        final TumorHeterogeneityModeller clonalModeller = new TumorHeterogeneityModeller(data, priorsClonal, NUM_POPULATIONS_CLONAL, rng);
        clonalModeller.fitMCMC(numSamplesClonal, numBurnInClonal);
        clonalModeller.outputSamples(samplesFileClonal);
        clonalModeller.outputSummary(summaryFileClonal);

        logger.info("Tumor heterogeneity clonal run complete.");
        logger.info("Result samples output to " + samplesFileClonal + ".");
        logger.info("Result summary output to " + summaryFileClonal + ".");

//        final TumorHeterogeneityState initialState = TumorHeterogeneityState.initializeStateFromClonalResult(priors, clonalModeller, maxNumPopulations);
//        final File resultFile = new File(outputPrefix + RESULT_FILE_SUFFIX);
//        final TumorHeterogeneityModeller modeller = new TumorHeterogeneityModeller(data, initialState, rng);
//        modeller.fitMCMC(numSamples, numBurnIn);
//        modeller.output(resultFile);
//
//        logger.info("SUCCESS: Tumor heterogeneity full run complete and result output to " + resultFile + ".");
    }

    private static PloidyStatePrior calculatePloidyStatePrior(final double ploidyStatePriorCompleteDeletionPenalty,
                                                              final double ploidyStatePriorChangePenalty,
                                                              final int maxAllelicCopyNumber) {
        final Function<PloidyState, Double> ploidyLogPDF = ps -> Math.log(ploidyStatePriorCompleteDeletionPenalty) * ((ps.m() == 0 ? 1 : 0) + (ps.n() == 0 ? 1 : 0))
                + Math.log(ploidyStatePriorChangePenalty) * (Math.abs(NORMAL_PLOIDY_STATE.m() - ps.m()) + Math.abs(NORMAL_PLOIDY_STATE.n() - ps.n()));
        final Map<PloidyState, Double> unnormalizedLogProbabilityMassFunctionMap = new LinkedHashMap<>();
        for (int n = 0; n <= maxAllelicCopyNumber; n++) {
            for (int m = 0; m <= n; m++) {
                unnormalizedLogProbabilityMassFunctionMap.put(new PloidyState(m, n), ploidyLogPDF.apply(new PloidyState(m, n)));
            }
        }
        return new PloidyStatePrior(unnormalizedLogProbabilityMassFunctionMap);
    }

    //validate CLI arguments
    private void validateArguments() {
        //TODO UPDATE THIS
        Utils.validateArg(maxAllelicCopyNumberClonal > 0, MAX_ALLELIC_COPY_NUMBER_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(maxAllelicCopyNumber > 0, MAX_ALLELIC_COPY_NUMBER_LONG_NAME + " must be positive.");
        Utils.validateArg(maxNumPopulations >= 2, MAX_NUM_POPULATIONS_LONG_NAME + " must be greater than or equal to 2.");
        Utils.validateArg(numSamplesClonal > 0, NUM_SAMPLES_CLONAL_LONG_NAME + " must be positive.");
        Utils.validateArg(numBurnInClonal >= 0 && numBurnInClonal < numSamplesClonal, NUM_BURN_IN_CLONAL_LONG_NAME + " must be non-negative and strictly less than " + NUM_SAMPLES_CLONAL_LONG_NAME);
        Utils.validateArg(numSamples > 0, NUM_SAMPLES_LONG_NAME + " must be positive.");
        Utils.validateArg(numBurnIn >= 0 && numBurnIn < numSamples, NUM_BURN_IN_LONG_NAME + " must be non-negative and strictly less than " + NUM_SAMPLES_LONG_NAME);
        Utils.validateArg(concentrationPriorAlpha > 0, CONCENTRATION_PRIOR_ALPHA_LONG_NAME + " must be positive.");
        Utils.validateArg(concentrationPriorBeta > 0, CONCENTRATION_PRIOR_BETA_LONG_NAME + " must be positive.");
        Utils.validateArg(TumorHeterogeneityModeller.CONCENTRATION_MIN < concentrationPriorAlpha / concentrationPriorBeta &&
                concentrationPriorAlpha / concentrationPriorBeta < TumorHeterogeneityModeller.CONCENTRATION_MAX,
                CONCENTRATION_PRIOR_ALPHA_LONG_NAME + " / " + CONCENTRATION_PRIOR_BETA_LONG_NAME + " must be in (" +
                        TumorHeterogeneityModeller.CONCENTRATION_MIN + ", " + TumorHeterogeneityModeller.CONCENTRATION_MAX + ").");
        Utils.validateArg(ploidyStatePriorCompleteDeletionPenalty > 0, PLOIDY_STATE_PRIOR_COMPLETE_DELETION_PENALTY_LONG_NAME + " must be positive.");
        Utils.validateArg(ploidyStatePriorChangePenalty > 0, PLOIDY_STATE_PRIOR_CHANGE_PENALTY_LONG_NAME + " must be positive.");
    }
}