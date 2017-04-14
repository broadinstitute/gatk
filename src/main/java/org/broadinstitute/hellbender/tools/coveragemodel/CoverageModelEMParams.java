package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;

/**
 * This class represents the hyper-parameters and parameters for {@link CoverageModelEMAlgorithm}
 *
 * Works both as an argument collection as well as a getter/setter.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelEMParams {

    public enum PsiUpdateMode {
        /**
         * Target-resolved unexplained variance
         */
        PSI_TARGET_RESOLVED,

        /**
         * Isotropic unexplained variance (i.e. one value for all targets)
         */
        PSI_ISOTROPIC
    }

    public enum WSolverStrategy {
        /**
         * Perform the M-step update of bias covariates locally
         */
        W_SOLVER_LOCAL,

        /**
         * Perform the M-step update of bias covariates using Spark
         */
        W_SOLVER_SPARK
    }

    public enum CopyRatioHMMType {
        /**
         * Update the copy ratio (or copy number) locally
         */
        COPY_RATIO_HMM_LOCAL,

        /**
         * Update the copy ratio (or copy number) using Spark
         */
        COPY_RATIO_HMM_SPARK
    }

    public enum ComputeNodeCommunicationPolicy {
        /**
         * Push new data to the compute nodes via broadcast-hash-join method
         */
        BROADCAST_HASH_JOIN,

        /**
         * Make an RDD from the new data and perform a joint-map with the compute note RDDs
         */
        RDD_JOIN
    }

    public enum ModelInitializationStrategy {
        /**
         * Random factor loadings, zero mean bias, zero target-specific unexplained variance
         */
        RANDOM,

        /**
         * Estimate factor loadings, mean bias, and target-specific unexplained variance by performing PCA
         */
        PCA
    }

    public static final ModelInitializationStrategy DEFAULT_MODEL_INITIALIZATION_STRATEGY = ModelInitializationStrategy.PCA;
    public static final String MODEL_INITIALIZATION_STRATEGY_SHORT_NAME = "MIS";
    public static final String MODEL_INITIALIZATION_STRATEGY_LONG_NAME = "modelInitializationStrategy";

    public static final double DEFAULT_LOG_LIKELIHOOD_TOL = 1e-5;
    public static final String LOG_LIKELIHOOD_TOL_SHORT_NAME = "LLT";
    public static final String LOG_LIKELIHOOD_TOL_LONG_NAME = "logLikelihoodTol";

    public static final double DEFAULT_PARAM_ABS_TOL = 1e-5;
    public static final String PARAM_ABS_TOL_SHORT_NAME = "PMAT";
    public static final String PARAM_ABS_TOL_LONG_NAME = "paramAbsoluteTolerance";

    public static final double DEFAULT_POSTERIOR_ABS_TOL = 1e-2;
    public static final String POSTERIOR_ABS_TOL_SHORT_NAME = "PSAT";
    public static final String POSTERIOR_ABS_TOL_LONG_NAME = "posteriorAbsoluteTolerance";

    public static final double DEFAULT_MEAN_FIELD_ADMIXING_RATIO = 1.0;
    public static final String MEAN_FIELD_ADMIXING_RATIO_SHORT_NAME = "MFAR";
    public static final String MEAN_FIELD_ADMIXING_RATIO_LONG_NAME = "meanFieldAdmixingRatio";

    public static final int DEFAULT_MAX_M_STEP_CYCLES = 1;
    public static final String MAX_M_STEP_CYCLES_SHORT_NAME = "MMC";
    public static final String MAX_M_STEP_CYCLES_LONG_NAME = "maximumMStepCycles";

    public static final int DEFAULT_MAX_E_STEP_CYCLES = 1;
    public static final String MAX_E_STEP_CYCLES_SHORT_NAME = "MES";
    public static final String MAX_E_STEP_CYCLES_LONG_NAME = "maximumEStepCycles";

    public static final int DEFAULT_MAX_EM_ITERATIONS = 50;
    public static final String MAX_EM_ITERATIONS_SHORT_NAME = "MEMI";
    public static final String MAX_EM_ITERATIONS_LONG_NAME = "maximumEMIterations";

    public static final double DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING = 1e-2;
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING_SHORT_NAME = "LLTTCRC";
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING_LONG_NAME = "logLikelihoodTolThresholdCopyRatioCalling";

    public static final double DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_PSI_SWITCHING = 1e-2;
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_PSI_SWITCHING_SHORT_NAME = "LLTTPS";
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_PSI_SWITCHING_LONG_NAME = "logLikelihoodTolThresholdPsiSwitching";

    public static final double DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_ARD_UPDATE = 1e-2;
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_ARD_UPDATE_SHORT_NAME = "LLTTARDU";
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_ARD_UPDATE_LONG_NAME = "logLikelihoodTolThresholdARDUpdate";

    public static final int DEFAULT_NUM_LATENTS = 5;
    public static final String NUM_LATENTS_SHORT_NAME = "NL";
    public static final String NUM_LATENTS_LONG_NAME = "numLatents";

    public static final int DEFAULT_MIN_LEARNING_READ_COUNT = 5;
    public static final String MIN_LEARNING_READ_COUNT_SHORT_NAME = "MLRC";
    public static final String MIN_LEARNING_READ_COUNT_LONG_NAME = "minimumLearningReadCount";

    public static final int DEFAULT_MIN_PCA_INIT_READ_COUNT = 10;
    public static final String MIN_PCA_INIT_READ_COUNT_SHORT_NAME = "MPCAIRC";
    public static final String MIN_PCA_INIT_READ_COUNT_LONG_NAME = "minimumPCAInitializationReadCount";

    public static final double DEFAULT_MAPPING_ERROR_RATE = 1e-3;
    public static final String MAPPING_ERROR_RATE_SHORT_NAME = "MER";
    public static final String MAPPING_ERROR_RATE_LONG_NAME = "mappingErrorRate";

    /* ARD related */

    public static final boolean DEFAULT_ARD_ENABLED = true;
    public static final String ARD_ENABLED_SHORT_NAME = "ARDU";
    public static final String ARD_ENABLED_LONG_NAME = "ARDUpdate";

    public static final double DEFAULT_INITIAL_ARD_PRECISION_ABSOLUTE = 1e-8;
    public static final String INITIAL_ARD_PRECISION_ABSOLUTE_SHORT_NAME = "IARDPA";
    public static final String INITIAL_ARD_PRECISION_ABSOLUTE_LONG_NAME = "initialARDPrecisionAbsolute";

    public static final double DEFAULT_INITIAL_ARD_PRECISION_RELATIVE_TO_NOISE = 1e-3;
    public static final String INITIAL_ARD_PRECISION_RELATIVE_TO_NOISE_SHORT_NAME = "IARDPRN";
    public static final String INITIAL_ARD_PRECISION_RELATIVE_TO_NOISE_LONG_NAME = "initialARDPrecisionRelativeToNoise";

    public static final double DEFAULT_MAX_ARD_PRECISION = 1e100;
    public static final String MAX_ARD_PRECISION_SHORT_NAME = "MARDP";
    public static final String MAX_ARD_PRECISION_LONG_NAME = "maxARDPrecision";

    public static final boolean DEFAULT_INCLUDE_ARD_IN_LOG_LIKELIHOOD = false;
    public static final String INCLUDE_ARD_IN_LOG_LIKELIHOOD_SHORT_NAME = "IARDILL";
    public static final String INCLUDE_ARD_IN_LOG_LIKELIHOOD_LONG_NAME = "includeARDInLogLikelihood";

    /* Psi related */

    public static final boolean DEFAULT_PSI_UPDATE_ENABLED = true;
    public static final String PSI_UPDATE_ENABLED_SHORT_NAME = "PU";
    public static final String PSI_UPDATE_ENABLED_LONG_NAME = "psiUpdate";

    public static final PsiUpdateMode DEFAULT_PSI_SOLVER_MODE = PsiUpdateMode.PSI_TARGET_RESOLVED;
    public static final String PSI_SOLVER_MODE_SHORT_NAME = "PSM";
    public static final String PSI_SOLVER_MODE_LONG_NAME = "psiUpdateMode";

    public static final boolean DEFAULT_ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED = false;
    public static final String ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_SHORT_NAME = "APSMS";
    public static final String ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_LONG_NAME = "adaptivePsiSolverModeSwitching";

    public static final double DEFAULT_PSI_UPPER_LIMIT = 5.0;
    public static final String PSI_UPPER_LIMIT_SHORT_NAME = "PUL";
    public static final String PSI_UPPER_LIMIT_LONG_NAME = "psiUpperLimit";

    public static final int DEFAULT_MAX_PSI_ITERATIONS = 200;
    public static final String MAX_PSI_ITERATIONS_SHORT_NAME = "MPI";
    public static final String MAX_PSI_ITERATIONS_LONG_NAME = "maxPsiIterations";

    public static final double DEFAULT_PSI_ABS_TOL = 1e-8;
    public static final String PSI_ABS_TOL_SHORT_NAME = "PAT";
    public static final String PSI_ABS_TOL_LONG_NAME = "psiAbsoluteTolerance";

    public static final double DEFAULT_PSI_REL_TOL = 1e-6;
    public static final String PSI_REL_TOL_SHORT_NAME = "PRT";
    public static final String PSI_REL_TOL_LONG_NAME = "psiRelativeTolerance";

    public static final int DEFAULT_PSI_SOLVER_NUM_BISECTIONS = 10;
    public static final String PSI_SOLVER_NUM_BISECTIONS_SHORT_NAME = "PSNB";
    public static final String PSI_SOLVER_NUM_BISECTIONS_LONG_NAME = "psiSolverNumBisections";

    public static final int DEFAULT_PSI_SOLVER_REFINEMENT_DEPTH = 3;
    public static final String PSI_SOLVER_REFINEMENT_DEPTH_SHORT_NAME = "SSRD";
    public static final String PSI_SOLVER_REFINEMENT_DEPTH_LONG_NAME = "psiSolverRefinementDepth";

    public static final int DEFAULT_PSI_SOLVER_NUM_THREADS = 1;
    public static final String PSI_SOLVER_NUM_THREADS_SHORT_NAME = "PSNT";
    public static final String PSI_SOLVER_NUM_THREADS_LONG_NAME = "psiSolverNumThreads";


    /* W related */

    public static final WSolverStrategy DEFAULT_W_SOLVER_TYPE = WSolverStrategy.W_SOLVER_SPARK;
    public static final String W_SOLVER_TYPE_SHORT_NAME = "WST";
    public static final String W_SOLVER_TYPE_LONG_NAME = "wSolverStrategy";

    public static final ComputeNodeCommunicationPolicy DEFAULT_W_COMMUNICATION_POLICY = ComputeNodeCommunicationPolicy.RDD_JOIN;
    public static final String W_COMMUNICATION_POLICY_SHORT_NAME = "WCP";
    public static final String W_COMMUNICATION_POLICY_LONG_NAME = "wCommunicationPolicy";

    public static final int DEFAULT_MAX_W_ITERATIONS = 200;
    public static final String MAX_W_ITERATIONS_SHORT_NAME = "MWI";
    public static final String MAX_W_ITERATIONS_LONG_NAME = "maxWInterations";

    public static final double DEFAULT_W_ABS_TOL = 1e-8;
    public static final String W_ABS_TOL_SHORT_NAME = "WAT";
    public static final String W_ABS_TOL_LONG_NAME = "wAbsoluteTolerance";

    public static final double DEFAULT_W_REL_TOL = 1e-6;
    public static final String W_REL_TOL_SHORT_NAME = "WRT";
    public static final String W_REL_TOL_LONG_NAME = "wRelativeTolerance";

    /* gamma related */

    public static final boolean DEFAULT_GAMMA_UPDATE_ENABLED = true;
    public static final String GAMMA_UPDATE_ENABLED_SHORT_NAME = "GU";
    public static final String GAMMA_UPDATE_ENABLED_LONG_NAME = "gammaUpdate";

    public static final double DEFAULT_GAMMA_UPPER_LIMIT = 1.0;
    public static final String GAMMA_UPPER_LIMIT_SHORT_NAME = "GUL";
    public static final String GAMMA_UPPER_LIMIT_LONG_NAME = "gammaUpperLimit";

    public static final int DEFAULT_MAX_GAMMA_ITERATIONS = 200;
    public static final String MAX_GAMMA_ITERATIONS_SHORT_NAME = "MGI";
    public static final String MAX_GAMMA_ITERATIONS_LONG_NAME = "maxGammaIterations";

    public static final double DEFAULT_GAMMA_ABS_TOL = 1e-8;
    public static final String GAMMA_ABS_TOL_SHORT_NAME = "GAT";
    public static final String GAMMA_ABS_TOL_LONG_NAME = "gammaAbsoluteTolerance";

    public static final double DEFAULT_GAMMA_REL_TOL = 1e-6;
    public static final String GAMMA_REL_TOL_SHORT_NAME = "GRT";
    public static final String GAMMA_REL_TOL_LONG_NAME = "gammaRelativeTolerance";

    public static final int DEFAULT_GAMMA_SOLVER_NUM_BISECTIONS = 10;
    public static final String GAMMA_SOLVER_NUM_BISECTIONS_SHORT_NAME = "GSNB";
    public static final String GAMMA_SOLVER_NUM_BISECTIONS_LONG_NAME = "gammaSolverNumBisections";

    public static final int DEFAULT_GAMMA_SOLVER_REFINEMENT_DEPTH = 3;
    public static final String GAMMA_SOLVER_REFINEMENT_DEPTH_SHORT_NAME = "GSRD";
    public static final String GAMMA_SOLVER_REFINEMENT_DEPTH_LONG_NAME = "gammaSolverRefinementDepth";

    /* CNV-avoiding penalty related */

    public static final boolean DEFAULT_FOURIER_REGULARIZATION_ENABLED = false;
    public static final String FOURIER_REGULARIZATION_ENABLED_SHORT_NAME = "FR";
    public static final String FOURIER_REGULARIZATION_ENABLED_LONG_NAME = "fourierRegularization";

    public static final int DEFAULT_MIN_CNV_LENGTH = 10;
    public static final String MIN_CNV_LENGTH_SHORT_NAME = "MNCL";
    public static final String MIN_CNV_LENGTH_LONG_NAME = "minimumCNVLength";

    public static final int DEFAULT_MAX_CNV_LENGTH = 1000;
    public static final String MAX_CNV_LENGTH_SHORT_NAME = "MXCL";
    public static final String MAX_CNV_LENGTH_LONG_NAME = "maximumCNVLength";

    public static final boolean DEFAULT_ZERO_PAD_FFT = false;
    public static final String ZERO_PAD_FFT_SHORT_NAME = "ZPFFT";
    public static final String ZERO_PAD_FFT_LONG_NAME = "zeroPadFFT";

    public static final double DEFAULT_FOURIER_REGULARIZATION_STRENGTH = 10_000;
    public static final String FOURIER_REGULARIZATION_STRENGTH_SHORT_NAME = "FRS";
    public static final String FOURIER_REGULARIZATION_STRENGTH_LONG_NAME = "fourierRegularizationStrength";


    public static final boolean DEFAULT_CR_UPDATE_ENABLED = true;
    public static final String CR_UPDATE_ENABLED_SHORT_NAME = "CRU";
    public static final String CR_UPDATE_ENABLED_LONG_NAME = "copyRatioUpdate";

    public static final CopyRatioHMMType DEFAULT_CR_HMM_TYPE = CopyRatioHMMType.COPY_RATIO_HMM_SPARK;
    public static final String CR_HMM_TYPE_SHORT_NAME = "CRHMM";
    public static final String CR_HMM_TYPE_LONG_NAME = "copyRatioHMMType";


    public static final int DEFAULT_NUMBER_OF_TARGET_SPACE_PARTITIONS = 1;
    public static final String NUMBER_OF_TARGET_SPACE_PARTITIONS_SHORT_NAME = "NTSP";
    public static final String NUMBER_OF_TARGET_SPACE_PARTITIONS_LONG_NAME = "numTargetSpacePartitions";

    public static final int DEFAULT_RDD_CHECKPOINTING_INTERVAL = 10;
    public static final String RDD_CHECKPOINTING_INTERVAL_SHORT_NAME = "RDDCPI";
    public static final String RDD_CHECKPOINTING_INTERVAL_LONG_NAME = "rddCheckpointingInterval";

    public static final boolean DEFAULT_RDD_CHECKPOINTING_ENABLED = true;
    public static final String RDD_CHECKPOINTING_ENABLED_SHORT_NAME = "RDDCP";
    public static final String RDD_CHECKPOINTING_ENABLED_LONG_NAME = "rddCheckpointing";

    public static final String DEFAULT_RDD_CHECKPOINTING_PATH = "/dev/null";
    public static final String RDD_CHECKPOINTING_PATH_SHORT_NAME = "RDDCPP";
    public static final String RDD_CHECKPOINTING_PATH_LONG_NAME = "rddCheckpointingPath";


    public static final int DEFAULT_RUN_CHECKPOINTING_INTERVAL = 1;
    public static final String RUN_CHECKPOINTING_INTERVAL_SHORT_NAME = "RCPI";
    public static final String RUN_CHECKPOINTING_INTERVAL_LONG_NAME = "runCheckpointingInterval";

    public static final boolean DEFAULT_RUN_CHECKPOINTING_ENABLED = false;
    public static final String RUN_CHECKPOINTING_ENABLED_SHORT_NAME = "RCP";
    public static final String RUN_CHECKPOINTING_ENABLED_LONG_NAME = "runCheckpointing";

    public static final String DEFAULT_RUN_CHECKPOINTING_PATH = "/dev/null";
    public static final String RUN_CHECKPOINTING_PATH_SHORT_NAME = "RCPP";
    public static final String RUN_CHECKPOINTING_PATH_LONG_NAME = "runCheckpointingPath";

    public static final boolean DEFAULT_EXTENDED_POSTERIOR_OUTPUT_ENABLED = true;
    public static final String EXTENDED_POSTERIOR_OUTPUT_ENABLED_SHORT_NAME = "XPO";
    public static final String EXTENDED_POSTERIOR_OUTPUT_ENABLED_LONG_NAME = "extendedPosteriorOutputEnabled";

    @Advanced
    @Argument(
            doc = "Psi upper limit",
            shortName = PSI_UPPER_LIMIT_SHORT_NAME,
            fullName = PSI_UPPER_LIMIT_LONG_NAME,
            optional = true
    )
    protected double psiUpperLimit = DEFAULT_PSI_UPPER_LIMIT;

    @Advanced
    @Argument(
            doc = "Gamma upper limit",
            shortName = GAMMA_UPPER_LIMIT_SHORT_NAME,
            fullName = GAMMA_UPPER_LIMIT_LONG_NAME,
            optional = true
    )
    protected double gammaUpperLimit = DEFAULT_GAMMA_UPPER_LIMIT;

    @Advanced
    @Argument(
            doc = "Gamma solver number of bisections",
            shortName = GAMMA_SOLVER_NUM_BISECTIONS_SHORT_NAME,
            fullName = GAMMA_SOLVER_NUM_BISECTIONS_LONG_NAME,
            optional = true
    )
    protected int gammaSolverNumBisections = DEFAULT_GAMMA_SOLVER_NUM_BISECTIONS;

    @Advanced
    @Argument(
            doc = "Psi solver number of bisections",
            shortName = PSI_SOLVER_NUM_BISECTIONS_SHORT_NAME,
            fullName = PSI_SOLVER_NUM_BISECTIONS_LONG_NAME,
            optional = true
    )
    protected int psiSolverNumBisections = DEFAULT_PSI_SOLVER_NUM_BISECTIONS;

    @Advanced
    @Argument(
            doc = "Gamma solver grid refinement depth",
            shortName = GAMMA_SOLVER_REFINEMENT_DEPTH_SHORT_NAME,
            fullName = GAMMA_SOLVER_REFINEMENT_DEPTH_LONG_NAME,
            optional = true
    )
    protected int gammaSolverRefinementDepth = DEFAULT_GAMMA_SOLVER_REFINEMENT_DEPTH;

    @Advanced
    @Argument(
            doc = "Psi solver grid refinement depth",
            shortName = PSI_SOLVER_REFINEMENT_DEPTH_SHORT_NAME,
            fullName = PSI_SOLVER_REFINEMENT_DEPTH_LONG_NAME,
            optional = true
    )
    protected int psiSolverRefinementDepth = DEFAULT_PSI_SOLVER_REFINEMENT_DEPTH;


    @Argument(
            doc = "Maximum EM iterations",
            shortName = MAX_EM_ITERATIONS_SHORT_NAME,
            fullName = MAX_EM_ITERATIONS_LONG_NAME,
            optional = true
    )
    protected int maxEMIterations = DEFAULT_MAX_EM_ITERATIONS;

    @Argument(
            doc = "Dimension of the latent space",
            shortName = NUM_LATENTS_SHORT_NAME,
            fullName = NUM_LATENTS_LONG_NAME,
            optional = true
    )
    protected int numLatents = DEFAULT_NUM_LATENTS;

    @Argument(
            doc = "Log likelihood absolute error tolerance stopping criterion",
            shortName = LOG_LIKELIHOOD_TOL_SHORT_NAME,
            fullName = LOG_LIKELIHOOD_TOL_LONG_NAME,
            optional = true
    )
    protected double logLikelihoodTol = DEFAULT_LOG_LIKELIHOOD_TOL;

    @Argument(
            doc = "Model parameters absolute error tolerance stopping criterion",
            shortName = PARAM_ABS_TOL_SHORT_NAME,
            fullName = PARAM_ABS_TOL_LONG_NAME,
            optional = true
    )
    protected double paramAbsTol = DEFAULT_PARAM_ABS_TOL;

    @Argument(
            doc = "Posterior absolute error tolerance E-step cycle termination criterion",
            shortName = POSTERIOR_ABS_TOL_SHORT_NAME,
            fullName = POSTERIOR_ABS_TOL_LONG_NAME,
            optional = true
    )
    protected double posteriorAbsTol = DEFAULT_POSTERIOR_ABS_TOL;

    @Argument(
            doc = "Number of sequential maximization in the M-step",
            shortName = MAX_M_STEP_CYCLES_SHORT_NAME,
            fullName = MAX_M_STEP_CYCLES_LONG_NAME,
            optional = true
    )
    protected int maxMStepCycles = DEFAULT_MAX_M_STEP_CYCLES;

    @Argument(
            doc = "Maximum number of E-step cycles",
            shortName = MAX_E_STEP_CYCLES_SHORT_NAME,
            fullName = MAX_E_STEP_CYCLES_LONG_NAME,
            optional = true
    )
    protected int maxEStepCycles = DEFAULT_MAX_E_STEP_CYCLES;

    @Advanced
    @Argument(
            doc = "Log likelihood absolute change tolerance before updating copy ratio posteriors",
            shortName = LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING_SHORT_NAME,
            fullName = LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING_LONG_NAME,
            optional = true
    )
    protected double logLikelihoodTolThresholdCRCalling = DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING;

    @Advanced
    @Argument(
            doc = "Log likelihood absolute change tolerance before switching to target-resolved unexplained variance updates",
            shortName = LOG_LIKELIHOOD_TOL_THRESHOLD_PSI_SWITCHING_SHORT_NAME,
            fullName = LOG_LIKELIHOOD_TOL_THRESHOLD_PSI_SWITCHING_LONG_NAME,
            optional = true
    )
    protected double logLikelihoodTolThresholdPsiSwitching = DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_PSI_SWITCHING;

    @Advanced
    @Argument(
            doc = "Log likelihood absolute change tolerance before updating ARD coefficients",
            shortName = LOG_LIKELIHOOD_TOL_THRESHOLD_ARD_UPDATE_SHORT_NAME,
            fullName = LOG_LIKELIHOOD_TOL_THRESHOLD_ARD_UPDATE_LONG_NAME,
            optional = true
    )
    protected double logLikelihoodTolThresholdARDUpdate = DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_ARD_UPDATE;

    @Advanced
    @Argument(
            doc = "Minimum read count on a target to use it for learning",
            shortName = MIN_LEARNING_READ_COUNT_SHORT_NAME,
            fullName = MIN_LEARNING_READ_COUNT_LONG_NAME,
            optional = true
    )
    protected int minLearningReadCount = DEFAULT_MIN_LEARNING_READ_COUNT;

    @Advanced
    @Argument(
            doc = "Minimum read count on a target to use it for PCA initialization (if chosen)",
            shortName = MIN_PCA_INIT_READ_COUNT_SHORT_NAME,
            fullName = MIN_PCA_INIT_READ_COUNT_LONG_NAME,
            optional = true
    )
    protected int minPCAInitializationReadCount = DEFAULT_MIN_PCA_INIT_READ_COUNT;

    @Argument(
            doc = "Typical mapping error rate (used for safeguarding CNV calls on low-coverage targets)",
            shortName = MAPPING_ERROR_RATE_SHORT_NAME,
            fullName = MAPPING_ERROR_RATE_LONG_NAME,
            optional = true
    )
    protected double mappingErrorRate = DEFAULT_MAPPING_ERROR_RATE;

    @Argument(
            doc = "Enable CNV-avoiding regularization",
            shortName = FOURIER_REGULARIZATION_ENABLED_SHORT_NAME,
            fullName = FOURIER_REGULARIZATION_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean useFourierRegularization = DEFAULT_FOURIER_REGULARIZATION_ENABLED;

    @Argument(
            doc = "Minimum length of CNV event (for regularization)",
            shortName = MIN_CNV_LENGTH_SHORT_NAME,
            fullName = MIN_CNV_LENGTH_LONG_NAME,
            optional = true
    )
    protected int minCNVLength = DEFAULT_MIN_CNV_LENGTH;

    @Argument(
            doc = "Maximum length of CNV event (for regularization)",
            shortName = MAX_CNV_LENGTH_SHORT_NAME,
            fullName = MAX_CNV_LENGTH_LONG_NAME,
            optional = true
    )
    protected int maxCNVLength = DEFAULT_MAX_CNV_LENGTH;

    @Advanced
    @Argument(
            doc = "Zero pad data before FFT and promote the size to powers of 2",
            shortName = ZERO_PAD_FFT_SHORT_NAME,
            fullName = ZERO_PAD_FFT_LONG_NAME,
            optional = true
    )
    protected boolean zeroPadFFT = DEFAULT_ZERO_PAD_FFT;

    @Advanced
    @Argument(
            doc = "CNV-avoidance regularization strength",
            shortName = FOURIER_REGULARIZATION_STRENGTH_SHORT_NAME,
            fullName = FOURIER_REGULARIZATION_STRENGTH_LONG_NAME,
            optional = true
    )
    protected double fourierRegularizationStrength = DEFAULT_FOURIER_REGULARIZATION_STRENGTH;

    @Advanced
    @Argument(
            doc = "Psi solver mode",
            shortName = PSI_SOLVER_MODE_SHORT_NAME,
            fullName = PSI_SOLVER_MODE_LONG_NAME,
            optional = true
    )
    protected PsiUpdateMode psiUpdateMode = DEFAULT_PSI_SOLVER_MODE;

    @Argument(
            doc = "W solver type (local vs. spark)",
            shortName = W_SOLVER_TYPE_SHORT_NAME,
            fullName = W_SOLVER_TYPE_LONG_NAME,
            optional = true
    )
    protected WSolverStrategy wSolverStrategy = DEFAULT_W_SOLVER_TYPE;

    @Argument(
            doc = "Copy ratio HMM type (local vs. spark)",
            shortName = CR_HMM_TYPE_SHORT_NAME,
            fullName = CR_HMM_TYPE_LONG_NAME,
            optional = true
    )
    protected CopyRatioHMMType copyRatioHMMType = DEFAULT_CR_HMM_TYPE;

    @Argument(
            doc = "Maximum Psi solver iterations",
            shortName = MAX_PSI_ITERATIONS_SHORT_NAME,
            fullName = MAX_PSI_ITERATIONS_LONG_NAME,
            optional = true
    )
    protected int psiMaxIterations = DEFAULT_MAX_PSI_ITERATIONS;

    @Argument(
            doc = "Maximum W solver iterations (if regularization is enabled)",
            shortName = MAX_W_ITERATIONS_SHORT_NAME,
            fullName = MAX_W_ITERATIONS_LONG_NAME,
            optional = true
    )
    protected int wMaxIterations = DEFAULT_MAX_W_ITERATIONS;

    @Argument(
            doc = "Maximum gamma solver iterations",
            shortName = MAX_GAMMA_ITERATIONS_SHORT_NAME,
            fullName = MAX_GAMMA_ITERATIONS_LONG_NAME,
            optional = true
    )
    protected int gammaMaxIterations = DEFAULT_MAX_GAMMA_ITERATIONS;

    @Argument(
            doc = "W solver absolute error tolerance termination criterion (if regularization is enabled)",
            shortName = W_ABS_TOL_SHORT_NAME,
            fullName = W_ABS_TOL_LONG_NAME,
            optional = true
    )
    protected double wAbsTol = DEFAULT_W_ABS_TOL;

    @Argument(
            doc = "W solver relative error tolerance termination criterion (if regularization is enabled)",
            shortName = W_REL_TOL_SHORT_NAME,
            fullName = W_REL_TOL_LONG_NAME,
            optional = true
    )
    protected double wRelTol = DEFAULT_W_REL_TOL;

    @Argument(
            doc = "Psi solver absolute error tolerance termination criterion",
            shortName = PSI_ABS_TOL_SHORT_NAME,
            fullName = PSI_ABS_TOL_LONG_NAME,
            optional = true
    )
    protected double psiAbsTol = DEFAULT_PSI_ABS_TOL;

    @Argument(
            doc = "Psi solver relative error tolerance termination criterion",
            shortName = PSI_REL_TOL_SHORT_NAME,
            fullName = PSI_REL_TOL_LONG_NAME,
            optional = true
    )
    protected double psiRelTol = DEFAULT_PSI_REL_TOL;

    @Argument(
            doc = "Gamma solver absolute error tolerance termination criterion",
            shortName = GAMMA_ABS_TOL_SHORT_NAME,
            fullName = GAMMA_ABS_TOL_LONG_NAME,
            optional = true
    )
    protected double gammaAbsTol = DEFAULT_GAMMA_ABS_TOL;

    @Argument(
            doc = "Gamma solver relative error tolerance termination criterion",
            shortName = GAMMA_REL_TOL_SHORT_NAME,
            fullName = GAMMA_REL_TOL_LONG_NAME,
            optional = true
    )
    protected double gammaRelTol = DEFAULT_GAMMA_REL_TOL;

    @Argument(
            doc = "RDD checkpointing interval (in spark mode)",
            shortName = RDD_CHECKPOINTING_INTERVAL_SHORT_NAME,
            fullName = RDD_CHECKPOINTING_INTERVAL_LONG_NAME,
            optional = true
    )
    protected int rddCheckpointingInterval = DEFAULT_RDD_CHECKPOINTING_INTERVAL;

    @Argument(
            doc = "Model checkpointing interval",
            shortName = RUN_CHECKPOINTING_INTERVAL_SHORT_NAME,
            fullName = RUN_CHECKPOINTING_INTERVAL_LONG_NAME,
            optional = true
    )
    protected int runCheckpointingInterval = DEFAULT_RUN_CHECKPOINTING_INTERVAL;

    @Advanced
    @Argument(
            doc = "W communication policy (in spark mode)",
            shortName = W_COMMUNICATION_POLICY_SHORT_NAME,
            fullName = W_COMMUNICATION_POLICY_LONG_NAME,
            optional = true
    )
    protected ComputeNodeCommunicationPolicy principalMapComputeNodeCommunicationPolicy = DEFAULT_W_COMMUNICATION_POLICY;


    @Advanced
    @Argument(
            doc = "Admixing ratio for E-step mean-field equations",
            shortName = MEAN_FIELD_ADMIXING_RATIO_SHORT_NAME,
            fullName = MEAN_FIELD_ADMIXING_RATIO_LONG_NAME,
            optional = true
    )
    protected double meanFieldAdmixingRatio = DEFAULT_MEAN_FIELD_ADMIXING_RATIO;

    @Argument(
            doc = "Enable RDD checkpointing (in spark mode; recommended for stability)",
            shortName = RDD_CHECKPOINTING_ENABLED_SHORT_NAME,
            fullName = RDD_CHECKPOINTING_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean rddCheckpointingEnabled = DEFAULT_RDD_CHECKPOINTING_ENABLED;

    @Argument(
            doc = "Enable model checkpointing",
            shortName = RUN_CHECKPOINTING_ENABLED_SHORT_NAME,
            fullName = RUN_CHECKPOINTING_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean runCheckpointingEnabled = DEFAULT_RUN_CHECKPOINTING_ENABLED;

    @Argument(
            doc = "Enable gamma updates",
            shortName = GAMMA_UPDATE_ENABLED_SHORT_NAME,
            fullName = GAMMA_UPDATE_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean gammaUpdateEnabled = DEFAULT_GAMMA_UPDATE_ENABLED;

    @Argument(
            doc = "Enable psi updates",
            shortName = PSI_UPDATE_ENABLED_SHORT_NAME,
            fullName = PSI_UPDATE_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean psiUpdateEnabled = DEFAULT_PSI_UPDATE_ENABLED;

    @Advanced
    @Argument(
            doc = "Enable adaptive switiching of Psi solver modes (isotropic to target-resolved)",
            shortName = ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_SHORT_NAME,
            fullName = ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean adaptivePsiSolverModeSwitchingEnabled = DEFAULT_ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED;

    @Argument(
            doc = "Enable copy ratio updates",
            shortName = CR_UPDATE_ENABLED_SHORT_NAME,
            fullName = CR_UPDATE_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean copyRatioUpdateEnabled = DEFAULT_CR_UPDATE_ENABLED;

    @Argument(
            doc = "RDD checkpointing path (required if RDD checkpointing is enabled)",
            shortName = RDD_CHECKPOINTING_PATH_SHORT_NAME,
            fullName = RDD_CHECKPOINTING_PATH_LONG_NAME,
            optional = true
    )
    protected String rddCheckpointingPath = DEFAULT_RDD_CHECKPOINTING_PATH;

    @Argument(
            doc = "Model checkpointing path (required if model checkpointing is enabled)",
            shortName = RUN_CHECKPOINTING_PATH_SHORT_NAME,
            fullName = RUN_CHECKPOINTING_PATH_LONG_NAME,
            optional = true
    )
    protected String runCheckpointingPath = DEFAULT_RUN_CHECKPOINTING_PATH;

    @Advanced
    @Argument(
            doc = "Enable extended posterior output",
            shortName = EXTENDED_POSTERIOR_OUTPUT_ENABLED_SHORT_NAME,
            fullName = EXTENDED_POSTERIOR_OUTPUT_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean extendedPosteriorOutputEnabled = DEFAULT_EXTENDED_POSTERIOR_OUTPUT_ENABLED;

    @Advanced
    @Argument(
            doc = "Number of target space partitions (for spark mode)",
            shortName = NUMBER_OF_TARGET_SPACE_PARTITIONS_SHORT_NAME,
            fullName = NUMBER_OF_TARGET_SPACE_PARTITIONS_LONG_NAME,
            optional = true
    )
    protected int numTargetSpacePartitions = DEFAULT_NUMBER_OF_TARGET_SPACE_PARTITIONS;

    @Argument(
            doc = "Enable automatic relevance determination (ARD) of bias covariates",
            shortName = ARD_ENABLED_SHORT_NAME,
            fullName = ARD_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean ardEnabled = DEFAULT_ARD_ENABLED;

    @Advanced
    @Argument(
            doc = "Initial absolute ARD precision (if RANDOM model initialization is selected)",
            shortName = INITIAL_ARD_PRECISION_ABSOLUTE_SHORT_NAME,
            fullName = INITIAL_ARD_PRECISION_ABSOLUTE_LONG_NAME,
            optional = true
    )
    protected double initialARDPrecisionAbsolute = DEFAULT_INITIAL_ARD_PRECISION_ABSOLUTE;

    @Advanced
    @Argument(
            doc = "Initial absolute ARD precision relative to noise (if PCA model initialization is selected)",
            shortName = INITIAL_ARD_PRECISION_RELATIVE_TO_NOISE_SHORT_NAME,
            fullName = INITIAL_ARD_PRECISION_RELATIVE_TO_NOISE_LONG_NAME,
            optional = true
    )
    protected double initialARDPrecisionRelativeToNoise = DEFAULT_INITIAL_ARD_PRECISION_RELATIVE_TO_NOISE;

    @Advanced
    @Argument(
            doc = "Maximum ARD precision",
            shortName = MAX_ARD_PRECISION_SHORT_NAME,
            fullName = MAX_ARD_PRECISION_LONG_NAME,
            optional = true
    )
    protected double maxARDPrecision = DEFAULT_MAX_ARD_PRECISION;

    @Advanced
    @Argument(
            doc = "Number of threads for concurrent target-resolved unexpalained variance solver",
            shortName = PSI_SOLVER_NUM_THREADS_SHORT_NAME,
            fullName = PSI_SOLVER_NUM_THREADS_LONG_NAME,
            optional = true
    )
    protected int psiSolverNumThreads = DEFAULT_PSI_SOLVER_NUM_THREADS;

    @Argument(
            doc = "Model initialization strategy (only used in learning tasks)",
            shortName = MODEL_INITIALIZATION_STRATEGY_SHORT_NAME,
            fullName = MODEL_INITIALIZATION_STRATEGY_LONG_NAME,
            optional = true
    )
    protected ModelInitializationStrategy modelInitializationStrategy = DEFAULT_MODEL_INITIALIZATION_STRATEGY;

    @Argument(
            doc = "Include the contribution of ARD in the reported log likelihood values",
            shortName = INCLUDE_ARD_IN_LOG_LIKELIHOOD_SHORT_NAME,
            fullName = INCLUDE_ARD_IN_LOG_LIKELIHOOD_LONG_NAME,
            optional = true
    )
    protected boolean includeARDInLogLikelihood = DEFAULT_INCLUDE_ARD_IN_LOG_LIKELIHOOD;

    /*
     * setters and getters
     */

    public int getMaxEMIterations() { return maxEMIterations; }

    public int getNumLatents() { return numLatents; }

    public double getLogLikelihoodTolerance() { return logLikelihoodTol; }

    public int getMaxMStepCycles() { return maxMStepCycles; }

    public boolean fourierRegularizationEnabled() { return useFourierRegularization; }

    public double getFourierRegularizationStrength() { return fourierRegularizationStrength; }

    public double getPsiAbsoluteTolerance() { return psiAbsTol; }

    public double getPsiRelativeTolerance() { return psiRelTol; }

    public int getPsiMaxIterations() { return psiMaxIterations; }

    public double getWAbsoluteTolerance() { return this.wAbsTol; }

    public double getWRelativeTolerance() { return this.wRelTol; }

    public int getWMaxIterations() { return wMaxIterations; }

    public boolean zeroPadFFT() {
        return zeroPadFFT;
    }

    public double getParameterEstimationAbsoluteTolerance() { return this.paramAbsTol; }

    public PsiUpdateMode getPsiUpdateMode() {
        return psiUpdateMode;
    }

    public CoverageModelEMParams setPsiUpdateMode(@Nonnull final PsiUpdateMode psiUpdateMode) {
        this.psiUpdateMode = Utils.nonNull(psiUpdateMode, "Psi solver mode must be non-null");
        return this;
    }

    public WSolverStrategy getWSolverType() {
        return wSolverStrategy;
    }

    public int getMinimumCNVLength() { return minCNVLength; }

    public int getMaximumCNVLength() { return maxCNVLength; }

    public boolean isRDDCheckpointingEnabled() {
        return rddCheckpointingEnabled;
    }

    public int getRDDCheckpointingInterval() {
        return rddCheckpointingInterval;
    }

    public CopyRatioHMMType getCopyRatioHMMType() {
        return copyRatioHMMType;
    }

    public double getLogLikelihoodTolThresholdCRCalling() {
        return logLikelihoodTolThresholdCRCalling;
    }

    public double getLogLikelihoodTolThresholdPsiSwitching() {
        return logLikelihoodTolThresholdPsiSwitching;
    }

    public double getLogLikelihoodTolThresholdARDUpdate() {
        return logLikelihoodTolThresholdARDUpdate;
    }

    public double getPosteriorAbsTol() {
        return posteriorAbsTol;
    }

    public int getMaxEStepCycles() {
        return maxEStepCycles;
    }

    public ComputeNodeCommunicationPolicy getBiasCovariatesComputeNodeCommunicationPolicy() {
        return principalMapComputeNodeCommunicationPolicy;
    }

    public double getMeanFieldAdmixingRatio() {
        return meanFieldAdmixingRatio;
    }

    public int getRunCheckpointingInterval() {
        return runCheckpointingInterval;
    }

    public double getGammaAbsoluteTolerance() {
        return gammaAbsTol;
    }

    public double getGammaRelativeTolerance() {
        return gammaRelTol;
    }

    public int getGammaMaximumIterations() {
        return gammaMaxIterations;
    }

    public boolean gammaUpdateEnabled() {
        return gammaUpdateEnabled;
    }

    public boolean psiUpdateEnabled() {
        return psiUpdateEnabled;
    }

    public boolean copyRatioUpdateEnabled() {
        return copyRatioUpdateEnabled;
    }

    public boolean adaptivePsiSolverModeSwitchingEnabled() {
        return adaptivePsiSolverModeSwitchingEnabled;
    }

    public double getPsiUpperLimit() {
        return psiUpperLimit;
    }

    public double getGammaUpperLimit() {
        return gammaUpperLimit;
    }

    public boolean isRunCheckpointingEnabled() {
        return runCheckpointingEnabled;
    }

    public String getRunCheckpointingPath() {
        return runCheckpointingPath;
    }

    public String getRDDCheckpointingPath() {
        return rddCheckpointingPath;
    }

    public boolean extendedPosteriorOutputEnabled() {
        return extendedPosteriorOutputEnabled;
    }

    public int getNumTargetSpacePartitions() {
        return numTargetSpacePartitions;
    }

    public int getGammaSolverRefinementDepth() {
        return gammaSolverRefinementDepth;
    }

    public int getGammaSolverNumBisections() {
        return gammaSolverNumBisections;
    }

    public int getPsiSolverRefinementDepth() {
        return psiSolverRefinementDepth;
    }

    public int getPsiSolverNumBisections() {
        return psiSolverNumBisections;
    }

    public int getMinLearningReadCount() {
        return minLearningReadCount;
    }

    public int getMinPCAInitializationReadCount() {
        return minPCAInitializationReadCount;
    }

    public double getMappingErrorRate() {
        return mappingErrorRate;
    }

    public boolean isARDEnabled() {
        return ardEnabled;
    }

    public double getInitialARDPrecisionAbsolute() {
        return initialARDPrecisionAbsolute;
    }

    public double getInitialARDPrecisionRelativeToNoise() {
        return initialARDPrecisionRelativeToNoise;
    }

    public double getMaxARDPrecision() {
        return maxARDPrecision;
    }

    public int getPsiSolverNumThreads() { return psiSolverNumThreads; }

    public ModelInitializationStrategy getModelInitializationStrategy() {
        return modelInitializationStrategy;
    }

    public boolean includeARDInLogLikelihood() {
        return includeARDInLogLikelihood;
    }

    /**
     * Validate parameters.
     *
     * TODO github/gatk-protected issue #843 -- more validations
     *      positive/negative values
     *      Maxes greater than mins
     *
     */
    public void validate() {

        ParamUtils.isPositive(maxEMIterations, "Maximum EM iterations must be positive");
        ParamUtils.isPositiveOrZero(numLatents, "Number of latent variables must be non-negative");
        ParamUtils.isPositive(logLikelihoodTol, "Convergence tolerance on log likelihood must be positive");
        ParamUtils.isPositive(maxMStepCycles, "The number of sequential partial maximization steps must be positive");
        ParamUtils.isPositive(fourierRegularizationStrength, "The Fourier regularization strength must be positive");
        ParamUtils.isPositive(psiAbsTol, "The absolute tolerance for maximization of Psi must be positive");
        ParamUtils.isPositive(psiRelTol, "The relative tolerance for maximization of Psi must be positive");
        ParamUtils.isPositive(psiMaxIterations, "The maximum number of iterations for M-step of Psi must be positive");
        ParamUtils.isPositive(wAbsTol, "The absolute tolerance for maximization of Psi must be positive");
        ParamUtils.isPositive(wRelTol, "The relative tolerance for maximization of Psi must be positive");
        ParamUtils.isPositive(wMaxIterations, "The maximum number of iterations for M-step of W must be positive.");
        ParamUtils.isPositive(paramAbsTol, "The required tolerance on parameter change must be positive.");
        ParamUtils.isPositive(minCNVLength, "Minimum CNV length must be positive");
        ParamUtils.inRange(maxCNVLength, minCNVLength, Integer.MAX_VALUE, "Maximum CNV length must be greater than the minimum");
        ParamUtils.inRange(rddCheckpointingInterval, 1, Integer.MAX_VALUE, "RDD checkpointing interval must be >= 1");
        Utils.nonNull(copyRatioHMMType, "Copy ratio HMM type must be non-null");
        ParamUtils.isPositive(logLikelihoodTolThresholdCRCalling, "Log likelihood change threshold before updating" +
                " copy ratio posteriors must be positive");
        ParamUtils.isPositive(logLikelihoodTolThresholdPsiSwitching, "Log likelihood change threshold before switching" +
                " to target-resolved unexplained variance updates must be positive");
        ParamUtils.isPositive(logLikelihoodTolThresholdARDUpdate, "Log likelihood change threshold before updating" +
                " ARD coefficients must be positive");
        ParamUtils.isPositive(posteriorAbsTol, "Posterior absolute error tolerance must be positive");
        ParamUtils.inRange(maxEStepCycles, 1, Integer.MAX_VALUE, "Maximum number of E-step cycles  must be positive");
        Utils.nonNull(principalMapComputeNodeCommunicationPolicy);
        ParamUtils.inRange(meanFieldAdmixingRatio, 0, 1, "The mean-field admixing ratio must be between 0 and 1");
        ParamUtils.inRange(runCheckpointingInterval, 1, Integer.MAX_VALUE, "Model checkpointing intervals must be >= 1");
        ParamUtils.isPositive(gammaAbsTol, "Gamma absolute error tolerance must be positive");
        ParamUtils.isPositive(gammaRelTol, "Gamma relative error tolerance must be positive");
        ParamUtils.isPositive(gammaMaxIterations, "Gamma solver maximum iterations must be positive");
        ParamUtils.isPositive(psiUpperLimit, "Psi upper limit must be positive");
        ParamUtils.isPositive(gammaUpperLimit, "Gamma upper limit must be positive");
        Utils.nonNull(runCheckpointingPath, "Run checkpointing path must be non-null");
        Utils.nonNull(rddCheckpointingPath, "RDD checkpointing path must be non-null");
        ParamUtils.isPositive(numTargetSpacePartitions, "Number of target space partitions must be positive");
        ParamUtils.isPositive(minLearningReadCount, "The minimum learning read count must be positive");
        ParamUtils.isPositive(minPCAInitializationReadCount, "The minimum PCA initialization read count must be positive");
        ParamUtils.isPositiveOrZero(mappingErrorRate, "The mapping error rate must be non-negative");
        ParamUtils.isPositive(initialARDPrecisionAbsolute, "The absolute initial ARD precision must be positive");
        ParamUtils.isPositive(initialARDPrecisionRelativeToNoise, "The relative initial ARD precision must be positive");
        ParamUtils.isPositive(maxARDPrecision, "The maximum ARD precision must be positive");
        ParamUtils.isPositive(psiSolverNumThreads, "Number of psi solver threads must be positive");

        Utils.validateArg(!isRunCheckpointingEnabled() || !runCheckpointingPath.equals("/dev/null"),
                "Run checkpointing is enabled but checkpointing path is not set properly");
        Utils.validateArg(!fourierRegularizationEnabled(), "Fourier regularization is not properly" +
                " implemented yet");
        Utils.validateArg(numLatents > 0 || !ardEnabled, "ARD must be disabled if the dimension of the" +
                " bias latent space is zero");

    }
}
