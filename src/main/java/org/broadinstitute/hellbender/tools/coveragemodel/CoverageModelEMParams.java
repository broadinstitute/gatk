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

    public static final double DEFAULT_LOG_LIKELIHOOD_TOL = 1e-5;
    public static final String LOG_LIKELIHOOD_TOL_SHORT_NAME = "LLT";
    public static final String LOG_LIKELIHOOD_TOL_LONG_NAME = "logLikelihoodTol";

    public static final double DEFAULT_PARAM_ABS_TOL = 1e-4;
    public static final String PARAM_ABS_TOL_SHORT_NAME = "PMAT";
    public static final String PARAM_ABS_TOL_LONG_NAME = "paramAbsoluteTolerance";

    public static final double DEFAULT_POSTERIOR_ABS_TOL = 1e-2;
    public static final String POSTERIOR_ABS_TOL_SHORT_NAME = "PSAT";
    public static final String POSTERIOR_ABS_TOL_LONG_NAME = "posteriorAbsoluteTolerance";

    public static final double DEFAULT_MEAN_FIELD_ADMIXING_RATIO = 0.75;
    public static final String MEAN_FIELD_ADMIXING_RATIO_SHORT_NAME = "MFAR";
    public static final String MEAN_FIELD_ADMIXING_RATIO_LONG_NAME = "meanFieldAdmixingRatio";

    public static final int DEFAULT_MAX_M_STEP_CYCLES = 1;
    public static final String MAX_M_STEP_CYCLES_SHORT_NAME = "MMC";
    public static final String MAX_M_STEP_CYCLES_LONG_NAME = "maximumMStepCycles";

    public static final int DEFAULT_MAX_E_STEP_CYCLES = 4;
    public static final String MAX_E_STEP_CYCLES_SHORT_NAME = "MES";
    public static final String MAX_E_STEP_CYCLES_LONG_NAME = "maximumEStepCycles";

    public static final int DEFAULT_MAX_EM_ITERATIONS = 20;
    public static final String MAX_EM_ITERATIONS_SHORT_NAME = "MEMI";
    public static final String MAX_EM_ITERATIONS_LONG_NAME = "maximumEMIterations";

    public static final double DEFAULT_LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING = 5e-2;
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING_SHORT_NAME = "LLTTCRC";
    public static final String LOG_LIKELIHOOD_TOL_THRESHOLD_CR_CALLING_LONG_NAME = "logLikelihoodTolThresholdCopyRatioCalling";

    /**
     * Number of bias continuous latent variables
     */
    public static final int DEFAULT_NUM_LATENTS = 5;
    public static final String NUM_LATENTS_SHORT_NAME = "NL";
    public static final String NUM_LATENTS_LONG_NAME = "numLatents";

    /**
     * Minimum read count on a target to perform learning
     */
    public static final int DEFAULT_MIN_LEARNING_READ_COUNT = 20;
    public static final String MIN_LEARNING_READ_COUNT_SHORT_NAME = "MLRC";
    public static final String MIN_LEARNING_READ_COUNT_LONG_NAME = "minimumLearningReadCount";

    public static final double DEFAULT_MAPPING_ERROR_RATE = 1e-3;
    public static final String MAPPING_ERROR_RATE_SHORT_NAME = "MER";
    public static final String MAPPING_ERROR_RATE_LONG_NAME = "mappingErrorRate";




    public static final boolean DEFAULT_PSI_UPDATE_ENABLED = true;
    public static final String PSI_UPDATE_ENABLED_SHORT_NAME = "PU";
    public static final String PSI_UPDATE_ENABLED_LONG_NAME = "psiUpdate";

    public static final PsiUpdateMode DEFAULT_PSI_SOLVER_MODE = PsiUpdateMode.PSI_TARGET_RESOLVED;
    public static final String PSI_SOLVER_MODE_SHORT_NAME = "PSM";
    public static final String PSI_SOLVER_MODE_LONG_NAME = "psiUpdateMode";

    public static final boolean DEFAULT_ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED = true;
    public static final String ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_SHORT_NAME = "APSMS";
    public static final String ADAPTIVE_PSI_SOLVER_MODE_SWITCHING_ENABLED_LONG_NAME = "adaptivePsiSolverModeSwitching";

    public static final double DEFAULT_PSI_UPPER_LIMIT = 1.0;
    public static final String PSI_UPPER_LIMIT_SHORT_NAME = "PUL";
    public static final String PSI_UPPER_LIMIT_LONG_NAME = "psiUpperLimit";

    public static final int DEFAULT_MAX_PSI_ITERATIONS = 200;
    public static final String MAX_PSI_ITERATIONS_SHORT_NAME = "MPI";
    public static final String MAX_PSI_ITERATIONS_LONG_NAME = "maxPsiInterations";

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




    public static final WSolverStrategy DEFAULT_W_SOLVER_TYPE = WSolverStrategy.W_SOLVER_SPARK;
    public static final String W_SOLVER_TYPE_SHORT_NAME = "WST";
    public static final String W_SOLVER_TYPE_LONG_NAME = "wSolverStrategy";

    public static final ComputeNodeCommunicationPolicy DEFAULT_W_COMMUNICATION_POLICY = ComputeNodeCommunicationPolicy.RDD_JOIN;
    public static final String W_COMMUNICATION_POLICY_SHORT_NAME = "WCP";
    public static final String W_COMMUNICATION_POLICY_LONG_NAME = "wCommunicationPolicy";

    public static final boolean DEFAULT_W_ORTHOGONALIZATION_ENABLED = true;
    public static final String W_ORTHOGONALIZATION_ENABLED_SHORT_NAME = "WOE";
    public static final String W_ORTHOGONALIZATION_ENABLED_LONG_NAME = "wOrthogonalizationEnabled";

    public static final int DEFAULT_MAX_W_ITERATIONS = 200;
    public static final String MAX_W_ITERATIONS_SHORT_NAME = "MWI";
    public static final String MAX_W_ITERATIONS_LONG_NAME = "maxWInterations";

    public static final double DEFAULT_W_ABS_TOL = 1e-8;
    public static final String W_ABS_TOL_SHORT_NAME = "WAT";
    public static final String W_ABS_TOL_LONG_NAME = "wAbsoluteTolerance";

    public static final double DEFAULT_W_REL_TOL = 1e-6;
    public static final String W_REL_TOL_SHORT_NAME = "WRT";
    public static final String W_REL_TOL_LONG_NAME = "wRelativeTolerance";




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
            doc = "Minimum read count on a target to use it for learning",
            shortName = MIN_LEARNING_READ_COUNT_SHORT_NAME,
            fullName = MIN_LEARNING_READ_COUNT_LONG_NAME,
            optional = true
    )
    protected int minLearningReadCount = DEFAULT_MIN_LEARNING_READ_COUNT;

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
            doc = "Enable orthogonalization of principal map",
            shortName = W_ORTHOGONALIZATION_ENABLED_SHORT_NAME,
            fullName = W_ORTHOGONALIZATION_ENABLED_LONG_NAME,
            optional = true
    )
    protected boolean orthogonalizeAndSortPrincipalMapEnabled = DEFAULT_W_ORTHOGONALIZATION_ENABLED;

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
            doc = "Number of target space paritions (for spark mode)",
            shortName = NUMBER_OF_TARGET_SPACE_PARTITIONS_SHORT_NAME,
            fullName = NUMBER_OF_TARGET_SPACE_PARTITIONS_LONG_NAME,
            optional = true
    )
    protected int numTargetSpaceParititions = DEFAULT_NUMBER_OF_TARGET_SPACE_PARTITIONS;

    /**********************************************************************************
     * setters and getters (used for on-the-fly update of parameters and for testing) *
     **********************************************************************************/

    public CoverageModelEMParams setMaxEMIterations(final int maxEMIterations) {
        this.maxEMIterations = ParamUtils.isPositive(maxEMIterations, "Maximum EM iterations must be positive.");
        return this;
    }

    public int getMaxEMIterations() { return maxEMIterations; }

    public CoverageModelEMParams setNumLatents(final int numLatents) {
        this.numLatents = ParamUtils.isPositive(numLatents, "Number of latent variables must be positive.");
        return this;
    }

    public int getNumLatents() { return numLatents; }

    public CoverageModelEMParams setLogLikelihoodTolerance(final double tol) {
        logLikelihoodTol = ParamUtils.isPositive(tol, "The required tolerance on log likelihood " +
                "must be positive.");
        return this;
    }

    public double getLogLikelihoodTolerance() { return logLikelihoodTol; }

    public CoverageModelEMParams setMaxMStepCycles(final int maxMStepCycles) {
        this.maxMStepCycles = ParamUtils.isPositive(maxMStepCycles, "The number of " +
                "sequential partial maximimization steps must be positive.");
        return this;
    }

    public int getMaxMStepCycles() { return maxMStepCycles; }

    public CoverageModelEMParams enableFourierRegularization() {
        useFourierRegularization = true;
        return this;
    }

    public CoverageModelEMParams disableFourierRegularization() {
        useFourierRegularization = false;
        return this;
    }

    public boolean fourierRegularizationEnabled() { return useFourierRegularization; }

    public CoverageModelEMParams setFourierRegularizationStrength(final double fourierRegularizationStrength) {
        this.fourierRegularizationStrength = ParamUtils.isPositive(fourierRegularizationStrength, "The Fourier " +
                "regularization strength must be positive");
        return this;
    }

    public double getFourierRegularizationStrength() { return fourierRegularizationStrength; }

    public CoverageModelEMParams setPsiAbsoluteTolerance(final double tol) {
        this.psiAbsTol = ParamUtils.isPositive(tol, "The absolute tolerance for maximization of Psi must be positive");
        return this;
    }

    public double getPsiAbsoluteTolerance() { return psiAbsTol; }

    public CoverageModelEMParams setPsiRelativeTolerance(final double tol) {
        this.psiRelTol = ParamUtils.isPositive(tol, "The relative tolerance for maximization of Psi must be positive");
        return this;
    }

    public double getPsiRelativeTolerance() { return psiRelTol; }

    public CoverageModelEMParams setPsiMaxIterations(final int psiMaxIterations) {
        this.psiMaxIterations = ParamUtils.isPositive(psiMaxIterations, "The maximum number of interations for M-step of Psi " +
                "must be positive.");
        return this;
    }

    public int getPsiMaxIterations() { return psiMaxIterations; }

    public CoverageModelEMParams setWAbsoluteTolerance(final double tol) {
        this.wAbsTol = ParamUtils.isPositive(tol, "The absolute tolerance for maximization of Psi must be positive");
        return this;
    }

    public double getWAbsoluteTolerance() { return this.wAbsTol; }

    public CoverageModelEMParams setWRelativeTolerance(final double tol) {
        this.wRelTol = ParamUtils.isPositive(tol, "The relative tolerance for maximization of Psi must be positive");
        return this;
    }

    public double getWRelativeTolerance() { return this.wRelTol; }

    public CoverageModelEMParams setWMaxIterations(final int wMaxIterations) {
        this.wMaxIterations = ParamUtils.isPositive(wMaxIterations, "The maximum number of interations for M-step of W " +
                "must be positive.");
        return this;
    }

    public int getWMaxIterations() { return wMaxIterations; }

    public CoverageModelEMParams enableZeroPadFFT() {
        zeroPadFFT = true;
        return this;
    }

    public CoverageModelEMParams disableZeroPadFFT() {
        zeroPadFFT = false;
        return this;
    }

    public boolean zeroPadFFT() {
        return zeroPadFFT;
    }

    public CoverageModelEMParams setParameterEstimationAbsoluteTolerance(final double val) {
        this.paramAbsTol = ParamUtils.isPositive(paramAbsTol, "The required tolerance on parameter change must be positive.");
        return this;
    }

    public double getParameterEstimationAbsoluteTolerance() { return this.paramAbsTol; }

    public PsiUpdateMode getPsiUpdateMode() {
        return psiUpdateMode;
    }

    public CoverageModelEMParams setPsiPsiolverType(@Nonnull final PsiUpdateMode psiUpdateMode) {
        this.psiUpdateMode = Utils.nonNull(psiUpdateMode, "Psi solver mode must be non-null");
        return this;
    }

    public WSolverStrategy getWSolverType() {
        return wSolverStrategy;
    }

    public CoverageModelEMParams setWSolverType(@Nonnull final WSolverStrategy wSolverStrategy) {
        this.wSolverStrategy = Utils.nonNull(wSolverStrategy, "W solver type must be non-null");
        return this;
    }

    public CoverageModelEMParams setMinimumCNVLength(final int minCNVLength) {
        this.minCNVLength = ParamUtils.isPositive(minCNVLength, "Minimum CNV length must be positive");
        return this;
    }

    public CoverageModelEMParams setMaximumCNVLength(final int maxCNVLength) {
        this.maxCNVLength = ParamUtils.inRange(maxCNVLength, minCNVLength, Integer.MAX_VALUE, "Maximum CNV length must be greater than" +
                " the minimum");
        return this;
    }

    public int getMinimumCNVLength() { return minCNVLength; }

    public int getMaximumCNVLength() { return maxCNVLength; }

    public CoverageModelEMParams enableRDDCheckpointing() {
        this.rddCheckpointingEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableRDDCheckpointing() {
        this.rddCheckpointingEnabled = false;
        return this;
    }

    public boolean isRDDCheckpointingEnabled() {
        return rddCheckpointingEnabled;
    }

    public int getRDDCheckpointingInterval() {
        return rddCheckpointingInterval;
    }

    public CoverageModelEMParams setRddCheckpointingInterval(final int rddCheckpointingInterval) {
        this.rddCheckpointingInterval = ParamUtils.inRange(rddCheckpointingInterval, 1, Integer.MAX_VALUE,
            "RDD checkpointing interval must be >= 1");
        return this;
    }

    public CopyRatioHMMType getCopyRatioHMMType() {
        return copyRatioHMMType;
    }

    public CoverageModelEMParams setCopyRatioHMMType(final CopyRatioHMMType copyRatioHMMType) {
        this.copyRatioHMMType = Utils.nonNull(copyRatioHMMType, "Copy ratio HMM type must be non-null");
        return this;
    }

    public double getLogLikelihoodTolThresholdCRCalling() {
        return logLikelihoodTolThresholdCRCalling;
    }

    public CoverageModelEMParams setLogLikelihoodTolThresholdCRCalling(final double logLikelihoodTolThresholdCRCalling) {
        this.logLikelihoodTolThresholdCRCalling = ParamUtils.isPositive(logLikelihoodTolThresholdCRCalling,
                "Log likelihood change threshold before updating copy ratio posteriors must be positive");
        return this;
    }

    public double getPosteriorAbsTol() {
        return posteriorAbsTol;
    }

    public CoverageModelEMParams setPosteriorAbsTol(final double posteriorAbsTol) {
        this.posteriorAbsTol = ParamUtils.isPositive(posteriorAbsTol, "Posterior absolute error tolerance must be" +
                " positive");
        return this;
    }

    public int getMaxEStepCycles() {
        return maxEStepCycles;
    }

    public CoverageModelEMParams setMaxEStepCycles(final int maxEStepCycles) {
        this.maxEStepCycles = ParamUtils.inRange(maxEStepCycles, 1, Integer.MAX_VALUE, "Maximum number of E-step cycles" +
                " must be positive");
        return this;
    }

    public ComputeNodeCommunicationPolicy getBiasCovariatesComputeNodeCommunicationPolicy() {
        return principalMapComputeNodeCommunicationPolicy;
    }

    public CoverageModelEMParams setPrincipalMapComputeNodeCommunicationPolicy(final ComputeNodeCommunicationPolicy principalMapComputeNodeCommunicationPolicy) {
        this.principalMapComputeNodeCommunicationPolicy = Utils.nonNull(principalMapComputeNodeCommunicationPolicy);
        return this;
    }

    public double getMeanFieldAdmixingRatio() {
        return meanFieldAdmixingRatio;
    }

    public CoverageModelEMParams setMeanFieldAdmixingRatio(double meanFieldAdmixingRatio) {
        this.meanFieldAdmixingRatio = ParamUtils.inRange(meanFieldAdmixingRatio, 0, 1,
            "The mean-field admixing ratio must be between 0 and 1");
        return this;
    }

    public boolean isOrthogonalizeAndSortBiasCovariatesEnabled() {
        return orthogonalizeAndSortPrincipalMapEnabled;
    }

    public CoverageModelEMParams enableOrthogonalizeAndSortPrincipalMap() {
        orthogonalizeAndSortPrincipalMapEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableOrthogonalizeAndSortPrincipalMap() {
        orthogonalizeAndSortPrincipalMapEnabled = false;
        return this;
    }

    public int getRunCheckpointingInterval() {
        return runCheckpointingInterval;
    }

    public CoverageModelEMParams setRunCheckpointingInterval(final int runCheckpointingInterval) {
        this.runCheckpointingInterval = ParamUtils.inRange(runCheckpointingInterval, 1, Integer.MAX_VALUE,
            "Model checkpointing intervals must be >= 1");
        return this;
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

    public CoverageModelEMParams setGammaAbsoluteTolerance(final double gammaAbsTol) {
        this.gammaAbsTol = ParamUtils.isPositive(gammaAbsTol, "Gamma absolute error tolerance must be positive");
        return this;
    }

    public CoverageModelEMParams setGammaRelativeTolerance(final double gammaRelTol) {
        this.gammaRelTol = ParamUtils.isPositive(gammaRelTol, "Gamma relative error tolerance must be positive");
        return this;
    }

    public CoverageModelEMParams setGammaMaximumIterations(final int gammaMaxIterations) {
        this.gammaMaxIterations = ParamUtils.isPositive(gammaMaxIterations, "Gamma solver maximum iterations must be positive");
        return this;
    }

    public boolean gammaUpdateEnabled() {
        return gammaUpdateEnabled;
    }

    public CoverageModelEMParams enableGammaUpdate() {
        gammaUpdateEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableGammaUpdate() {
        gammaUpdateEnabled = false;
        return this;
    }


    public boolean psiUpdateEnabled() {
        return psiUpdateEnabled;
    }

    public CoverageModelEMParams enablePsiUpdate() {
        psiUpdateEnabled = true;
        return this;
    }

    public CoverageModelEMParams disablePsiUpdate() {
        psiUpdateEnabled = false;
        return this;
    }



    public boolean copyRatioUpdateEnabled() {
        return copyRatioUpdateEnabled;
    }

    public CoverageModelEMParams enableCopyRatioUpdate() {
        copyRatioUpdateEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableCopyRatioUpdate() {
        copyRatioUpdateEnabled = false;
        return this;
    }

    public boolean adaptivePsiSolverModeSwitchingEnabled() {
        return adaptivePsiSolverModeSwitchingEnabled;
    }

    public CoverageModelEMParams enableAdaptivePsiSolverModeSwitching() {
        adaptivePsiSolverModeSwitchingEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableAdaptivePsiSolverModeSwitching() {
        adaptivePsiSolverModeSwitchingEnabled = false;
        return this;
    }

    public CoverageModelEMParams setPsiUpperLimit(final double psiUpperLimit) {
        this.psiUpperLimit = ParamUtils.isPositive(psiUpperLimit, "Psi upper limit must be positive");
        return this;
    }

    public double getPsiUpperLimit() {
        return psiUpperLimit;
    }

    public CoverageModelEMParams setGammaUpperLimit(final double gammaUpperLimit) {
        this.gammaUpperLimit = ParamUtils.isPositive(gammaUpperLimit, "Gamma upper limit must be positive");
        return this;
    }

    public double getGammaUpperLimit() {
        return gammaUpperLimit;
    }







    public CoverageModelEMParams enableRunCheckpointing() {
        runCheckpointingEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableRunCheckpointing() {
        runCheckpointingEnabled = false;
        return this;
    }

    public boolean isRunCheckpointingEnabled() {
        return runCheckpointingEnabled;
    }

    public String getRunCheckpointingPath() {
        return runCheckpointingPath;
    }

    public CoverageModelEMParams setRunCheckpointingPath(final String runCheckpointingPath) {
        this.runCheckpointingPath = Utils.nonNull(runCheckpointingPath, "Run checkpointing path must be non-null");
        return this;
    }

    public String getRDDCheckpointingPath() {
        return rddCheckpointingPath;
    }

    public CoverageModelEMParams setRDDCheckpointingPath(final String rddCheckpointingPath) {
        this.rddCheckpointingPath = Utils.nonNull(rddCheckpointingPath, "RDD checkpointing path must be non-null");
        return this;
    }

    public CoverageModelEMParams enableExtendedPosteriorOutput() {
        this.extendedPosteriorOutputEnabled = true;
        return this;
    }

    public CoverageModelEMParams disableExtendedPosteriorOutput() {
        this.extendedPosteriorOutputEnabled = false;
        return this;
    }

    public boolean extendedPosteriorOutputEnabled() {
        return extendedPosteriorOutputEnabled;
    }

    public int getNumTargetSpaceParititions() {
        return numTargetSpaceParititions;
    }

    public CoverageModelEMParams setNumTargetSpacePartitions(final int numTargetSpaceParititions) {
        this.numTargetSpaceParititions = ParamUtils.isPositive(numTargetSpaceParititions, "Number of target space" +
                " partitions must be positive");
        return this;
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

    public CoverageModelEMParams setGammaSolverNumBisections(final int gammaSolverNumBisections) {
        this.gammaSolverNumBisections = gammaSolverNumBisections;
        return this;
    }

    public CoverageModelEMParams setGammaSolverRefinementDepth(final int gammaSolverRefinementDepth) {
        this.gammaSolverRefinementDepth = gammaSolverRefinementDepth;
        return this;
    }

    public CoverageModelEMParams setPsiSolverNumBisections(final int psiSolverNumBisections) {
        this.psiSolverNumBisections = psiSolverNumBisections;
        return this;
    }

    public CoverageModelEMParams setPsiSolverRefinementDepth(final int psiSolverRefinementDepth) {
        this.psiSolverRefinementDepth = psiSolverRefinementDepth;
        return this;
    }

    public CoverageModelEMParams setMinLearningReadCount(final int minLearningReacCount) {
        this.minLearningReadCount = ParamUtils.isPositive(minLearningReacCount,
                "The minimum learning read count must be positive");
        return this;
    }

    public int getMinLearningReadCount() {
        return minLearningReadCount;
    }

    public CoverageModelEMParams setMappingErrorRate(final double mappingErrorRate) {
        this.mappingErrorRate = ParamUtils.isPositiveOrZero(mappingErrorRate,
                "The mapping error rate must be non-negative");
        return this;
    }

    public double getMappingErrorRate() {
        return mappingErrorRate;
    }

    /**
     * Validate parameters
     *
     * TODO github/gatk-protected issue #843 -- more validations
     *
     */
    public void validate() {
        Utils.validateArg(!isRunCheckpointingEnabled() || !runCheckpointingPath.equals("/dev/null"),
                "Run checkpointing is enabled but checkpointing path is not set properly");
        Utils.validateArg(!fourierRegularizationEnabled(), "Fourier regularization is not properly" +
                " implemented yet");
    }
}
