package org.broadinstitute.hellbender.tools.coveragemodel;

import java.io.Serializable;

/**
 * Global constants for coverage model package classes
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelGlobalConstants implements Serializable {

    private static final long serialVersionUID = 3185225539827967119L;

    /**
     * This value is arbitrary and immaterial
     */
    public static final double POISSON_STATISTICAL_VARIANCE_ON_MASKED_TARGETS = 1.0;

    /**
     * This value is arbitrary and immaterial
     */
    public static final double LOG_READ_COUNT_ON_MASKED_TARGETS = 0.0;

    /**
     * This value is arbitrary and immaterial
     */
    public static final double MEAN_LOG_COPY_RATIO_ON_MASKED_TARGETS = 0.0;

    /**
     * This value is arbitrary and immaterial
     */
    public static final double VAR_LOG_COPY_RATIO_ON_MASKED_TARGETS = 0.0;

    /**
     * Function evaluation accuracy (used in various root finders)
     */
    public static final double DEFAULT_FUNCTION_EVALUATION_ACCURACY = 1e-12;

    /**
     * Read counts on zero ploidy targets (such as targets on Y contig on XX samples) will be replaced by the
     * following value, regardless of the value provided by the user
     */
    public static final int READ_COUNT_ON_ZERO_PLOIDY_TARGETS = 0;

    /**
     * Minimum size of a target block
     */
    public static final int DEFAULT_MIN_TARGET_BLOCK_SIZE = 5;

    /**
     * Uniform random maximum for target-specific unexplained variance
     */
    public static final double RANDOM_UNEXPLAINED_VARIANCE_MAX = 0.0;

    /**
     * TODO -- either estimate this from data for expose to user
     *
     * Gaussian random standard deviation for bias covariates
     */
    public static final double RANDOM_BIAS_COVARIATES_STD = 0.1;

    /**
     * Gaussian random standard deviation for mean log bias
     */
    public static final double RANDOM_MEAN_LOG_BIAS_STD = 0.0;

    /**
     * Random generator seed for initial model parameters
     */
    public static final long RANDOM_MODEL_SEED = 1984;

    /**
     * Copy ratio max likelihood estimates output file name
     */
    public static final String COPY_RATIO_MAX_LIKELIHOOD_ESTIMATES_FILENAME = "copy_ratio_max_likelihood_estimate_matrix.tsv";

    /**
     * Copy ratio precision (= inverse total unexplained variance) output file name
     */
    public static final String COPY_RATIO_PRECISION_FILENAME = "copy_ratio_precision_matrix.tsv";

    /**
     * Copy ratio Viterbi hidden state chains output file name
     */
    public static final String COPY_RATIO_VITERBI_FILENAME = "copy_ratio_Viterbi_matrix.tsv";

    /**
     * Sample read depth posteriors output file name
     */
    public static final String SAMPLE_READ_DEPTH_POSTERIORS_FILENAME = "sample_read_depth_posteriors.tsv";

    /**
     * Sample model log likelihood output file name
     */
    public static final String SAMPLE_LOG_LIKELIHOODS_FILENAME = "sample_log_likelihoods.tsv";

    /**
     * Sample-specific unexplained variance output file name
     */
    public static final String SAMPLE_UNEXPLAINED_VARIANCE_FILENAME = "sample_specific_unexplained_variance.tsv";

    /**
     * Bias latent indicator variables posteriors output file name
     */
    public static final String SAMPLE_BIAS_LATENT_POSTERIORS_FILENAME = "sample_bias_latent_posteriors.tsv";

    /**
     * Total unexplained variance (sum of target- and sample-specific) output file name
     */
    public static final String TOTAL_UNEXPLAINED_VARIANCE_FILENAME = "total_unexplained_variance_matrix.tsv";

    /**
     * Total covariate bias (= \sum_{\mu} W_{t\mu} E[z_{s\mu}]) output file name
     */
    public static final String TOTAL_COVARIATE_BIAS_FILENAME = "total_covariate_bias_matrix.tsv";

    /**
     * Target-specific mean log bias (= m_t) output file name
     */
    public static final String TARGET_MEAN_LOG_BIAS_OUTPUT_FILE = "target_specific_mean_log_bias.tsv";

    /**
     * Target-specific unexplained variance (= \Psi_{t}) output file name
     */
    public static final String TARGET_UNEXPLAINED_VARIANCE_OUTPUT_FILE = "target_specific_unexplained_variance.tsv";

    /**
     * Bias covariates (= W_{t\mu} = "log bias principal components") output file name
     */
    public static final String MEAN_BIAS_COVARIATES_OUTPUT_FILE = "mean_bias_covariates_matrix.tsv";

    /**
     * Norm_2 of bias covariates (= \sum_t |W_{t\mu}|^2) output file name
     */
    public static final String MEAN_BIAS_COVARIATES_NORM2_OUTPUT_FILE = "mean_bias_covariates_norm2.tsv";

    /**
     * ARD coefficients of bias covariates output file
     */
    public static final String BIAS_COVARIATES_ARD_COEFFICIENTS_OUTPUT_FILE = "bias_covariates_ARD_coefficients.tsv";

    /**
     * History of ARD coefficients of bias covariates output file
     */
    public static final String BIAS_COVARIATES_ARD_COEFFICIENTS_HISTORY_OUTPUT_FILE = "bias_covariates_ARD_coefficients_history.tsv";

    /**
     * Log likelihood history output file
     */
    public static final String LOG_LIKELIHOODS_HISTORY_FILENAME = "log_likelihood_history.tsv";

    /**
     * Processed targets output file name (order may be different than original, and some may have been dropped)
     */
    public static final String TARGET_LIST_OUTPUT_FILE = "targets.tsv";

    /**
     * Copy ratio segments subdirectory
     */
    public static final String COPY_RATIO_SEGMENTS_SUBDIR = "segments";

    /**
     * Prefix for posteriors checkpointing output directories
     */
    public static final String POSTERIOR_CHECKPOINT_PATH_PREFIX = "posteriors_checkpoint";

    /**
     * Prefix for model checkpointing output directories
     */
    public static final String MODEL_CHECKPOINT_PATH_PREFIX = "model_checkpoint";
}
