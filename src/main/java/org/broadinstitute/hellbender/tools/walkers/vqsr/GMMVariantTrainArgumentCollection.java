package org.broadinstitute.hellbender.tools.walkers.vqsr;


import com.fasterxml.jackson.annotation.JsonAlias;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixture;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/*
 * A bunch of arguments for GMMVariantTrain.
 * TODO expose number of K-means iterations and symmetry threshold in BGMM
 */
final class GMMVariantTrainArgumentCollection {

    private static final Set<String> ALLOWED_HYPERPARAMETERS = new HashSet<>(Arrays.asList(
            "n_components", "tol", "reg_covar", "max_iter", "n_init", "init_params",
            "weight_concentration_prior", "mean_precision_prior", "degrees_of_freedom_prior", "mean_prior", "covariance_prior",
            "random_state", "warm_start", "verbose_interval"));

    /**
     * JSON file specifying GMM hyperparameters using a subset of the sklearn argument names; see @link #ALLOWED_HYPERPARAMETERS}
     * and https://scikit-learn.org/1.0/modules/generated/sklearn.mixture.BayesianGaussianMixture.html for argument
     * descriptions and types. A few special cases:
     *
     * init_params: "kmeans++", "random"
     * weight_concentration_prior: if unspecified, 1 / n_components will ultimately be used
     * degrees_of_freedom_prior: if unspecified, the number of features will ultimately be used
     * mean_prior: if a valid double, the prior will be taken as this value multiplied by the unit vector; if unspecified, the empirical mean will ultimately be used
     * covariance_prior: if a positive double, the prior will be taken as this value multiplied by the identity matrix; if unspecified, the empirical covariance will ultimately be used
     */
    @Argument(fullName="hyperparameters-json",
            doc="JSON file containing hyperparameters for the GMM.")
    public File hyperparametersJSONFile;

    /**
     * Use either SNP for recalibrating only SNPs (emitting indels untouched in the output VCF) or INDEL for indels (emitting SNPs untouched in the output VCF). There is also a BOTH option for recalibrating both SNPs and indels simultaneously, but this is meant for testing purposes only and should not be used in actual analyses.
     */
    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ", optional = false)
    public VariantRecalibratorArgumentCollection.Mode MODE = VariantRecalibratorArgumentCollection.Mode.SNP;

    /**
     * Generate a VQSR model using per-allele data instead of the default per-site data, assuming that the input VCF contains allele-specific annotations.
     * Annotations should be specified using their full names with AS_ prefix. Non-allele-specific (scalar) annotations will be applied to all alleles.
     */
    @Argument(fullName="use-allele-specific-annotations",
            shortName="AS",
            doc="If specified, the variant recalibrator will attempt to use the allele-specific versions of the specified annotations.", optional=true)
    public boolean useASannotations = false;

    /**
     * If a variant has annotations more than -std standard deviations away from mean, it won't be used for building
     * the Gaussian mixture model.
     */
    @Advanced
    @Argument(fullName = "standard-deviation-threshold", shortName = "std", doc = "Annotation value divergence threshold (number of standard deviations from the means) ", optional = true)
    public double STD_THRESHOLD = 10.0;

    /**
     * The number of variants to use in building the Gaussian mixture model. Training sets larger than this will be randomly downsampled.
     */
    @Advanced
    @Argument(fullName = "maximum-training-variants", doc = "Maximum number of training data", optional = true)
    protected int MAX_NUM_TRAINING_DATA = 2500000;

    /**
     * MQ is capped at a "max" value (60 for bwa-mem) when the alignment is considered perfect. Typically, a huge
     * proportion of the reads in a dataset are perfectly mapped, which yields a distribution of MQ values with a
     * blob below the max value and a huge peak at the max value. This does not conform to the expectations of the
     * Gaussian mixture model of VQSR and has been observed to yield a ROC curve with a jump.
     *
     * This argument aims to mitigate this problem. Using MQCap = X has 2 effects:  (1) MQs are transformed by a scaled
     * logit on [0,X] (+ epsilon to avoid division by zero) to make the blob more Gaussian-like and (2) the transformed
     * MQ=X are jittered to break the peak into a narrow Gaussian.
     *
     * If this option is not used, or if MQCap is set to 0, MQ will not be transformed.
     * Note that AS_MQ uses --mq-jitter instead to control artificial variance
     */
    @Advanced
    @Argument(fullName="mq-cap-for-logit-jitter-transform", shortName = "mq-cap", doc="Apply logit transform and jitter to MQ values", optional=true)
    public int MQ_CAP = 0;

    /**
     * Amount of jitter (as a multiplier to a Normal(0,1) distribution) to add to the AS_MQ and transformed MQ values
     */
    @Advanced
    @Argument(fullName="mq-jitter", doc="Amount of jitter (as a multiplier to a Normal(0,1) distribution) to add to the AS_MQ and transformed MQ values", optional = true)
    public double MQ_JITTER = 0.05;

    @Advanced
    @Argument(fullName = "debug-stdev-thresholding", doc="Output variants that fail standard deviation thresholding to the log for debugging purposes. Redirection of stdout to a file is recommended.", optional = true)
    public boolean debugStdevThresholding = false;

    Hyperparameters hyperparameters = null;

    void readHyperparameters() {
        IOUtils.canReadFile(hyperparametersJSONFile);
        try {
            final ObjectMapper mapper = new ObjectMapper();
            hyperparameters = mapper.readValue(hyperparametersJSONFile, Hyperparameters.class);
        } catch (final Exception e) {
            throw new UserException.BadInput("Could not read hyperparameters JSON.", e);
        }
    }

    // TODO perhaps just put this in the BGMM classes?
    // there is the complication of specifying mean and covariance priors differently here,
    // as well as the possible desire to have different defaults;
    // would also require slightly more code to invoke a BGMM if we use a builder for these instead
    static final class Hyperparameters {
        @JsonProperty("n_components")
        int nComponents = 6;
        @JsonProperty("tol")
        double tol = 1.;
        @JsonProperty("reg_covar")
        double regCovar = 1E-6;
        @JsonProperty("max_iter")
        int maxIter = 150;
        @JsonProperty("n_init")
        int nInit = 1;
        @JsonProperty("init_params")
        String initMethod = "kmeans++";
        @JsonProperty("weight_concentration_prior")
        Double weightConcentrationPrior = null;
        @JsonProperty("mean_precision_prior")
        double meanPrecisionPrior = 1.;
        @JsonProperty("mean_prior")
        Double meanPrior = null;
        @JsonProperty("degrees_of_freedom_prior")
        Double degreesOfFreedomPrior = null;
        @JsonProperty("covariance_prior")
        Double covariancePrior = null;
        @JsonProperty("random_state")
        int seed = 0;
        @JsonProperty("warm_start")
        boolean warmStart = false;
        @JsonProperty("verbose_interval")
        int verboseInterval = 5;

        Hyperparameters() {}

        BayesianGaussianMixture.InitMethod getInitMethod() {
            switch (initMethod) {
                case "kmeans++":
                    return BayesianGaussianMixture.InitMethod.K_MEANS_PLUS_PLUS;
                case "random":
                    return BayesianGaussianMixture.InitMethod.RANDOM;
                default:
                    throw new UserException.BadInput("Unknown initialization method specified.");
            }
        }

        @Override
        public String toString() {
            return "Hyperparameters{" +
                    "nComponents=" + nComponents +
                    ", tol=" + tol +
                    ", regCovar=" + regCovar +
                    ", maxIter=" + maxIter +
                    ", nInit=" + nInit +
                    ", initMethod='" + initMethod + '\'' +
                    ", weightConcentrationPrior=" + weightConcentrationPrior +
                    ", meanPrecisionPrior=" + meanPrecisionPrior +
                    ", meanPrior=" + meanPrior +
                    ", degreesOfFreedomPrior=" + degreesOfFreedomPrior +
                    ", covariancePrior=" + covariancePrior +
                    ", seed=" + seed +
                    ", warmStart=" + warmStart +
                    ", verboseInterval=" + verboseInterval +
                    '}';
        }
    }
}