package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.TrainVariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModelPosterior;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class wraps the {@link BayesianGaussianMixtureModeller} and adds warm-start functionality, while adding
 * further restrictions to the API. See documentation in {@link BGMMVariantAnnotationsModel.Hyperparameters}.
 * Note also that hyperparameters for regularization are not exposed and the default values are used.
 *
 * All methods assume valid input (all annotation names and dimensions correct, missing values imputed, etc.)
 * has been provided by calling code in {@link TrainVariantAnnotationsModel}.
 */
public final class BGMMVariantAnnotationsModel implements VariantAnnotationsModel {

    private static final Logger logger = LogManager.getLogger(BGMMVariantAnnotationsModel.class);
    private static final int RANDOM_SEED_FOR_WARM_START_SUBSAMPLING = 0;

    public static final String BGMM_FIT_HDF5_SUFFIX = ".bgmmFit.hdf5";

    private final Hyperparameters hyperparameters;

    public BGMMVariantAnnotationsModel(final File hyperparametersJSONFile) {
        hyperparameters = Hyperparameters.readHyperparameters(hyperparametersJSONFile);
    }

    @Override
    public void trainAndSerialize(final File trainingAnnotationsFile,
                                  final File unlabeledAnnotationsFile,
                                  final String outputPrefix) {
        // load annotation training data
        final List<String> annotationNames = LabeledVariantAnnotationsData.readAnnotationNames(trainingAnnotationsFile);
        final double[][] trainingAnnotations = LabeledVariantAnnotationsData.readAnnotations(trainingAnnotationsFile);

        final int nSamples = trainingAnnotations.length;
        final int nFeatures = trainingAnnotations[0].length;

        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();
        final BayesianGaussianMixtureModeller bgmm = new BayesianGaussianMixtureModeller.Builder()
                .nComponents(hyperparameters.nComponents)
                .tol(hyperparameters.tol)
                .regCovar(hyperparameters.regCovar)
                .maxIter(hyperparameters.maxIter)
                .nInit(hyperparameters.nInit)
                .initMethod(hyperparameters.getInitMethod())
                .weightConcentrationPrior(hyperparameters.weightConcentrationPrior)
                .meanPrecisionPrior(hyperparameters.meanPrecisionPrior)
                .meanPrior(
                        hyperparameters.meanPrior != null
                                ? Collections.nCopies(nFeatures, hyperparameters.meanPrior).stream().mapToDouble(x -> x).toArray()
                                : null)
                .degreesOfFreedomPrior(hyperparameters.degreesOfFreedomPrior)
                .covariancePrior(covariancePrior)
                .seed(hyperparameters.seed)
                .warmStart(hyperparameters.warmStart)
                .verboseInterval(hyperparameters.verboseInterval)
                .build();

        // preprocess
        final Preprocesser preprocesser = new Preprocesser();
        final double[][] preprocessedTrainingAnnotations = preprocesser.fitTransform(trainingAnnotations);

        if (!hyperparameters.warmStart && hyperparameters.warmStartSubsample != 0. && hyperparameters.warmStartSubsample != 1.) {
            throw new UserException.BadInput("Warm start was not specified, but a nontrivial value for the warm-start subsample fraction was specified.");
        } else if (!hyperparameters.warmStart || hyperparameters.warmStartSubsample == 0. || hyperparameters.warmStartSubsample == 1.) {
            bgmm.fit(preprocessedTrainingAnnotations);
        } else {
            final int nSubsamples = (int) (hyperparameters.warmStartSubsample * nSamples);
            logger.info(String.format("Performing warm start with a subsample fraction (%d / %d) of the data...",
                    nSubsamples, nSamples));
            final List<Integer> shuffledIndices = IntStream.range(0, nSamples).boxed().collect(Collectors.toList());
            Collections.shuffle(shuffledIndices, new Random(RANDOM_SEED_FOR_WARM_START_SUBSAMPLING));
            final double[][] preprocessedTrainingAnnotationsSubsample =
                    shuffledIndices.subList(0, nSubsamples).stream()
                            .map(i -> preprocessedTrainingAnnotations[i])
                            .toArray(double[][]::new);
            bgmm.fit(preprocessedTrainingAnnotationsSubsample); // perform hyperparameters.nInit initializations and fits with subsample
            bgmm.fit(preprocessedTrainingAnnotations);          // perform final fit with full sample
        }

        // serialize scorer = preprocesser + BGMM
        final BGMMVariantAnnotationsScorer scorer = new BGMMVariantAnnotationsScorer(annotationNames, preprocesser, bgmm);
        scorer.serialize(new File(outputPrefix + BGMMVariantAnnotationsScorer.BGMM_SCORER_SER_SUFFIX));

        // write model fit to HDF5
        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();
        fit.write(new File(outputPrefix + BGMM_FIT_HDF5_SUFFIX));
        // TODO append standardization fields to this HDF5 file?
    }

    /**
     * Hyperparameters closely mirror those for the scikit-learn BayesianGaussianMixture module
     * (see https://scikit-learn.org/1.0/modules/generated/sklearn.mixture.BayesianGaussianMixture.html),
     * with the following differences:
     *
     *  1) Only full covariance matrices are allowed (covariance_type = 'full' is fixed),
     *  2) Responsibilities can be initialized using k-means++ clustering or at random using the parameter
     *     init_params: {'kmeans++', 'random'}; in contrast, the scikit-learn module allows initialization using
     *     standard k-means clustering or at random,
     *  3) Only Dirichlet weight concentration priors are allowed (weight_concentration_prior_type = 'dirichlet_distribution' is fixed),
     *  4) The mean_prior parameter is specified by a single float value, which is multiplied by the unit array to give the array-like prior,
     *  5) The covariance_prior parameter is specified by a single float value, which is multiplied by the identity matrix to give the array-like prior,
     *  6) The null value should be used in place of None where appropriate,
     *  7) Boolean values should not be capitalized,
     *  8) An additional parameter warm_start_subsample (float in [0, 1]) is used to specify the fraction of data that
     *  will be used for warm-start fitting; a random subsample will be used to generate n_init initial fits, of which
     *  the best will be selected and used to initialize a final fit on the full data.
     *
     *  This documentation is replicated in the default JSON at
     *  src/main/resources/org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/bgmm-hyperparameters.json.
     */
    static final class Hyperparameters {
        @JsonProperty("comment")
        String comment;
        @JsonProperty("n_components")
        int nComponents;
        @JsonProperty("tol")
        double tol;
        @JsonProperty("reg_covar")
        double regCovar;
        @JsonProperty("max_iter")
        int maxIter;
        @JsonProperty("n_init")
        int nInit;
        @JsonProperty("init_params")
        String initMethod = "kmeans++";
        @JsonProperty("weight_concentration_prior")
        Double weightConcentrationPrior;
        @JsonProperty("mean_precision_prior")
        double meanPrecisionPrior;
        @JsonProperty("mean_prior")
        Double meanPrior;
        @JsonProperty("degrees_of_freedom_prior")
        Double degreesOfFreedomPrior;
        @JsonProperty("covariance_prior")
        Double covariancePrior;
        @JsonProperty("random_state")
        int seed;
        @JsonProperty("warm_start")
        boolean warmStart;
        @JsonProperty("verbose_interval")
        int verboseInterval;
        @JsonProperty("warm_start_subsample")
        double warmStartSubsample;

        Hyperparameters() {
        }

        static Hyperparameters readHyperparameters(final File hyperparametersJSONFile) {
            IOUtils.canReadFile(hyperparametersJSONFile);
            try {
                final ObjectMapper mapper = new ObjectMapper();
                final Hyperparameters hyperparameters = mapper.readValue(hyperparametersJSONFile, Hyperparameters.class);
                ParamUtils.inRange(hyperparameters.warmStartSubsample, 0., 1.,
                        "The hyperparameter warm_start_subsample must be in [0, 1].");
                return hyperparameters;
            } catch (final Exception e) {
                throw new UserException.BadInput(String.format("Could not read hyperparameters JSON: %s.", e));
            }
        }

        BayesianGaussianMixtureModeller.InitMethod getInitMethod() {
            switch (initMethod) {
                case "kmeans++":
                    return BayesianGaussianMixtureModeller.InitMethod.K_MEANS_PLUS_PLUS;
                case "random":
                    return BayesianGaussianMixtureModeller.InitMethod.RANDOM;
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
                    ", warmStartSubsample=" + warmStartSubsample +
                    '}';
        }
    }

    // TODO tighten this up?
    static final class Preprocesser implements Serializable {
        private static final long serialVersionUID = 1L;

        private double[] meansForStandardization;
        private double[] standardDeviationsForStandardization;
        private double[] standardizedMediansForImputation;

        Preprocesser() {
        }

        private static double standardize(final double value,
                                          final double mean,
                                          final double standardDeviation) {
            // if standard deviation is zero, no standardization is needed, but we must guard against divide-by-zero
            return (value - mean) / (standardDeviation < Double.MIN_VALUE ? 1. : standardDeviation);
        }

        double[][] fitTransform(final double[][] data) {
            // TODO validation
            final int nSamples = data.length;
            final int nFeatures = data[0].length;
            meansForStandardization = IntStream.range(0, nFeatures)
                    .mapToDouble(j -> new Mean().evaluate(
                            IntStream.range(0, nSamples).boxed()
                                    .mapToDouble(i -> data[i][j])
                                    .filter(Double::isFinite)
                                    .toArray()))
                    .toArray();
            standardDeviationsForStandardization = IntStream.range(0, nFeatures)
                    .mapToDouble(j -> new Variance().evaluate(
                            IntStream.range(0, nSamples).boxed()
                                    .mapToDouble(i -> data[i][j])
                                    .filter(Double::isFinite)
                                    .toArray()))
                    .map(Math::sqrt)
                    .toArray();

            standardizedMediansForImputation = IntStream.range(0, nFeatures)
                    .mapToDouble(j -> new Median().evaluate(
                            IntStream.range(0, nSamples).boxed()
                                    .mapToDouble(i -> data[i][j])
                                    .filter(Double::isFinite)
                                    .map(x -> standardize(x, meansForStandardization[j], standardDeviationsForStandardization[j]))
                                    .toArray()))
                    .toArray();

            return transform(data);
        }

        double[][] transform(final double[][] data) {
            // TODO validation
            final int nSamples = data.length;
            final int nFeatures = data[0].length;
            double[][] preprocessedData = new double[nSamples][nFeatures];
            for (int i = 0; i < nSamples; i++) {
                for (int j = 0; j < nFeatures; j++) {
                    final double value = data[i][j];
                    final double preprocessedValue = Double.isNaN(value)
                            ? standardizedMediansForImputation[j]
                            : standardize(value, meansForStandardization[j], standardDeviationsForStandardization[j]);
                    if (!Double.isFinite(preprocessedValue)) {
                        throw new GATKException("Encountered non-finite value during preprocessing of annotations.");
                    }
                    preprocessedData[i][j] = preprocessedValue;
                }
            }
            return preprocessedData;
        }
    }
}