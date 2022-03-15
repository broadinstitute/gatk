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

public final class BGMMVariantAnnotationsModel implements VariantAnnotationsModel {

    private static final Logger logger = LogManager.getLogger(BGMMVariantAnnotationsModel.class);
    private static final int RANDOM_SEED_FOR_WARM_START_SUBSAMPLING = 0;

    public static final String BGMM_FIT_HDF5_SUFFIX = ".bgmmFit.hdf5";
    public static final String BGMM_SCORER_SER_SUFFIX = ".bgmmScorer.ser";

    private final Hyperparameters hyperparameters;

    public BGMMVariantAnnotationsModel(final File hyperparametersJSONFile) {
        hyperparameters = Hyperparameters.readHyperparameters(hyperparametersJSONFile);
    }

    @Override
    public void trainAndSerialize(final File trainingAnnotationsFile,
                                  final String outputPrefix) {
        // load annotation training data
        // TODO validate annotation names and trainingData size
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

        if (!hyperparameters.warmStart || hyperparameters.warmStartSubsample == 0. || hyperparameters.warmStartSubsample == 1.) {
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
        // TODO early fail if warmStart == false and warmStartSubsample nontrivial

        // serialize scorer = preprocesser + BGMM
        // TODO fix up output paths and validation, logging
        final BGMMVariantAnnotationsScorer scorer = new BGMMVariantAnnotationsScorer(annotationNames, preprocesser, bgmm);
        scorer.serialize(new File(outputPrefix + BGMM_SCORER_SER_SUFFIX));

        // write model fit to HDF5
        // TODO fix up output paths and validation, logging
        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();
        fit.write(new File(outputPrefix + BGMM_FIT_HDF5_SUFFIX), "/bgmm");
    }

    // TODO document differences from the sklearn API
    // mean and covariance priors are each specified by a single number here;
    // init_params and warm_start_subsample also require specific documentation
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
        @JsonProperty("warm_start_subsample")
        double warmStartSubsample = 1.;

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
                throw new UserException.BadInput("Could not read hyperparameters JSON.", e);
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