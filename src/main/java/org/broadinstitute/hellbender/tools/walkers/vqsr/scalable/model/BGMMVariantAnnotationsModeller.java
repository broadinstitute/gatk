package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.model;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.collect.ImmutableList;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.VariantAnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModelPosterior;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

final class BGMMVariantAnnotationsModeller implements VariantAnnotationsModeller {

    private static final Logger logger = LogManager.getLogger(BGMMVariantAnnotationsModeller.class);

    private static final String SCORER_SER_SUFFIX = ".scorer.ser";
    private static final String BGMM_HDF5_SUFFIX = ".bgmm.hdf5";

    private final List<String> annotationNames;
    private final BayesianGaussianMixtureModeller bgmm;

    public BGMMVariantAnnotationsModeller(final File hyperparametersJSONFile,
                                          final File inputAnnotationsFile) {
        final Hyperparameters hyperparameters = Hyperparameters.readHyperparameters(hyperparametersJSONFile);

        annotationNames = ImmutableList.copyOf(LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile));

        final int nFeatures = annotationNames.size();

        final double[][] covariancePrior = MatrixUtils.createRealIdentityMatrix(nFeatures).getData();
        bgmm = new BayesianGaussianMixtureModeller.Builder()
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
        logger.info(String.format("Constructed BayesianGaussianMixtureModeller for %d annotations...", nFeatures));
    }

    @Override
    public List<String> getAnnotationNames() {
        return annotationNames;
    }

    @Override
    public void fitAndWriteTrainingAndTruthScoresAndSerialize(final File inputAnnotationsFile,
                                                              final String outputPrefix) {
        // load annotation training data
        logger.info("Reading annotations...");
        // TODO validate annotation names
        final double[][] allData = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);
        final List<Boolean> isTraining = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRAINING_LABEL);
        final List<Boolean> isTruth = LabeledVariantAnnotationsData.readLabel(inputAnnotationsFile, LabeledVariantAnnotationsData.TRUTH_LABEL); // TODO make truth optional
        final double[][] trainingData = IntStream.range(0, isTraining.size()).boxed()
                .filter(isTraining::get)
                .map(i -> allData[i])
                .toArray(double[][]::new);
        final double[][] truthData = IntStream.range(0, isTruth.size()).boxed()
                .filter(isTruth::get)
                .map(i -> allData[i])
                .toArray(double[][]::new);

        logger.info(String.format("Training BayesianGaussianMixtureModeller with %d training sites x %d annotations...", nSamples, nFeatures));

        // preprocess
        final VariantAnnotationUtils.Preprocesser preprocesser = new VariantAnnotationUtils.Preprocesser();
        final double[][] preprocessedTrainingData = preprocesser.fitTransform(trainingData);
        final double[][] preprocessedTruthData = preprocesser.transform(truthData);

        bgmm.fit(preprocessedTrainingData);

        // TODO clean this up, method for inputAnnotationsFile -> outputPreprocessedAnnotationsFile
//        // write preprocessed annotations
//        final File outputPreprocessedAnnotationsFile = new File(outputPrefix + ".annot.pre.hdf5");
//        try (final HDF5File hdf5File = new HDF5File(outputPreprocessedAnnotationsFile, HDF5File.OpenMode.CREATE)) { // TODO allow appending
//            IOUtils.canReadFile(hdf5File.getFile());
//
//            hdf5File.makeStringArray("/annotations/names", annotationNames);
//            HDF5Utils.writeChunkedDoubleMatrix(hdf5File, "/annotations", preprocessedTrainingData, maximumChunkSize);
//            hdf5File.makeDoubleArray("/labels/training", isTraining.stream().mapToDouble(x -> x ? 1 : 0).toArray());
//            hdf5File.makeDoubleArray("/labels/truth", isTraining.stream().mapToDouble(x -> x ? 1 : 0).toArray());
//        } catch (final HDF5LibException exception) {
//            throw new GATKException(String.format("Exception encountered during writing of preprocessed annotations (%s). Output file at %s may be in a bad state.",
//                    exception, outputPreprocessedAnnotationsFile.getAbsolutePath()));
//        }
//        logger.info(String.format("Preprocessed annotations written to %s.", outputPreprocessedAnnotationsFile.getAbsolutePath()));

        // serialize scorer = preprocesser + BGMM
        // TODO fix up output paths and validation
        VariantAnnotationUtils.serialize(
                new File(outputPrefix + SCORER_SER_SUFFIX), new VariantAnnotationUtils.Scorer(preprocesser, bgmm)
        );
        final BayesianGaussianMixtureModelPosterior fit = bgmm.getBestFit();
        fit.write(new File(outputPrefix + BGMM_HDF5_SUFFIX), "/bgmm");

        // generate scores and write to HDF5
        final double[] trainingScores = bgmm.scoreSamples(preprocessedTrainingData);
        VariantAnnotationUtils.writeScores(outputTrainingScoresFile, trainingScores);

        final double[] truthScores = bgmm.scoreSamples(preprocessedTruthData);
        VariantAnnotationUtils.writeScores(outputTruthScoresFile, truthScores);
    }

    @Override
    public void scoreSamplesAndWriteScores(final String inputModelPrefix, final File outputScoresFile) {

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

        Hyperparameters() {
        }

        static Hyperparameters readHyperparameters(final File hyperparametersJSONFile) {
            IOUtils.canReadFile(hyperparametersJSONFile);
            try {
                final ObjectMapper mapper = new ObjectMapper();
                return mapper.readValue(hyperparametersJSONFile, Hyperparameters.class);
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
}
