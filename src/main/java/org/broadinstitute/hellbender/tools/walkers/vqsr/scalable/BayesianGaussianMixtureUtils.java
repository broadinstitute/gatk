package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

final class BayesianGaussianMixtureUtils {

    private static final Logger logger = LogManager.getLogger(BayesianGaussianMixtureUtils.class);

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

    static void serializeBGMM(final BayesianGaussianMixtureModeller bgmm,
                              final File outputFile) {
        try (final FileOutputStream fileOutputStream = new FileOutputStream(outputFile);
             final ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream)) {
            objectOutputStream.writeObject(bgmm);
        } catch (final IOException e) {
            throw new GATKException(String.format("Exception encountered during serialization of BayesianGaussianMixtureModeller to %s: %s",
                    outputFile.getAbsolutePath(), e));
        }
    }

    static BayesianGaussianMixtureModeller deserializeBGMM(final File inputFile) {
        try (final FileInputStream fileInputStream = new FileInputStream(inputFile);
             final ObjectInputStream objectInputStream = new ObjectInputStream(fileInputStream)) {
            return (BayesianGaussianMixtureModeller) objectInputStream.readObject();
        } catch (final IOException | ClassNotFoundException e) {
            throw new GATKException(String.format("Exception encountered during deserialization of BayesianGaussianMixtureModeller from %s: %s",
                    inputFile.getAbsolutePath(), e));
        }
    }
}
