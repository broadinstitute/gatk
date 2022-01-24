package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a variational posterior on Gaussian Mixture Model parameters, as given in Bishop 10.55-10.63.
 * Should be treated as immutable.
 */
public final class BayesianGaussianMixtureModelPosterior implements Serializable {

    private static final long serialVersionUID = 1L;

    private static final Logger logger = LogManager.getLogger(BayesianGaussianMixtureModelPosterior.class);

    private static final String WEIGHTS_SUBPATH = "/weights";                                // double array
    private static final String WEIGHT_CONCENTRATION_SUBPATH = "/weight_concentration";      // double array
    private static final String MEAN_PRECISION_SUBPATH = "/mean_precision";                  // double array
    private static final String MEANS_SUBPATH = "/means";                                    // double arrays in means/0, means/1, etc. subpaths for each component
    private static final String PRECISIONS_CHOLESKY_SUBPATH = "/precisions_cholesky";        // double matrices in precisions_cholesky/0, precisions_cholesky/1, etc. subpaths for each component
    private static final String COVARIANCES_SUBPATH = "/covariances";                        // double matrices in covariances/0, covariances/1, etc. subpaths for each component
    private static final String DEGREES_OF_FREEDOM_SUBPATH = "/degrees_of_freedom";          // double array
    private static final String NUMBER_OF_COMPONENTS_SUBPATH = "/number_of_components";      // double (can't write int with our HDF5 library)
    private static final String NUMBER_OF_FEATURES_SUBPATH = "/number_of_features";          // double (can't write int with our HDF5 library)

    private final RealVector weightConcentration;
    private final RealVector meanPrecision;
    private final List<RealVector> means;
    private final List<RealMatrix> precisionsCholesky;
    private final List<RealMatrix> covariances;
    private final RealVector degreesOfFreedom;
    private final int nComponents;
    private final int nFeatures;

    public BayesianGaussianMixtureModelPosterior(final RealVector weightConcentration,
                                                 final RealVector meanPrecision,
                                                 final List<RealVector> means,
                                                 final List<RealMatrix> precisionsCholesky,
                                                 final List<RealMatrix> covariances,
                                                 final RealVector degreesOfFreedom) {
        // TODO validation
        this.weightConcentration = weightConcentration;
        this.meanPrecision = meanPrecision;
        this.means = means;
        this.precisionsCholesky = precisionsCholesky;
        this.covariances = covariances;
        this.degreesOfFreedom = degreesOfFreedom;
        nComponents = weightConcentration.getDimension();
        nFeatures = means.get(0).getDimension();
    }

    public RealVector getWeights() {
        return weightConcentration.copy().mapDivideToSelf(BayesianGaussianMixtureUtils.sum(weightConcentration));
    }

    public RealVector getWeightConcentration() {
        return weightConcentration.copy();
    }

    public RealVector getMeanPrecision() {
        return meanPrecision.copy();
    }

    public List<RealVector> getMeans() {
        return means.stream().map(RealVector::copy).collect(Collectors.toList());
    }

    public List<RealMatrix> getPrecisionsCholesky() {
        return precisionsCholesky.stream().map(RealMatrix::copy).collect(Collectors.toList());
    }

    public List<RealMatrix> getCovariances() {
        return covariances.stream().map(RealMatrix::copy).collect(Collectors.toList());
    }

    public RealVector getDegreesOfFreedom() {
        return degreesOfFreedom.copy();
    }

    public int getNumberOfComponents() {
        return nComponents;
    }

    public int getNumberOfFeatures() {
        return nFeatures;
    }

    public static BayesianGaussianMixtureModelPosterior read(final File file,
                                                             final String path) {
        IOUtils.canReadFile(file);
        Utils.validateArg(path.startsWith("/"), "The path should start with / and specify the path within the HDF5 file for the model.");

        try (final HDF5File hdf5File = new HDF5File(file, HDF5File.OpenMode.READ_ONLY)) {
            final int nComponents = (int) hdf5File.readDouble(path + NUMBER_OF_COMPONENTS_SUBPATH);
            final RealVector weightConcentration = new ArrayRealVector(hdf5File.readDoubleArray(path + WEIGHT_CONCENTRATION_SUBPATH));
            final RealVector meanPrecision = new ArrayRealVector(hdf5File.readDoubleArray(path + MEAN_PRECISION_SUBPATH));
            final List<RealVector> means = IntStream.range(0, nComponents).boxed()
                    .map(k -> new ArrayRealVector(hdf5File.readDoubleArray(path + MEANS_SUBPATH + '/' + k)))
                    .collect(Collectors.toList());
            final List<RealMatrix> precisionsCholesky = IntStream.range(0, nComponents).boxed()
                    .map(k -> new Array2DRowRealMatrix(hdf5File.readDoubleMatrix(path + PRECISIONS_CHOLESKY_SUBPATH + '/' + k)))
                    .collect(Collectors.toList());
            final List<RealMatrix> covariances = IntStream.range(0, nComponents).boxed()
                    .map(k -> new Array2DRowRealMatrix(hdf5File.readDoubleMatrix(path + COVARIANCES_SUBPATH + '/' + k)))
                    .collect(Collectors.toList());
            final RealVector degreesOfFreedom = new ArrayRealVector(hdf5File.readDoubleArray(path + DEGREES_OF_FREEDOM_SUBPATH));

            // TODO validate nComponents and nFeatures

            return new BayesianGaussianMixtureModelPosterior(
                    weightConcentration,
                    meanPrecision,
                    means,
                    precisionsCholesky,
                    covariances,
                    degreesOfFreedom);
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during reading of BayesianGaussianMixtureModelPosterior (%s) from the input file at %s.",
                    exception, file.getAbsolutePath()));
        }
    }

    public void write(final File file,
                      final String path) {
        try (final HDF5File hdf5File = new HDF5File(file, HDF5File.OpenMode.CREATE)) { // TODO allow appending
            IOUtils.canReadFile(hdf5File.getFile());
            Utils.validateArg(path.startsWith("/"), "The path should start with / and specify the path within the HDF5 file for the model.");

            hdf5File.makeDoubleArray(path + WEIGHTS_SUBPATH, getWeights().toArray());
            hdf5File.makeDoubleArray(path + WEIGHT_CONCENTRATION_SUBPATH, weightConcentration.toArray());
            hdf5File.makeDoubleArray(path + MEAN_PRECISION_SUBPATH, meanPrecision.toArray());
            IntStream.range(0, nComponents).forEach(
                    k -> hdf5File.makeDoubleArray(path + MEANS_SUBPATH + '/' + k, means.get(k).toArray()));
            IntStream.range(0, nComponents).forEach(
                    k -> hdf5File.makeDoubleMatrix(path + PRECISIONS_CHOLESKY_SUBPATH + '/' + k, precisionsCholesky.get(k).getData()));
            IntStream.range(0, nComponents).forEach(
                    k -> hdf5File.makeDoubleMatrix(path + COVARIANCES_SUBPATH + '/' + k, covariances.get(k).getData()));
            hdf5File.makeDoubleArray(path + DEGREES_OF_FREEDOM_SUBPATH, degreesOfFreedom.toArray());
            hdf5File.makeDouble(path + NUMBER_OF_COMPONENTS_SUBPATH, nComponents);
            hdf5File.makeDouble(path + NUMBER_OF_FEATURES_SUBPATH, nFeatures);

            logger.info(String.format("BayesianGaussianMixtureModelPosterior written to %s.", file.getAbsolutePath()));
        } catch (final RuntimeException exception) {
            throw new GATKException(String.format("Exception encountered during writing of BayesianGaussianMixtureModelPosterior (%s). Output file at %s may be in a bad state.",
                    exception, file.getAbsolutePath()));
        }
    }
}