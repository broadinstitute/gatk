package org.broadinstitute.hellbender.utils.clustering;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a variational posterior on Gaussian Mixture Model parameters, as given in Bishop 10.55-10.63.
 * Should be treated as immutable.
 */
public final class BayesianGaussianMixtureModelPosterior {

    private final RealVector weightConcentration;
    private final RealVector meanPrecision;
    private final List<RealVector> means;
    private final List<RealMatrix> precisionsCholesky;
    private final List<RealMatrix> covariances;
    private final RealVector degreesOfFreedom;

    public BayesianGaussianMixtureModelPosterior(final RealVector weightConcentration,
                                                 final RealVector meanPrecision,
                                                 final List<RealVector> means,
                                                 final List<RealMatrix> precisionsCholesky,
                                                 final List<RealMatrix> covariances,
                                                 final RealVector degreesOfFreedom) {
        this.weightConcentration = weightConcentration;
        this.meanPrecision = meanPrecision;
        this.means = means;
        this.precisionsCholesky = precisionsCholesky;
        this.covariances = covariances;
        this.degreesOfFreedom = degreesOfFreedom;
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

    // TODO static read/write HDF5
}