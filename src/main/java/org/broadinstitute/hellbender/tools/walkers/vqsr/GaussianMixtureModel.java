package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import Jama.Matrix;

class GaussianMixtureModel {

    protected final static Logger logger = LogManager.getLogger(GaussianMixtureModel.class);

    private final List<MultivariateGaussian> gaussians;
    private final double shrinkage;
    private final double dirichletParameter;
    private final double priorCounts;
    private final double[] empiricalMu;
    private final Matrix empiricalSigma;
    public boolean isModelReadyForEvaluation;
    public boolean failedToConverge = false;

    public GaussianMixtureModel( final int numGaussians, final int numVariantData, final int numAnnotations,
                                 final double shrinkage, final double dirichletParameter, final double priorCounts ) {

        gaussians = new ArrayList<>( numGaussians );
        for( int iii = 0; iii < numGaussians; iii++ ) {
            final MultivariateGaussian gaussian = new MultivariateGaussian( numVariantData, numAnnotations );
            gaussians.add( gaussian );
        }
        this.shrinkage = shrinkage;
        this.dirichletParameter = dirichletParameter;
        this.priorCounts = priorCounts;
        empiricalMu = new double[numAnnotations];
        empiricalSigma = new Matrix(numAnnotations, numAnnotations);
        isModelReadyForEvaluation = false;
        Arrays.fill(empiricalMu, 0.0);
        empiricalSigma.setMatrix(0, empiricalMu.length - 1, 0, empiricalMu.length - 1, Matrix.identity(empiricalMu.length, empiricalMu.length).times(200.0).inverse());
    }

    //this is used for the model output unit test
    protected GaussianMixtureModel(final List<MultivariateGaussian> gaussians, final double shrinkage, final double dirichletParameter, final double priorCounts ) {
        this.gaussians = gaussians;
        final int numAnnotations = gaussians.get(0).mu.length;
        this.shrinkage = shrinkage;
        this.dirichletParameter = dirichletParameter;
        this.priorCounts = priorCounts;
        empiricalMu = new double[numAnnotations];
        empiricalSigma = new Matrix(numAnnotations, numAnnotations);
        isModelReadyForEvaluation = false;
        Arrays.fill(empiricalMu, 0.0);
        empiricalSigma.setMatrix(0, empiricalMu.length - 1, 0, empiricalMu.length - 1, Matrix.identity(empiricalMu.length, empiricalMu.length).times(200.0).inverse());

    }

    public void initializeRandomModel( final List<VariantDatum> data, final int numKMeansIterations ) {

        // initialize random Gaussian means // BUGBUG: this is broken up this way to match the order of calls to rand.nextDouble() in the old code
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.initializeRandomMu( Utils.getRandomGenerator() );
        }

        // initialize means using K-means algorithm
        logger.info( "Initializing model with " + numKMeansIterations + " k-means iterations..." );
        initializeMeansUsingKMeans( data, numKMeansIterations );

        // initialize uniform mixture coefficients, random covariance matrices, and initial hyperparameters
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.pMixtureLog10 = Math.log10( 1.0 / ((double) gaussians.size()) );
            gaussian.sumProb = 1.0 / ((double) gaussians.size());
            gaussian.initializeRandomSigma( Utils.getRandomGenerator() );
            gaussian.hyperParameter_a = priorCounts;
            gaussian.hyperParameter_b = shrinkage;
            gaussian.hyperParameter_lambda = dirichletParameter;
        }
    }

    private void initializeMeansUsingKMeans( final List<VariantDatum> data, final int numIterations ) {

        int ttt = 0;
        while( ttt++ < numIterations ) {
            // E step: assign each variant to the nearest cluster
            for( final VariantDatum datum : data ) {
                double minDistance = Double.MAX_VALUE;
                MultivariateGaussian minGaussian = null;
                datum.assignment = minGaussian;
                for( final MultivariateGaussian gaussian : gaussians ) {
                    final double dist = gaussian.calculateDistanceFromMeanSquared( datum );
                    if( dist < minDistance ) {
                        minDistance = dist;
                        minGaussian = gaussian;
                    }
                }
                datum.assignment = minGaussian;
            }

            // M step: update gaussian means based on assigned variants
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.zeroOutMu();
                int numAssigned = 0;

                for( final VariantDatum datum : data ) {
                    if( datum.assignment.equals(gaussian) ) {
                        numAssigned++;
                        gaussian.incrementMu( datum );
                    }
                }
                if( numAssigned != 0 ) {
                    gaussian.divideEqualsMu( ((double) numAssigned) );
                } else {
                    gaussian.initializeRandomMu( Utils.getRandomGenerator() );
                }
            }
        }
    }

    public void expectationStep( final List<VariantDatum> data ) {

        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.precomputeDenominatorForVariationalBayes( getSumHyperParameterLambda() );
        }

        for( final VariantDatum datum : data ) {
            final double[] pVarInGaussianLog10 = gaussians.stream().mapToDouble(g -> g.evaluateDatumLog10(datum)).toArray();
            final double[] pVarInGaussianNormalized = MathUtils.normalizeLog10DeleteMePlease( pVarInGaussianLog10, false);
            int gaussianIndex = 0;
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.assignPVarInGaussian( pVarInGaussianNormalized[gaussianIndex++] );
            }
        }
    }

    public void maximizationStep( final List<VariantDatum> data ) {
        gaussians.forEach(g -> g.maximizeGaussian( data, empiricalMu, empiricalSigma, shrinkage, dirichletParameter, priorCounts));
    }

    private double getSumHyperParameterLambda() {
        return gaussians.stream().mapToDouble(g -> g.hyperParameter_lambda).sum();
    }

    public void evaluateFinalModelParameters( final List<VariantDatum> data ) {
        gaussians.forEach(g -> g.evaluateFinalModelParameters(data));
        normalizePMixtureLog10();
    }

    public double normalizePMixtureLog10() {
        double sumDiff = 0.0;
        final double sumPK = gaussians.stream().mapToDouble(g -> g.sumProb).sum();

        final double log10SumPK = Math.log10(sumPK);
        final double[] pGaussianLog10 = gaussians.stream().mapToDouble(g -> Math.log10(g.sumProb) - log10SumPK).toArray();
        MathUtils.normalizeLog10DeleteMePlease( pGaussianLog10, true);

        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sumDiff += Math.abs( pGaussianLog10[gaussianIndex] - gaussian.pMixtureLog10 );
            gaussian.pMixtureLog10 = pGaussianLog10[gaussianIndex++];
        }
        return sumDiff;
    }

    public void precomputeDenominatorForEvaluation() {
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.precomputeDenominatorForEvaluation();
        }

        isModelReadyForEvaluation = true;
    }

    /**
     * A version of Log10SumLog10 that tolerates NaN values in the array
     *
     * In the case where one or more of the values are NaN, this function returns NaN
     *
     * @param values a non-null vector of doubles
     * @return log10 of the sum of the log10 values, or NaN
     */
    private double nanTolerantLog10SumLog10(final double[] values) {
        for ( final double value : values ) {
            if ( Double.isNaN(value) ) {
                return Double.NaN;
            }
        }
        return MathUtils.log10sumLog10(values);
    }

    public double evaluateDatum( final VariantDatum datum ) {
        for( final boolean isNull : datum.isNull ) {
            if( isNull ) {
                return evaluateDatumMarginalized( datum );
            }
        }
        // Fill an array with the log10 probability coming from each Gaussian and then use MathUtils to sum them up correctly
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + gaussian.evaluateDatumLog10( datum );
        }
        return nanTolerantLog10SumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
    }

    // Used only to decide which covariate dimension is most divergent in order to report in the culprit info field annotation
    public Double evaluateDatumInOneDimension( final VariantDatum datum, final int iii ) {
        if(datum.isNull[iii]) { return null; }

        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + MathUtils.normalDistributionLog10(gaussian.mu[iii], gaussian.sigma.get(iii, iii), datum.annotations[iii]);
        }
        return nanTolerantLog10SumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
    }

    public double evaluateDatumMarginalized( final VariantDatum datum ) {
        int numRandomDraws = 0;
        double sumPVarInGaussian = 0.0;
        final int numIterPerMissingAnnotation = 20; // Trade off here between speed of computation and accuracy of the marginalization
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        // for each dimension
        for( int iii = 0; iii < datum.annotations.length; iii++ ) {
            // if it is missing marginalize over the missing dimension by drawing X random values for the missing annotation and averaging the lod
            if( datum.isNull[iii] ) {
                for( int ttt = 0; ttt < numIterPerMissingAnnotation; ttt++ ) {
                    datum.annotations[iii] = Utils.getRandomGenerator().nextGaussian(); // draw a random sample from the standard normal distribution

                    // evaluate this random data point
                    int gaussianIndex = 0;
                    for( final MultivariateGaussian gaussian : gaussians ) {
                        pVarInGaussianLog10[gaussianIndex++] = gaussian.pMixtureLog10 + gaussian.evaluateDatumLog10( datum );
                    }

                    // add this sample's probability to the pile in order to take an average in the end
                    sumPVarInGaussian += Math.pow(10.0, nanTolerantLog10SumLog10(pVarInGaussianLog10)); // p = 10 ^ Sum(pi_k * p(v|n,k))
                    numRandomDraws++;
                }
            }
        }
        return Math.log10( sumPVarInGaussian / ((double) numRandomDraws) );
    }

    protected List<MultivariateGaussian> getModelGaussians() {return Collections.unmodifiableList(gaussians);}

    protected int getNumAnnotations() {return empiricalMu.length;}
}