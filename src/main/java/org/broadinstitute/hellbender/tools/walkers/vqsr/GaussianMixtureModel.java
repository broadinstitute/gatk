package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.annotations.VisibleForTesting;
import joptsimple.internal.Strings;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.ExpandingArrayList;

import java.util.*;

import static java.lang.Double.isInfinite;
import static java.lang.Math.log;
import static java.lang.Math.log10;
import static java.util.Comparator.comparing;
import static java.util.stream.Collectors.minBy;
import static org.apache.commons.math3.special.Gamma.digamma;
import static org.broadinstitute.hellbender.utils.MathUtils.*;
import static org.broadinstitute.hellbender.utils.Utils.getRandomGenerator;

/**
 * Defines functionality for fitting gaussian mixture models to data and computing likelihood of data points.
 */
final class GaussianMixtureModel {

    public final static double MIN_ACCEPTABLE_LOD_SCORE = -20000.0;
    public final static double MIN_PROB_CONVERGENCE = 2E-3;

    protected final static Logger logger = LogManager.getLogger(GaussianMixtureModel.class);

    //The gaussians we're fitting. Order is undefined.
    private final Collection<MultivariateGaussian> gaussians;

    private final double hyperPriorParameter_beta;
    private final double hyperPriorParameter_alpha;  //parameter of the symmetric dirichlet prior on the mixture weigths
    private final double hyperPriorParameter_nu;

    private final double[] empiricalMu;      //Note: stays fixed as 0 vector
    private final RealMatrix empiricalSigma; //Note: stays fixed as diagonal matrix of 1/200.

    private boolean isModelReadyForEvaluation;    //set to true after the EM iterations converge and the precomputeDenominatorForEvaluation method is called
    private boolean failedToConverge;             //set to true if something went wrong in the fitting

    private GaussianMixtureModel(final int numGaussians, final int numAnnotations,
                                                  final double hyperPriorParameter_beta, final double hyperPriorParameter_alpha, final double hyperPriorParameter_nu) {
        gaussians = new HashSet<>( numGaussians );
        for( int i = 0; i < numGaussians; i++ ) {
            gaussians.add( new MultivariateGaussian( numAnnotations ) );
        }
        this.hyperPriorParameter_beta = hyperPriorParameter_beta;
        this.hyperPriorParameter_alpha = hyperPriorParameter_alpha;
        this.hyperPriorParameter_nu = hyperPriorParameter_nu;
        this.empiricalMu = new double[numAnnotations];
        this.empiricalSigma = MultivariateGaussian.inverse(identityMatrix(numAnnotations).scalarMultiply(200.0));

        this.isModelReadyForEvaluation = false;
        this.failedToConverge = false;
    }

    /**
     * Makes an empty gaussian mixture model.
     * This is useful mostly in testing, where we want to set the gaussians rather than infer them.
     */
    public static GaussianMixtureModel makeEmptyModel(final int numGaussians, final int numAnnotations,
                                final double shrinkage, final double dirichletParameter, final double priorCounts) {
        return new GaussianMixtureModel(numGaussians, numAnnotations, shrinkage, dirichletParameter, priorCounts);
    }

    /**
     * Fits a model to the provided data and returns a fitted model.
     */
    public static GaussianMixtureModel makeFittedModel( final List<VariantDatum> data, final int maxGaussians, final VariantRecalibratorArgumentCollection VRAC) {
        if( data == null || data.isEmpty() ) { throw new IllegalArgumentException("No data found."); }
        if( maxGaussians <= 0 ) { throw new IllegalArgumentException("maxGaussians must be a positive integer but found: " + maxGaussians); }

        final GaussianMixtureModel model = GaussianMixtureModel.makeEmptyModel(maxGaussians, data.get(0).annotations.length, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER, VRAC.PRIOR_COUNTS );
        model.initializeRandomModel(data, VRAC.NUM_KMEANS_ITERATIONS, Utils.getRandomGenerator());
        model.variationalBayesExpectationMaximization(data, VRAC.MAX_ITERATIONS, MIN_PROB_CONVERGENCE);
        return model;
    }

    /**
     * Sets the gaussians to the arguments.
     * Note that the actual data is deep copied, for safety. For testing only.
     */
    @VisibleForTesting
    void setGaussians(Collection<MultivariateGaussian> mvgs) {
        gaussians.clear();
        for(MultivariateGaussian mvg : mvgs ){
            gaussians.add(mvg.deepCopy());
        }
    }

    /**
     * Returns the (deep copies of the) gaussians. Mostly for testing but not limited to it.
     */
    @VisibleForTesting
    List<MultivariateGaussian> getGaussians() {
        List<MultivariateGaussian> result = new ArrayList<>(gaussians.size());
        for(MultivariateGaussian mvg : gaussians ){
            result.add(mvg.deepCopy());
        }
        return result;
    }

    /**
     * Fits the Gaussian Mixture model to the data.
     * The model is assumed to have been initialized (eg by k-means).
     */
    public void variationalBayesExpectationMaximization( final List<VariantDatum> data, int maxIterations, double minProbConvergence ) {

        // The VBEM loop
        normalizePMixtureLog10(data.size());
        expectationStep(data);
        double currentChangeInMixtureCoefficients;
        int iteration = 0;
        logger.info("Finished iteration " + iteration + ".");
        while( iteration < maxIterations ) {
            iteration++;
            maximizationStep(data);
            currentChangeInMixtureCoefficients = normalizePMixtureLog10(data.size());
            expectationStep(data);
            if( iteration % 5 == 0 ) { // cut down on the number of output lines so that users can read the warning messages
                logger.info("Finished iteration " + iteration + ". \tCurrent change in mixture coefficients = " + String.format("%.5f", currentChangeInMixtureCoefficients));
            }
            if( iteration > 2 && currentChangeInMixtureCoefficients < minProbConvergence ) {
                logger.info("Convergence after " + iteration + " iterations!");
                break;
            }
        }

        evaluateFinalModelParameters(data);
    }

    /**
     * Initializes the model by running k-means on the data. K is the same as the number of gaussians.
     */
    public void initializeRandomModel( final List<VariantDatum> data, final int numKMeansIterations, Random rand) {

        gaussians.forEach(g -> g.initializeRandomMu(rand));

        logger.info( "Initializing model with " + numKMeansIterations + " k-means iterations..." );
        initializeMeansUsingKMeans(data, numKMeansIterations, rand);

        // initialize uniform mixture coefficients, random covariance matrices, and initial hyperparameters
        for( final MultivariateGaussian gaussian : gaussians ) {
            gaussian.setpMixtureLog10(log10(1.0 / ((double) gaussians.size())));
            gaussian.setN_k((1.0 * data.size()) / ((double) gaussians.size()));
            gaussian.initializeRandomSigma(rand);
            gaussian.setHyperParameter_nu(hyperPriorParameter_nu);
            gaussian.setHyperParameter_beta(hyperPriorParameter_beta);
            gaussian.setHyperParameter_alpha(hyperPriorParameter_alpha);
        }
    }

    //performs k-means clustering for a fixed number of iterations to initialize the gaussian means.
    private void initializeMeansUsingKMeans( final List<VariantDatum> data, final int numIterations, final Random rand ) {

        final Map<VariantDatum, MultivariateGaussian> assignment = new HashMap<>(data.size());

        int iter = 0;
        while( iter++ < numIterations ) {
            // E step: assign each variant to the nearest cluster
            data.forEach(d -> assignment.put(d, gaussians.stream().collect(minBy(comparing(mvg -> mvg.calculateDistanceFromMeanSquared(d)))).get()));

            // M step: update gaussian means based on assigned variants
            //NOTE: this could be sped up by traversing data only once
            gaussians.forEach(g -> g.zeroOutMu());
            final Map<MultivariateGaussian, Integer> numAssigned = new HashMap<>();
            gaussians.forEach(g -> numAssigned.put(g, 0));

            for( Map.Entry<VariantDatum, MultivariateGaussian> kv : assignment.entrySet() ) {
                MultivariateGaussian mvg = kv.getValue();
                mvg.incrementMu(kv.getKey());
                numAssigned.put(mvg, 1 + numAssigned.get(mvg));
            }
            gaussians.forEach(g -> {
                final int n = numAssigned.get(g);
                if( n != 0 ) {
                    g.divideEqualsMu((double) n);
                } else {
                    g.initializeRandomMu(rand);//no data item assigned to this gaussian - reset the mean randomly
                }
            });
        }
    }

    private void expectationStep( final List<VariantDatum> data ) {

        gaussians.forEach(g -> g.precomputeDenominatorForVariationalBayes(getSumHyperParameterAlpha()));

        for( final VariantDatum datum : data ) {
            final double[] pVarInGaussianLog10 = new double[gaussians.size()];
            int gaussianIndex = 0;
            for( final MultivariateGaussian gaussian : gaussians ) {
                pVarInGaussianLog10[gaussianIndex++] = gaussian.evaluateDatumLog10(datum);
            }
            final double[] pVarInGaussianNormalized = normalizeFromLog10(pVarInGaussianLog10);
            gaussianIndex = 0;
            for( final MultivariateGaussian gaussian : gaussians ) {
                gaussian.assignPVarInGaussian( pVarInGaussianNormalized[gaussianIndex++] );
            }
        }
    }

    private void maximizationStep( final List<VariantDatum> data ) {
        gaussians.forEach(g -> g.maximizeGaussian( data, empiricalMu, empiricalSigma, hyperPriorParameter_beta, hyperPriorParameter_alpha, hyperPriorParameter_nu));
    }

    private double getSumHyperParameterAlpha() {
        return gaussians.stream().mapToDouble(MultivariateGaussian::getHyperParameter_alpha).sum();
    }

    private void evaluateFinalModelParameters( final List<VariantDatum> data ) {
        gaussians.forEach(g -> g.evaluateFinalModelParameters(data));
        normalizePMixtureLog10(data.size());
    }

    private double normalizePMixtureLog10(int dataSize) {
        final double sumPK = gaussians.stream().mapToDouble(MultivariateGaussian::getN_k).sum();  //Note: sumPK should be == data size here
        if (Math.abs(sumPK - (double)dataSize) > 0.1) {
            throw new IllegalStateException("wrong sum of responsibilities " + sumPK + " dataSize:" + dataSize);
        }
        final double log10SumPK = log10(sumPK);

        int gaussianIndex = 0;
        final double[] pGaussianLog10 = new double[gaussians.size()];
        for( final MultivariateGaussian gaussian : gaussians ) {
            pGaussianLog10[gaussianIndex++] = log10(gaussian.getN_k()) - log10SumPK;
        }
        final double[] pGaussianLog10Normalized = normalizeFromLog10(pGaussianLog10, true);

        double sumDiff = 0.0;
        gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            sumDiff += Math.abs(pGaussianLog10Normalized[gaussianIndex] - gaussian.getpMixtureLog10());
            gaussian.setpMixtureLog10( pGaussianLog10Normalized[gaussianIndex++] );
        }
        return sumDiff;
    }

    private static RealMatrix identityMatrix(int n){
        double[] ones = new double[n];
        Arrays.fill(ones, 1);
        return new DiagonalMatrix(ones);
    }

    /**
     * Before data is evaluated, this methods needs to be called to precompute some numbers that are used for all values.
     */
    public void precomputeDenominatorForEvaluation() {
        gaussians.forEach(MultivariateGaussian::precomputeDenominatorForEvaluation);

        isModelReadyForEvaluation = true;
    }

    /**
     *  Return the likelihood of the data under the current model.
     *  If any dimension has isNull set to true, the missing dimensions will be marginalized
     *  over (by using random sampling).
     */
    public double evaluateDatum( final VariantDatum datum ) {
        if( !readyForEvaluation() ) {
            precomputeDenominatorForEvaluation();
        }

        for( final boolean isNull : datum.isNull ) {
            if( isNull ) {
                return evaluateDatumMarginalized( datum );
            }
        }
        // Fill an array with the log10 probability coming from each Gaussian and then use MathUtils to sum them up correctly
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            pVarInGaussianLog10[gaussianIndex++] = gaussian.getpMixtureLog10() + gaussian.evaluateDatumLog10( datum );
        }
        final double result = nanTolerantLog10SumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
        if (Double.isNaN(result)){
            failedToConverge = true;
        }
        return result;
    }

    /**
     * Return the likelihood of the data under the current model but using only 1 dimension of the gaussians.
     * Used only to decide which covariate dimension is most divergent in order to report in the culprit info field annotation
     */
    public Double evaluateDatumInOneDimension( final VariantDatum datum, final int iii ) {
        if( !readyForEvaluation() ) {
            precomputeDenominatorForEvaluation();
        }

        if(datum.isNull[iii]) {
            return null;
        }

        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        int gaussianIndex = 0;
        for( final MultivariateGaussian gaussian : gaussians ) {
            pVarInGaussianLog10[gaussianIndex++] = gaussian.getpMixtureLog10() + MathUtils.normalDistributionLog10(gaussian.getMu(iii), gaussian.getSigma(iii, iii), datum.annotations[iii]);
        }
        return nanTolerantLog10SumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
    }

    private double evaluateDatumMarginalized( final VariantDatum datum ) {
        int numRandomDraws = 0;
        double sumPVarInGaussian = 0.0;
        final int numIterPerMissingAnnotation = 1000; // Trade off here between speed of computation and accuracy of the marginalization
        final double[] pVarInGaussianLog10 = new double[gaussians.size()];
        // for each dimension
        for( int i = 0; i < datum.annotations.length; i++ ) {
            // if it is missing then marginalize over the missing dimension by drawing X random values for the missing annotation and averaging the lod
            if( datum.isNull[i] ) {
                for( int j = 0; j < numIterPerMissingAnnotation; j++ ) {
                    datum.annotations[i] = Utils.getRandomGenerator().nextGaussian(); // draw a random sample from the standard normal distribution

                    // evaluate this random data point
                    int gaussianIndex = 0;
                    for( final MultivariateGaussian gaussian : gaussians ) {
                        pVarInGaussianLog10[gaussianIndex++] = gaussian.getpMixtureLog10() + gaussian.evaluateDatumLog10( datum );
                    }

                    // add this sample's probability to the pile in order to take an average in the end
                    sumPVarInGaussian += Math.pow(10.0, nanTolerantLog10SumLog10(pVarInGaussianLog10)); // p = 10 ^ Sum(pi_k * p(v|n,k))
                    numRandomDraws++;
                }
            }
        }
        if (numRandomDraws == 0){
            failedToConverge = true;
        }
        return log10(sumPVarInGaussian / ((double) numRandomDraws));
    }

    @Override
    public String toString() {
        List<String> sb = new ArrayList<>();
        sb.add("hyperPriorParameter_beta:" + hyperPriorParameter_beta);
        sb.add("hyperPriorParameter_alpha:" + hyperPriorParameter_alpha);
        sb.add("hyperPriorParameter_nu:" + hyperPriorParameter_nu);
        sb.add("empiricalMu:" + Arrays.toString(empiricalMu));
        sb.add("empiricalSigma:" + Arrays.deepToString(empiricalSigma.getData()));

        sb.add("numGaussians:" +  gaussians.size());
        for (MultivariateGaussian mvg: gaussians){
            sb.add("------------");
            sb.add(mvg.toString());
        }
        return Utils.join("\n", sb);
    }

    /**
     * Returns whether the model has been fit and it ready for evaluation.
     * You need to call precomputeDenominatorForEvaluation to make it so.
     */
    public boolean readyForEvaluation() {
        return isModelReadyForEvaluation;
    }

    /**
     * Returns whether the model fitting failed to converge.
     */
    public boolean failedToConverge() {
        return failedToConverge;
    }

    /**
     * Computes the lod score for each data item and assigns it to the item.
     */
    public void setLodFromModel( final List<VariantDatum> data, final boolean isPositiveModel) {
        if( !readyForEvaluation() ) {
            precomputeDenominatorForEvaluation();
        }
        if (failedToConverge()){
            return;
        }

        logger.info("Evaluating full set of " + data.size() + " variants on positive model...");
        for( final VariantDatum datum : data ) {
            datum.lod = computeLod(datum, isPositiveModel);
            if( failedToConverge() ) {
                return;
            }
        }
    }

    private double computeLod(VariantDatum datum, boolean isPositiveModel){
        final double thisLod = evaluateDatum(datum);
        if (isPositiveModel){
            return thisLod; // positive model only so set the lod and return
        }

        final double positiveLod = datum.lod;
        final double negativeLod = thisLod;
        if (isInfinite(positiveLod)) {
            return (MIN_ACCEPTABLE_LOD_SCORE + getRandomGenerator().nextDouble() * MIN_ACCEPTABLE_LOD_SCORE); // Negative infinity lod values are possible when covariates are extremely far away from their tight Gaussians
        }

        return datum.prior + positiveLod - negativeLod;
    }

    @VisibleForTesting
    /**
     * This class represents one gaussian in the mixture model.
     */
    static final class MultivariateGaussian {
        private double pMixtureLog10;  //log10 of the mixture weight for this gaussian
        private double N_k;        //sum of the fractional weigths assigned to this gaussian from all the data points
        private double[] mu;           //mean vector
        private RealMatrix sigma;      //covariance matrix
        private double hyperParameter_nu;
        private double hyperParameter_beta;
        private double hyperParameter_alpha;
        private final ExpandingArrayList<Double> pVarInGaussian;

        //caches
        private double cachedDenomLog10;
        private RealMatrix cachedSigmaInverse;

        MultivariateGaussian(double pMixtureLog10, double N_k, double[] mu, double[][] sigmaData, double hyperParameter_nu, double hyperParameter_beta, double hyperParameter_alpha) {
            this.pMixtureLog10 = pMixtureLog10;
            this.N_k = N_k;
            this.mu = Arrays.copyOf(mu, mu.length);
            this.sigma = new Array2DRowRealMatrix(sigmaData);   //makes a copy of the array
            this.hyperParameter_nu = hyperParameter_nu;
            this.hyperParameter_beta = hyperParameter_beta;
            this.hyperParameter_alpha = hyperParameter_alpha;
            this.pVarInGaussian = new ExpandingArrayList<>();
            //TODO: add checks for the parameter values
        }

        MultivariateGaussian(final int numAnnotations) {
            this(0, 0, new double[numAnnotations], new double[numAnnotations][numAnnotations], 0, 0, 0);
        }

        MultivariateGaussian deepCopy() {
            return new MultivariateGaussian(
                    pMixtureLog10,
                    N_k,
                    mu,              //data copies by the constructor so pass the array here (no copy)
                    sigma.getData(), //data copies by the constructor so pass the array here (no copy)
                    hyperParameter_nu,
                    hyperParameter_beta,
                    hyperParameter_alpha
            );
        }

        @Override
        public String toString() {
            List<String> sb = new ArrayList<>();
            sb.add("pMixtureLog10:" + pMixtureLog10);
            sb.add("N_k:" + N_k);
            sb.add("mu:" + Arrays.toString(mu));
            sb.add("sigma:" + Arrays.deepToString(sigma.getData()));
            sb.add("hyperParameter_nu:" + hyperParameter_nu);
            sb.add("hyperParameter_beta:" + hyperParameter_beta);
            sb.add("hyperParameter_alpha:" + hyperParameter_alpha);
            sb.add("cachedDenomLog10:" + cachedDenomLog10);
            sb.add("pVarInGaussian:" + pVarInGaussian);

            return Strings.join(sb, "\n");
        }
        void zeroOutMu() {
            Arrays.fill(mu, 0.0);
        }

        void zeroOutSigma() {
            sigma = new Array2DRowRealMatrix(new double[getNumDimensions()][getNumDimensions()]);
        }

        public static final double MU_MIN = -4.0;
        public static final double MU_MAX = 4.0;
        public static final double MU_SPAN = MU_MAX - MU_MIN;

        void initializeRandomMu(final Random rand) {
            for (int i = 0; i < getNumDimensions(); i++) {
                mu[i] = MU_MIN + MU_SPAN * rand.nextDouble();
            }
        }

        void initializeRandomSigma(final Random rand) {
            //This is equiv to drawing from Wishart.
            //Note: maybe we want empirical cov from kmeans clusters
            final double[][] randSigma = new double[getNumDimensions()][getNumDimensions()];
            for (int i = 0; i < getNumDimensions(); i++) {
                for (int j = i; j < getNumDimensions(); j++) {
                    randSigma[j][i] = 0.55 + 1.25 * rand.nextDouble();
                    if (rand.nextBoolean()) {
                        randSigma[j][i] *= -1.0;
                    }
                }
            }
            // Sigma is a symmetric, positive-definite matrix created by taking a lower triangular
            // matrix and multiplying it by its transpose

            RealMatrix tmp = new Array2DRowRealMatrix(randSigma);
            sigma = tmp.multiply(tmp.transpose());
        }

        double calculateDistanceFromMeanSquared(final VariantDatum datum) {
            return distanceSquared(datum.annotations, mu);
        }

        void incrementMu(final VariantDatum datum) {
            incrementMu(datum, 1.0);
        }

        void incrementMu(final VariantDatum datum, final double prob) {
            for (int i = 0; i < getNumDimensions(); i++) {
                mu[i] += prob * datum.annotations[i];
            }
        }

        void divideEqualsMu(final double x) {
            for (int i = 0; i < getNumDimensions(); i++) {
                mu[i] /= x;
            }
        }

        private void precomputeInverse() {
            try {
                cachedSigmaInverse = inverse(sigma);
            } catch (MathIllegalArgumentException e) {
                throw new UserException("Error during clustering. Most likely there are too few variants used during Gaussian mixture modeling. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example).");
            }
        }

        static RealMatrix inverse(RealMatrix m) {
            return new LUDecomposition(m).getSolver().getInverse();
        }

        private static double determinant(RealMatrix m) {
            return new LUDecomposition(m).getDeterminant();
        }

        void precomputeDenominatorForEvaluation() {
            precomputeInverse();
            cachedDenomLog10 = log10(Math.pow(2.0 * Math.PI, -1.0 * ((double) getNumDimensions()) / 2.0)) + log10(Math.pow(determinant(sigma), -0.5));
        }

        void precomputeDenominatorForVariationalBayes(final double sumHyperParameterAlpha) {
            precomputeInverse();

            //Note: from now on,  cachedSigmaInverse actually means cachedSigmaInverse * nu
            cachedSigmaInverse = cachedSigmaInverse.scalarMultiply(hyperParameter_nu);

            //Murphy eq. 21.129
            final double logOfPiTilde = digamma(hyperParameter_alpha) - digamma(sumHyperParameterAlpha);

            //Murphy eq. 21.131
            double logOfLambdaTilde = 0.0;
            for (int j = 1; j <= getNumDimensions(); j++) {    //note: using j to conform with Murphy
                logOfLambdaTilde += digamma((hyperParameter_nu + 1.0 - j) / 2.0);
            }
            logOfLambdaTilde -= log(determinant(sigma));
            logOfLambdaTilde += log(2.0) * getNumDimensions();

            //Murphy eq. 21.134 (multiplicative factor that does not depend on i -- index over data)
            final double x = (-1.0 * getNumDimensions()) / (2.0 * hyperParameter_beta);
            cachedDenomLog10 = lnToLog10(logOfPiTilde + 0.5 * logOfLambdaTilde + x);
        }

        private static double lnToLog10(final double x){
            return x / NATURAL_LOG_OF_TEN;
        }

        double evaluateDatumLog10(final VariantDatum datum) {
            double sumKernel = 0.0;
            final double[] crossProdTmp = new double[getNumDimensions()];
            for (int i = 0; i < getNumDimensions(); i++) {
                for (int j = 0; j < getNumDimensions(); j++) {
                    crossProdTmp[i] += (datum.annotations[j] - mu[j]) * cachedSigmaInverse.getEntry(j, i);
                }
            }
            for (int j = 0; j < getNumDimensions(); j++) {
                sumKernel += crossProdTmp[j] * (datum.annotations[j] - mu[j]);
            }

            return ((-0.5 * sumKernel) / NATURAL_LOG_OF_TEN) + cachedDenomLog10; // This is the definition of a Gaussian PDF Log10
        }

        void assignPVarInGaussian(final double pVar) {
            pVarInGaussian.add(pVar);
        }

        void resetPVarInGaussian() {
            pVarInGaussian.clear();
        }

        void maximizeGaussian(final List<VariantDatum> data, final double[] empiricalMu, final RealMatrix empiricalSigma,
                                     final double SHRINKAGE, final double DIRICHLET_PARAMETER, final double DEGREES_OF_FREEDOM) {
            N_k = 1E-10;
            final RealMatrix wishart = new Array2DRowRealMatrix(getNumDimensions(), getNumDimensions());
            zeroOutMu();
            zeroOutSigma();

            updateMu(data);

            final double shrinkageFactor = (SHRINKAGE * N_k) / (SHRINKAGE + N_k);
            for (int i = 0; i < getNumDimensions(); i++) {
                for (int j = 0; j < getNumDimensions(); j++) {
                    wishart.setEntry(i, j, shrinkageFactor * mu[i] * mu[j]);
//                    wishart.setEntry(i, j, shrinkageFactor * (mu[i] - empiricalMu[i]) * (mu[j] - empiricalMu[j]));
                }
            }

            updateSigma(data);

            sigma = sigma.add(empiricalSigma);
            sigma = sigma.add(wishart);

            for (int i = 0; i < getNumDimensions(); i++) {
                mu[i] = (N_k * mu[i] + SHRINKAGE * empiricalMu[i]) / (N_k + SHRINKAGE);
            }

            hyperParameter_nu = N_k + DEGREES_OF_FREEDOM;
            hyperParameter_beta = N_k + SHRINKAGE;
            hyperParameter_alpha = N_k + DIRICHLET_PARAMETER;

            resetPVarInGaussian(); // clean up some memory
        }

        void evaluateFinalModelParameters(final List<VariantDatum> data) {
            N_k = 0.0;
            zeroOutMu();
            zeroOutSigma();
            updateMu(data);
            updateSigma(data);
            sigma = sigma.scalarMultiply(1.0 / N_k);

            resetPVarInGaussian(); // clean up some memory
        }

        private void updateMu(List<VariantDatum> data) {
            int datumIndex = 0;
            for (final VariantDatum datum : data) {
                final double prob = pVarInGaussian.get(datumIndex++);
                N_k += prob;
                incrementMu(datum, prob);
            }
            divideEqualsMu(N_k);
        }

        private void updateSigma(List<VariantDatum> data) {
            int datumIndex = 0;
            final int dims = getNumDimensions();
            final RealMatrix pVarSigma = new Array2DRowRealMatrix(dims, dims);
            for (final VariantDatum datum : data) {
                final double prob = pVarInGaussian.get(datumIndex++);
                for (int i = 0; i < dims; i++) {
                    for (int j = 0; j < dims; j++) {
                        pVarSigma.setEntry(i, j, prob * (datum.annotations[i] - mu[i]) * (datum.annotations[j] - mu[j]));
                    }
                }
                sigma = sigma.add(pVarSigma);
            }
        }

        double getpMixtureLog10() {
            return pMixtureLog10;
        }

        int getNumDimensions() {
            return mu.length;
        }

        double getMu(int i) {
            return mu[i];
        }

        /**
         * Returns (a copy of) the mean vector.
         * A copy is made to avoid representation exposure.
         */
        double[] getMu() {
            return Arrays.copyOf(mu, mu.length);
        }

        /**
         * Returns (a copy of) the covariance matrix.
         * A copy is made to avoid representation exposure.
         */
        RealMatrix getSigma() {
            return sigma.copy();
        }

        double getSigma(int i, int j) {
            return sigma.getEntry(i, j);
        }

        void setpMixtureLog10(double pMixtureLog10) {
            this.pMixtureLog10 = pMixtureLog10;
        }

        double getN_k() {
            return N_k;
        }

        double getHyperParameter_alpha() {
            return hyperParameter_alpha;
        }

        void setHyperParameter_alpha(double hyperParameter_alpha) {
            this.hyperParameter_alpha = hyperParameter_alpha;
        }

        void setHyperParameter_beta(double hyperParameter_beta) {
            this.hyperParameter_beta = hyperParameter_beta;
        }

        void setHyperParameter_nu(double hyperParameter_nu) {
            this.hyperParameter_nu = hyperParameter_nu;
        }

        void setN_k(double n_k) {
            this.N_k = n_k;
        }

    }
}