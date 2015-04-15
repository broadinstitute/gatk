package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ml.distance.ChebyshevDistance;
import org.apache.commons.math3.ml.distance.DistanceMeasure;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.vqsr.GaussianMixtureModel.MultivariateGaussian;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static java.lang.Double.isFinite;
import static java.lang.Math.*;
import static org.broadinstitute.hellbender.utils.MathUtils.identityMatrix;
import static org.testng.Assert.*;

public class GaussianMixtureModelUnitTest extends BaseTest{

    @BeforeMethod
    private void resetRandom(){
        Utils.resetRandomGenerator();
    }

    @Test
    public void test1(){
        List<VariantDatum> data = readData(new File(CommandLineProgramTest.getTestDataDir() + "/GaussianMixtureModelUnitTest/" + "gmmTest1.goodData.10k.csv.gz"));

        VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.SHRINKAGE = 1.0;
        VRAC.MAX_GAUSSIANS = 2;
        VRAC.DIRICHLET_PARAMETER = 0.001;
        VRAC.PRIOR_COUNTS = 19.0;
        final GaussianMixtureModel gmm = GaussianMixtureModel.makeFittedModel(data, VRAC.MAX_GAUSSIANS, VRAC);
        assertTrue(!gmm.failedToConverge());
        assertTrue(!gmm.readyForEvaluation());
    }

    private List<VariantDatum> readData(File f) {
        List<VariantDatum> vd = new ArrayList<>();
        for ( String line : new XReadLines(f, true) ) {
            String[] parts = line.split(",");
            VariantDatum datum = new VariantDatum();
            datum.annotations = new ArrayRealVector(parts.length);
            for(int i = 0; i < parts.length; i++){
                datum.annotations.setEntry(i, Double.valueOf(parts[i]));
            }
            vd.add(datum);
        }

        return vd;
    }

    private static class TestConfig{
        public final double[][] mus;
        public final double[][][] sigmas;
        public final double[] weigths;
        private final int ndataPoints;
        public final double zeroZeroEval; //value of the likelihood at the 0,0 point
        public final double[] zeroZeroEvalPerDim;
        public final double[] zeroZeroEvalMarginalized;

        TestConfig(double[][] mus, double[][][] sigmas, double[] weigths, int ndataPoints, double zeroZeroEval, double[] zeroZeroEvalPerDim, double[] zeroZeroEvalMarginalized){
            this.mus = mus;
            this.sigmas = sigmas;
            this.weigths = weigths;
            this.ndataPoints = ndataPoints;

            //These are the expected values of the likelihood for the zero point.
            this.zeroZeroEval = zeroZeroEval;
            this.zeroZeroEvalPerDim = zeroZeroEvalPerDim;
            this.zeroZeroEvalMarginalized = zeroZeroEvalMarginalized;
        }

        @Override
        public String toString() {
            return "mus:" + Arrays.deepToString(mus) + " \n"
                    + "sigmas:" + Arrays.deepToString(sigmas) + "\n"
                    +  "weigths:" + Arrays.toString(weigths) + "\n"
                    + "ndatapoints:" + ndataPoints + "\n"
                    + "zeroZeroEval:" + zeroZeroEval + "\n"
                    + "zeroZeroEvalPerDim:" + Arrays.toString(zeroZeroEvalPerDim) + "\n"
                    + "zeroZeroEvalMarginalized:" + Arrays.toString(zeroZeroEvalMarginalized);
        }
    }

    @DataProvider(name="gaussians")
    public Object[][] gaussians() {
        List<Object[]> tests = new ArrayList<>();
        //Note: using local scopes to reuse variable names

        {
            final double[][] mus = {{1, 1}, {-1, -1}};
            final double[][][] sigmas = {{{1, 0}, {0, 1}}, {{0.5, 0}, {0, 2}}};
            final double[] weigths = {0.3, 0.7};
            final int ndataPoints = 100000;
            final double zeroZeroEval = -1.3055351341487855;
            final double[] zeroZeroEvalPerDim = {-0.8292145993103814, -0.7081560796687412};
            final double[] zeroZeroEvalMarginalized = {-1.2863457273966425, -1.162161863302764};
            tests.add(new Object[]{new TestConfig(mus, sigmas, weigths, ndataPoints, zeroZeroEval, zeroZeroEvalPerDim, zeroZeroEvalMarginalized)});
        }

        {
            final double[][] mus = {{1, 1}, {-1, -1}};
            final double[][][] sigmas = {{{1, 0.3}, {0.3, 0.6}}, {{1, 0}, {0, 2}}};
            final double[] weigths = {0.2, 0.8};
            final int ndataPoints = 100000;
            final double zeroZeroEval = -1.227312735746487;
            final double[] zeroZeroEvalPerDim = {-0.6162371751306832, -0.7594877507655053};
            final double[] zeroZeroEvalMarginalized = {-1.2999591422769814, -1.2576454723687285};
            tests.add(new Object[]{new TestConfig(mus, sigmas, weigths, ndataPoints, zeroZeroEval, zeroZeroEvalPerDim, zeroZeroEvalMarginalized)});
        }

        {
            //try 2 gaussians in 3 dimensions
            final double[][] mus = {{1, 1, -1}, {-1, -1, 1}};
            final double[][][] sigmas = {{{1, 0.3, 0.5}, {0.3, 0.6, 0.1}, {0.5, 0.1, 1.5}}, {{1, 0, 0}, {0, 2, 0}, {0, 0, 4}}};
            final double[] weigths = {0.2, 0.8};
            final int ndataPoints = 100000;
            final double zeroZeroEval = -1.9951548539045016;
            final double[] zeroZeroEvalPerDim = {-0.6162371751306832, -0.7594877507655053, -0.9210841500010781};
            final double[] zeroZeroEvalMarginalized = {-2.0364676215685167, -1.9772127396065897, -1.9209800981138132};
            tests.add(new Object[]{new TestConfig(mus, sigmas, weigths, ndataPoints, zeroZeroEval, zeroZeroEvalPerDim, zeroZeroEvalMarginalized)});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "gaussians")
    public void testDataEvaluationUnderModel(TestConfig tc){

        double[][] mus = tc.mus;
        double[][][] sigmas = tc.sigmas;
        double[] weights = tc.weigths;
        final int ngaussians = mus.length;
        final int ndimensions = mus[0].length;

        VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.SHRINKAGE = 1.0;
        VRAC.MAX_GAUSSIANS = 2;
        VRAC.DIRICHLET_PARAMETER = 0.001;
        VRAC.PRIOR_COUNTS = 20.0;
        final GaussianMixtureModel gmm = GaussianMixtureModel.makeEmptyModel(VRAC.MAX_GAUSSIANS, mus[0].length, VRAC.SHRINKAGE, VRAC.DIRICHLET_PARAMETER, VRAC.PRIOR_COUNTS);
        List<MultivariateGaussian> mvns = new ArrayList<>();
        for(int i = 0; i< ngaussians; i++){
            double pMixtureLog10  = log10(weights[i]);
            double[] mu = mus[i];
            double[][] sigmaData = sigmas[i];

            final RealVector prior_m = new ArrayRealVector(ndimensions);
            final RealMatrix prior_L = identityMatrix(ndimensions).scalarMultiply(200.0);

            MultivariateGaussian g = new MultivariateGaussian(ndimensions, VRAC.DIRICHLET_PARAMETER, VRAC.SHRINKAGE, VRAC.PRIOR_COUNTS, prior_m, prior_L);
            g.setpMixtureLog10(pMixtureLog10);
            g.setMu(new ArrayRealVector(mu));
            g.setParam_S(new Array2DRowRealMatrix(sigmaData));
            mvns.add(g);
        }

        final double prior = 17.0;
        gmm.setGaussians(mvns);
        VariantDatum vdZero= new VariantDatum();
        vdZero.annotations = new ArrayRealVector(ndimensions);
        vdZero.isNull = new boolean[ndimensions];//init to false
        vdZero.prior = prior;

        gmm.setLodFromModel(Arrays.asList(vdZero), true);
        assertEquals(vdZero.lod, tc.zeroZeroEval, 1.0 / 100000);
        gmm.setLodFromModel(Arrays.asList(vdZero), false);
        assertEquals(vdZero.lod, prior,1.0 / 100000);  //after this it should be prior + positive - negative = prior (because here pos = neg).

        assertEquals(gmm.evaluateDatum(vdZero), tc.zeroZeroEval, 1.0 / 100000);
        for(int dim = 0; dim < ndimensions; dim++){
            vdZero.isNull = new boolean[ndimensions];  //init to false
            assertEquals(gmm.evaluateDatumInOneDimension(vdZero, dim), tc.zeroZeroEvalPerDim[dim], 1.0 / 100000, "dim:" + dim + " diff:" + abs(gmm.evaluateDatumInOneDimension(vdZero, dim) - tc.zeroZeroEvalPerDim[dim]));
            vdZero.isNull[dim] = true;
            assertTrue(abs((gmm.evaluateDatum(vdZero) - tc.zeroZeroEvalMarginalized[dim]) / tc.zeroZeroEvalMarginalized[dim]) < 0.1, " marginalized dim:" + dim + " actual:" + gmm.evaluateDatum(vdZero) + " expected:" + tc.zeroZeroEvalMarginalized[dim]);
            vdZero.isNull[dim] = false;
        }

        //testing the infinity case for positive model - lod should be finite again after negative model
        vdZero.lod = Double.NEGATIVE_INFINITY;
        gmm.setLodFromModel(Arrays.asList(vdZero), false);
        assertTrue(isFinite(vdZero.lod));

        vdZero.setWorstPerformingAnnotation(gmm, gmm);
        assertNotEquals(vdZero.worstAnnotation, -1, "worst annotation");//using same model for pos and neg so it's arbitrary which annotation is selected

        Assert.assertNotNull(gmm.toString());//just checking that it does not blow up
    }

    @Test(dataProvider = "gaussians")
    public void testFitModelToSimulatedData(TestConfig tc){
        double[][] mus = tc.mus;
        double[][][] sigmas = tc.sigmas;
        double[] weights = tc.weigths;
        int ndata = tc.ndataPoints;

        final List<double[]> samples = sampleFromMixtureOfTwoGaussians(ndata, mus, sigmas, weights);
        List<VariantDatum> data = samples.stream().map(GaussianMixtureModelUnitTest::makeDatum).collect(Collectors.toList());

        final GaussianMixtureModel gmm = fitModel(data);

        double okWeightsDistance = 0.1;
        double okMeansDistance = 0.1;
        double okCovarianceDistance = 0.1;

        double[] okDistances = {okWeightsDistance, okMeansDistance, okCovarianceDistance};
        compareFitModelToExpected(mus, sigmas, weights, gmm, okDistances);
    }

    private GaussianMixtureModel fitModel(List<VariantDatum> data) {
        VariantRecalibratorArgumentCollection VRAC = new VariantRecalibratorArgumentCollection();
        VRAC.SHRINKAGE = 1.0;
        VRAC.MAX_GAUSSIANS = 2;
        VRAC.DIRICHLET_PARAMETER = 0.001;
        VRAC.PRIOR_COUNTS = 20.0;

        return GaussianMixtureModel.makeFittedModel(data, VRAC.MAX_GAUSSIANS, VRAC);
    }

    //XXX okDistances are in order: weights, means, covariance matrices
    private void compareFitModelToExpected(double[][] mus, double[][][] sigmas, double[] weights, GaussianMixtureModel gmm, double[] okDistances) {
        //Note: sort inferred gaussians because the model is not fully identifiable (up to reordering).
        //we sort them by the mixture weights
        final List<MultivariateGaussian> gaussians = new ArrayList<>(gmm.getGaussians());
        gaussians.sort((g1, g2) -> Double.compare(g1.getpMixtureLog10(), g2.getpMixtureLog10()));

        //now compare inferred results to truth
        MultivariateGaussian result1 = gaussians.get(0);
        MultivariateGaussian result2 = gaussians.get(1);
        assertEquals(gaussians.size(), 2);

        compareWeights(weights, result1, result2, okDistances[0]);
        compareMeans(mus, result1, result2, okDistances[1]);
        compareCovariances(sigmas, result1, result2, okDistances[2]);
    }

    private void compareWeights(double[] weights, MultivariateGaussian result1, MultivariateGaussian result2, double okDistance) {
        final double[] resultWeigths = {pow(10, result1.getpMixtureLog10()), pow(10, result2.getpMixtureLog10())};
        for (int i = 0 ; i < resultWeigths.length; i++) {
            assertTrue(distance(resultWeigths[i], weights[i]) < okDistance, "expected:" + weights[i] + " actual:" + resultWeigths[i]) ;
        }
    }

    private void compareMeans(double[][] mus, MultivariateGaussian result1, MultivariateGaussian result2, double okDistance) {
        final double[][] resultMeans = {result1.getXBar().toArray(), result2.getXBar().toArray()};
        for (int i = 0 ; i < resultMeans.length; i++) {
            double dist = distance(resultMeans[i], mus[i]);
            assertTrue(dist < okDistance, "expected:" + Arrays.toString(mus[i]) + " actual:" + Arrays.toString(resultMeans[i]) + " distance:" + dist);
        }
    }

    private void compareCovariances(double[][][] sigmas, MultivariateGaussian result1, MultivariateGaussian result2, double okDistance) {
        final double[][][] resultCovars = {result1.getParam_S().getData(), result2.getParam_S().getData()};
        for (int i = 0 ; i < resultCovars.length; i++) {
            assertTrue(distance(resultCovars[i], sigmas[i]) < okDistance, "expected:" + Arrays.deepToString(sigmas[i]) + " actual:" + Arrays.deepToString(resultCovars[i])) ;
        }
    }

    //XXX it is bit lame to compare matrices this way. This returns a max of the abs difference, point-wise.
    private static double distance(double[][] a, double[][] b) {
        RealMatrix m1 = new Array2DRowRealMatrix(a);
        RealMatrix m2 = new Array2DRowRealMatrix(b);
        RealMatrix m1MinusM2 = m1.subtract(m2);
        final double[][] diff = m1MinusM2.getData();
        double max = Double.MIN_VALUE;
        for(int i = 0; i < diff.length; i++){
            for(int j = 0; j < diff[i].length; j++){
                if (abs(diff[i][j]) < max){
                    max = diff[i][j];
                }
            }
        }
        return max;
    }

    private static double distance(double[] a, double[] b) {
        DistanceMeasure dist = new ChebyshevDistance();    //max of abs diffs
        return dist.compute(a, b);
    }

    private static double distance(double a, double b) {
        return abs(a - b);
    }

    public static VariantDatum makeDatum(double[] ds){
        VariantDatum vd = new VariantDatum();
        vd.annotations = new ArrayRealVector(ds);
        vd.isNull = new boolean[ds.length];//initialized to false
        return vd;
    }

    public List<double[]> sampleFromMixtureOfTwoGaussians(int nSamples, double[][] mus, double[][][] sigmas, double[] weigths){
        int nGaussians = mus.length;

        //Note: This is the apache math Pair, not the apache lang Pair.
        List<Pair<MultivariateNormalDistribution, Double>> pmf = new ArrayList<>(nGaussians);
        for(int i = 0; i < nGaussians; i++){
            MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(Utils.getApacheRandomGenerator(), mus[i], sigmas[i]);
            pmf.add(new Pair<>(mvn, weigths[i]));
        }

        EnumeratedDistribution<MultivariateNormalDistribution> catDistro = new EnumeratedDistribution<>(Utils.getApacheRandomGenerator(), pmf);
        List<double[]> result = new ArrayList<>(nSamples);
        for (int i = 0; i < nSamples; i++) {
            MultivariateNormalDistribution mvn = catDistro.sample();
            result.add(mvn.sample());
        }
        return result;
    }
}
