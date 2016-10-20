package org.broadinstitute.hellbender.tools.exome.conversion.allelicbalancecaller;

import org.apache.commons.math3.util.Pair;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.utils.SerializationTestUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;


public class AllelicSplitCallerUnitTest extends BaseTest {

    private static final String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/conversion/allelicbalancecaller/";
    private static final File ACNV_SEG_FILE = new File(TEST_DIR, "cell_line-sim-final.seg");

    @Test
    public void testMakeCalls() {

        // This mostly just checks that the calling does not crash and does produce results.
        final CNLOHCaller cnlohCaller = new CNLOHCaller();
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final List<ACNVModeledSegment> segs = SegmentUtils.readACNVModeledSegmentFile(ACNV_SEG_FILE);
        SerializationTestUtils.roundTripInKryo(segs.get(0), ACNVModeledSegment.class, ctx.getConf());

        // Make sure the CNLOH Caller is serializable before making calls.
        SerializationTestUtils.roundTripInKryo(cnlohCaller, CNLOHCaller.class, ctx.getConf());

        final List<AllelicSplitCall> calls = cnlohCaller.makeCalls(segs, 2, ctx);
        Assert.assertNotNull(calls);
        Assert.assertTrue(calls.size() > 0);
        Assert.assertTrue(calls.stream().allMatch(c -> c.getBalancedCall() != null));
        Assert.assertTrue(calls.stream().allMatch(c -> c.getCnlohCall() != null));
        Assert.assertTrue(calls.stream().allMatch(c -> c.getAcnvSegment() != null));

        // Make sure the CNLOH Caller is serializable after making calls.
        SerializationTestUtils.roundTripInKryo(cnlohCaller, CNLOHCaller.class, ctx.getConf());
        SerializationTestUtils.roundTripInKryo(calls.get(0), AllelicSplitCall.class, ctx.getConf());
    }

    @Test(dataProvider = "mafValues")
    public void testCalculateMaf(double rho, int m, int n, double gt) {
        Assert.assertEquals(CNLOHCaller.calculateMaf(rho, m, n, HomoSapiensConstants.DEFAULT_PLOIDY), gt, 1e-4);
        Assert.assertEquals(CNLOHCaller.calculateMaf(rho, n, m, HomoSapiensConstants.DEFAULT_PLOIDY), gt, 1e-4);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCalculateMafNaN() {
        CNLOHCaller.calculateMaf(Double.NaN, 1, 1, HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCalculateMafInf() {
        CNLOHCaller.calculateMaf(Double.NEGATIVE_INFINITY, 1, 1, HomoSapiensConstants.DEFAULT_PLOIDY);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCalculateMafNegative1() {
        CNLOHCaller.calculateMaf(0.5, -5, 1, HomoSapiensConstants.DEFAULT_PLOIDY);
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCalculateMafNegative2() {
        CNLOHCaller.calculateMaf(0.5, 2, -1, HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    @DataProvider(name="mafValues")
    public Object[][] mafValues() {
        return new Object[][] {
                { 1.0, 0, 0, CNLOHCaller.MIN_L },
                { 0, 0, 5, 0.5},
                { .5, 0, 2, 0.25},
                { .5, 0, 1, 0.33333333333},
                { .5, 0, 5, 0.1429},
                { 1, 0, 5, CNLOHCaller.MIN_L},
                { .5, 0, 0, 0.5}
        };
    }

    @Test(dataProvider = "fmafValues")
    public void testBasicFMaf(final double rho, final int m, final int n, final double credibleMode,
                              final double credibleLow, final double credibleHigh, final double gt) {
        final double guess = CNLOHCaller.calculateFmaf(rho, m, n, credibleMode, credibleLow, credibleHigh, HomoSapiensConstants.DEFAULT_PLOIDY);
        final double guessCNSwitched = CNLOHCaller.calculateFmaf(rho, n, m, credibleMode, credibleLow, credibleHigh, HomoSapiensConstants.DEFAULT_PLOIDY);
        Assert.assertEquals(guess, gt, 1e-9);
        Assert.assertEquals(guessCNSwitched, gt, 1e-9);
    }
    @DataProvider(name="fmafValues")
    public Object[][] fmafValues() {
        return new Object[][] {
                // rho, m, n, mode, low, high

                // truth = Minimum value
                { 1.0, 0, 0, 0.4, 0.39, 0.41, CNLOHCaller.MIN_L },

                // Matlab: 2.1846e-20
                { .5, 0, 0, 0.25, 0.2, 0.3, 2.1846e-20 },

                // Matlab: 81.1065489973630
                { .5, 0, 0, 0.499, 0.49, 0.4999, 81.1065489973630 },

                // Matlab: 0.00134286628603481
                { .3, 0, 2, 0.333, 0.33, 0.34, 0.00134286628603481}

        };
    }


    @Test(dataProvider = "crValues")
    public void testCalculateCr(double rho, int m, int n, double lambda, double gt) {
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, m, n, lambda, HomoSapiensConstants.DEFAULT_PLOIDY), gt, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, n, m, lambda, HomoSapiensConstants.DEFAULT_PLOIDY), gt, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, m, n, lambda*2, HomoSapiensConstants.DEFAULT_PLOIDY), gt/2, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, n, m, lambda/2, HomoSapiensConstants.DEFAULT_PLOIDY), gt*2, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, m, n, lambda/2, HomoSapiensConstants.DEFAULT_PLOIDY), gt*2, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, n, m, lambda*2, HomoSapiensConstants.DEFAULT_PLOIDY), gt/2, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, m, n, lambda*4, HomoSapiensConstants.DEFAULT_PLOIDY), gt/4, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, n, m, lambda/4, HomoSapiensConstants.DEFAULT_PLOIDY), gt*4, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, m, n, lambda/4, HomoSapiensConstants.DEFAULT_PLOIDY), gt*4, 1e-6);
        Assert.assertEquals(CNLOHCaller.calculateCopyRatio(rho, n, m, lambda*4, HomoSapiensConstants.DEFAULT_PLOIDY), gt/4, 1e-6);
    }

    // calculateCopyRatio(final double rho, final int m, final int n, final double lambda), ground truth
    @DataProvider(name="crValues")
    public Object[][] crValues() {
        return new Object[][] {
                { 1.0, 0, 0, 2, 0},
                { 0, 0, 5, 2, 1}, // rho == 0 --> CR = 1
                { .5, 0, 2, 2, 1},
                { .5, 0, 1, 2, 0.75},
                { .5, 0, 5, 2, 1.75},
                { 1, 0, 5, 2, 2.5},
                { 1, 0, 0, 2, 0}, //hom del at CCF & purity of 1
                { 0.5, 0, 0, 2, 0.5}, //hom del at CCF * purity of 0.5
        };
    }

    @Test(dataProvider = "3dArray")
    public void testSumOnlyFirstDimension(double[][][] array3d, double[][] gt) {
        Assert.assertEquals(CNLOHCaller.sumOverFirstDimension(array3d), gt);
    }

    @DataProvider(name="3dArray")
    public Object[][] threeDArrayValues() {
        return new Object[][]{
                {
                        new double[][][]{{{100, 200, 300}, {10, 20, 30}}, {{400, 500, 600}, {40, 50, 60}}},
                        new double[][] {{500, 700, 900}, {50, 70, 90}}
                }
        };
    }

    @Test(dataProvider = "2dArray")
    public void testMax2DIndex(double [][] array, Pair<Integer, Integer> gt) {
        Assert.assertEquals(CNLOHCaller.max2dIndices(array), gt);
    }

    @DataProvider(name="2dArray")
    public Object[][] twoDArrayValues() {
        return new Object[][]{
                {
                        new double[][] {{500, 700, 900}, {50, 70, 90}}, new Pair<>(0,2)
                }
        };
    }

    @Test(dataProvider = "crMafSegDists")
    public void testCalcE_zsk_vsm_wsn(final double mafMode, final double mafLo, final double mafHi,
                                      final double crMode, final double crLow, final double crHigh,
                                      final double lambda, final double gt) {
        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(ACNV_SEG_FILE);
        AllelicBalanceCallerModelState state = AllelicBalanceCallerModelState.createInitialCNLOHCallerModelState(0.2, segments,
                HomoSapiensConstants.DEFAULT_PLOIDY, CNLOHCaller.NUM_RHOS);

        final CNLOHCaller cnlohCaller = new CNLOHCaller();
        final double[][][] responsibilities = cnlohCaller.calculateResponsibilities(state.getEffectivePhis(), state.getEffectivePis(), state.getRhos(),
                mafMode, mafLo, mafHi, crMode, crLow, crHigh, lambda, state.getmVals(), state.getnVals());

        // Will be slightly less than 1.0, but should be pretty close.  rho == 0, M == N == 1
        Assert.assertEquals(responsibilities[0][1][1], gt, 1e-4);
    }

    @DataProvider(name="crMafSegDists")
    public Object[][] crMafSegDists() {
        return new Object[][]{
                // maf mode, maf lo, maf hi, cr mode, cr low, cr high, lambda, gt
                {
                        0.499, 0.48, 0.4999, 1.0, 0.9, 1.1, 2.0, 1.0
                },
                {
                        0.1, 0.09, 0.014, 1.0, 0.9, 1.1, 2.0, 0.0
                }
        };
    }

    @Test(dataProvider = "doubleGaussian")
    public void testDoubleGaussian(final double val, final double low, final double mode, final double high, final double gt) {
        Assert.assertEquals(CNLOHCaller.calculateDoubleGaussian(val, mode, low, high), gt, 1e-10);
    }

    @DataProvider(name="doubleGaussian")
    public Object[][] doubleGaussian() {
        return new Object[][] {
                // val, mode, low, high, gt
                {0, .3, .5, .6, 2.39018153097550e-05},
                {.1, .3, .5, .6, 0.00180038248042408},
                {.2, .3, .5, .6, 0.0519041698800480},
                {.3, .3, .5, .6, 0.572721254467825},
                {.4, .3, .5, .6, 2.41873300755702},
                {.5, .3, .5, .6, 7.81926869586808},
                {.6, .3, .5, .6, 1.14544250893565},
                {.7, .3, .5, .6, 0.00360076496084816},
                {.8, .3, .5, .6, 2.42901708557338e-07},
                {.9, .3, .5, .6, 3.51625759909017e-13},
                {1.0, .3, .5, .6, 1.09230800445324e-20},
        };
    }
}
