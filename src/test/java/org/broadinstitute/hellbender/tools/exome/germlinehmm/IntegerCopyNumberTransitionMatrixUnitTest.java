package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberTransitionMatrixUnitTest extends BaseTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/germlinehmm";
    private static final File HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TABLE_FILE = new File(TEST_SUB_DIR,
            "TCGA_T_matrix_autosomal.tsv");
    private static final File HOMO_SAPIENS_COPY_NUMBER_TRANSITION_BAD_AUTOSOMAL_TABLE_FILE = new File(TEST_SUB_DIR,
            "TCGA_T_matrix_autosomal_bad.tsv");

    /**
     * Copied directly from the corresponding tsv file
     */
    protected static final RealMatrix HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TRUTH =
            new Array2DRowRealMatrix(new double[][]{
                    {0.9997389770672177, 7.032645704368021e-07, 0.00026031966821188487, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0},
                    {2.1467075095351557e-07, 0.99981467801052126, 0.00018510163208542175, 5.6866424093646516e-09,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                    {5.9100196666398515e-08, 1.376696014985822e-07, 0.99999972823891037, 3.4570453537193609e-08,
                            3.3588164155451664e-08, 2.8971721269891569e-09, 2.1348140599967544e-09,
                            7.8350530879524278e-10, 9.4105226022640503e-10, 5.9212277047953571e-11, 4.2294483605681122e-12},
                    {0.0, 1.1958483005083788e-08, 9.7745650462803608e-05, 0.99990217064015618, 5.6802794274147987e-08,
                            1.1958483005083788e-08, 2.989620751270947e-09, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.00010292959512329779, 6.1564638523662337e-08, 0.99989682090607845, 1.2960976531297335e-08,
                            1.7497318317251403e-07, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 8.0138884782014782e-05, 1.1699107267447413e-07, 1.1699107267447413e-07,
                            0.99991953938976808, 8.7743304505855591e-08, 0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 9.0232386606850211e-05, 4.469162288600803e-08, 2.4133476358444333e-06,
                            1.3407486865802407e-07, 0.99990713080764293, 4.469162288600803e-08, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 7.9777767777145916e-05, 0.0, 0.0, 0.0, 1.0766230469250462e-07, 0.99992011456991814,
                            0.0, 0.0, 0.0},
                    {0.0, 0.0, 7.5782475345150855e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99992421752465488, 0.0, 0.0},
                    {0.0, 0.0, 8.8396952830754563e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9999116030471692, 0.0},
                    {0.0, 0.0, 8.9130531663621367e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99991086946833641}}).transpose();

    /**
     * Copied directly from the corresponding tsv file
     */
    protected static final RealMatrix HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_X_TRUTH =
            new Array2DRowRealMatrix(new double[][] {
                    {0.99966751443861601, 0.0, 0.00033248556138398773, 0.0, 0.0},
                    {0.0, 0.9997899779920747, 0.00021002200792527603, 0.0, 0.0},
                    {4.6641935242897276e-08, 9.3423238818193651e-08, 0.99999985473174158, 5.0172599663674365e-09, 1.8582444319879394e-10},
                    {0.0, 0.0, 4.541929579905158e-05, 0.99995458070420096, 0.0},
                    {0.0, 0.0, 7.8833267638943636e-05, 0.0, 0.99992116673236109}}).transpose();

    /**
     * Copied directly from the corresponding tsv file
     */
    protected static final RealMatrix HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_Y_TRUTH =
            new Array2DRowRealMatrix(new double[][] {{1.0}});

    /**
     * Copied directly from the corresponding tsv file
     */
    protected static final RealMatrix HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_X_TRUTH =
            new Array2DRowRealMatrix(new double[][] {
                    {0.99971173098797461, 0.00028826901202540391, 0.0, 0.0},
                    {1.0067714836234777e-07, 0.99999989259309574, 6.5615120089096975e-09, 1.6824389766435122e-10},
                    {0.0, 7.456796468461193e-05, 0.99992504963549644, 3.8239981889544576e-07},
                    {0.0, 4.0420371867421184e-05, 8.0840743734842364e-06, 0.99995149555375906}
            }).transpose();

    /**
     * Copied directly from the corresponding tsv file
     */
    protected static final RealMatrix HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_Y_TRUTH =
            new Array2DRowRealMatrix(new double[][] {
                    {0.99966851990709416, 0.00033148009290586881, 0.0, 0.0, 0.0},
                    {5.9399783434370542e-08, 0.99999937404917871, 5.2251326738303193e-07, 3.4001255345191416e-08, 1.0036515132014333e-08},
                    {0.0, 0.00016831138093714932, 0.99983149068329535, 1.9793576746822735e-07, 0.0},
                    {0.0, 0.00035529148884256304, 6.4209305212511401e-06, 0.99963614727046246, 2.1403101737503797e-06},
                    {0.0, 0.00027047913446676971, 0.0, 5.519982336056525e-06, 0.99972400088319713}
            }).transpose();

    @Test
    public void testReadCopyNumberTransitionMatrixTable() {
        final IntegerCopyNumberTransitionMatrix loadedTransitionMatrix =
                IntegerCopyNumberTransitionMatrix.read(
                        HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TABLE_FILE, 0);
        MathObjectAsserts.assertRealMatrixEquals(loadedTransitionMatrix.getTransitionMatrix(),
                HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TRUTH);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadCopyNumberTransitionMatrixTableBadInput() {
        IntegerCopyNumberTransitionMatrix.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_BAD_AUTOSOMAL_TABLE_FILE, 0);
    }

    @Test
    public void testUnnormalizedProbability() {
        /* it should normalize unnormalized transition matrices and give a warning */
        final IntegerCopyNumberTransitionMatrix transitionMatrix = new IntegerCopyNumberTransitionMatrix(
                new Array2DRowRealMatrix(new double[][]{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}), 0);
        for (int i=0; i<3; i++) {
            final double[] col = transitionMatrix.getTransitionMatrix().getColumn(i);
            Assert.assertEquals(Arrays.stream(col).sum(), 1.0, 1e-12);
        }
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNegativeProbability() {
        new IntegerCopyNumberTransitionMatrix(
                new Array2DRowRealMatrix(new double[][]{{1, 2, 3}, {4, 5, 6}, {7, 8, -9}}), 0);
    }

    @Test
    public void testPadding() {
        final IntegerCopyNumberTransitionMatrix data = new IntegerCopyNumberTransitionMatrix(
                new Array2DRowRealMatrix(new double[][] {
                        {1, 2, 3},
                        {4, 5, 6},
                        {7, 8, 9}}), 2);
        final RealMatrix expected = new Array2DRowRealMatrix(
                new double[][]{
                        {1.0/12, 2.0/15, 3.0/18, 0, 0},
                        {4.0/12, 5.0/15, 6.0/18, 0, 0},
                        {7.0/12, 8.0/15, 9.0/18, 0, 0},
                        {0,      0,      0,      1, 0},
                        {0,      0,      0,      0, 1}});
        Assert.assertEquals(data.getTransitionMatrix().subtract(expected).getNorm(), 0, 1e-12);
    }
}
