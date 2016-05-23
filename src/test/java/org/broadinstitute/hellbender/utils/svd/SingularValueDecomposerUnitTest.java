package org.broadinstitute.hellbender.utils.svd;

import htsjdk.samtools.util.Log;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


public final class SingularValueDecomposerUnitTest extends BaseTest {

    @DataProvider(name="makeSVDers")
    private Object[][] makeSVDers(){
        return new Object[][]{
                {new ApacheSingularValueDecomposer()},
                {new OjAlgoSingularValueDecomposer()},
                {new SparkSingularValueDecomposer(SparkContextFactory.getTestSparkContext())}
        };
    }
    @Test(dataProvider = "makeSVDers")
    public void testBasicFunction(final SingularValueDecomposer svder){
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final Array2DRowRealMatrix m = new Array2DRowRealMatrix(new double[][]{{ 2.0, 0.0 }, { 0.0, -3.0 }, { 0.0, 0.0 }});

        final double[] expectedSingularValues = {3.0, 2.0};

        final SVD svd = svder.createSVD(m);
        Assert.assertEquals(svd.getU().multiply(MatrixUtils.createRealDiagonalMatrix(svd.getSingularValues())).multiply(svd.getV().transpose()), m);
        Assert.assertTrue(isUnitaryMatrix(svd.getV()), "V is unitary");
        Assert.assertTrue(isUnitaryMatrix(svd.getU()), "U is unitary");
        Assert.assertEquals(svd.getV().transpose().multiply(svd.getV()), MatrixUtils.createRealIdentityMatrix(2), "V is unitary");
        Assert.assertEquals(svd.getU().transpose().multiply(svd.getU()), MatrixUtils.createRealIdentityMatrix(2), "U is unitary");
        Assert.assertEquals(svd.getSingularValues(), expectedSingularValues, "singular value");
    }

    /**
     * Check that the given matrix is unitary.
     */
    public static boolean isUnitaryMatrix(final RealMatrix m){
        //Note can't use MatrixUtils.inverse because m may not be square
        final RealMatrix mInv = new QRDecomposition(m).getSolver().getInverse();
        final RealMatrix mT = m.transpose();

        for (int i = 0; i < mInv.getRowDimension(); i ++) {
            for (int j = 0; j < mInv.getColumnDimension(); j ++) {
                if (Math.abs(mInv.getEntry(i, j) - mT.getEntry(i, j)) > 1.0e-7){
                    return false;
                }
            }
        }
        return true;
    }
}

