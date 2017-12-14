package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.stream.IntStream;


public final class SingularValueDecomposerUnitTest extends GATKBaseTest {

    @DataProvider(name="makeSVDers")
    private Object[][] makeSVDers(){
        final RealMatrix m1  = new Array2DRowRealMatrix(new double[][]{{ 2.0, 0.0 }, { 0.0, -3.0 }, { 0.0, 0.0 }});
        return new Object[][]{
                {new ApacheSingularValueDecomposer(), m1},
                {new OjAlgoSingularValueDecomposer(), m1},
                {new SparkSingularValueDecomposer(SparkContextFactory.getTestSparkContext()), m1},
        };
    }
    @Test(dataProvider = "makeSVDers")
    public void testBasicFunction(final SingularValueDecomposer svder, final RealMatrix m){
        final SVD svd = svder.createSVD(m);
        assertSVD(svd, m);
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

    public void assertSVD(final SVD svd, final RealMatrix m) {
        final RealMatrix U = svd.getU();
        final RealMatrix V = svd.getV();
        final double[] s = svd.getSingularValues();

        final RealMatrix S = new Array2DRowRealMatrix(U.getColumnDimension(), V.getColumnDimension());
        IntStream.range(0, s.length).forEach(n -> S.setEntry(n, n, s[n]));

        final RealMatrix SPrime = new Array2DRowRealMatrix(V.getColumnDimension(), U.getColumnDimension());
        IntStream.range(0, s.length).forEach(n -> SPrime.setEntry(n, n, 1/s[n]));

        Assert.assertEquals(U.multiply(S).multiply(V.transpose()), m);
        Assert.assertEquals(V.transpose().multiply(V), MatrixUtils.createRealIdentityMatrix(2), "V is unitary");
        Assert.assertEquals(U.transpose().multiply(U), MatrixUtils.createRealIdentityMatrix(2), "U is unitary");

        Assert.assertTrue(isUnitaryMatrix(U), "matrix U");
        Assert.assertTrue(isUnitaryMatrix(V), "matrix V");
        Assert.assertEquals(U.multiply(S).multiply(V.transpose()), m);
        Assert.assertEquals(V.multiply(SPrime).multiply(U.transpose()), svd.getPinv());
        assertPseudoinverse(m, svd.getPinv());
    }

    private void assertPseudoinverse(final RealMatrix m, final RealMatrix p) {
        Assert.assertEquals(m.multiply(p).multiply(m), m);
        Assert.assertEquals(p.multiply(m).multiply(p), p);
        Assert.assertEquals(m.multiply(p).transpose(), m.multiply(p));
        Assert.assertEquals(p.multiply(m).transpose(), p.multiply(m));
    }
}

