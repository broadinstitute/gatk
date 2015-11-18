package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * A class to hold some static methods useful in unit testing the SVD functionality.
 */
final public class SVDTestUtils {

    /**
     * Assrt that the given matrix is unitary.
     * @param m
     */
    public static void assertUnitaryMatrix(final RealMatrix m){
        final RealMatrix mInv = new QRDecomposition(m).getSolver().getInverse();
        final RealMatrix mT = m.transpose();

        for (int i = 0; i < mInv.getRowDimension(); i ++) {
            for (int j = 0; j < mInv.getColumnDimension(); j ++) {
                Assert.assertEquals(mInv.getEntry(i, j), mT.getEntry(i, j), 1e-7);
            }
        }
    }

    /**
     *  Run several checks on a SVD for consistency and value matching.
     *
     * @param svd SVD to test
     * @param svGT ground truth values for the singular values.  These should be in a text file (plain text) with one value on each line.  No header.
     * @param vGT ground truth values for the V matrix.  These should be in a text file (plain text) with tab separated values on each line.  No header.
     * @param pinvGT ground truth values for the psuedoinverse matrix.  These should be in a text file (plain text) with tab separated values on each line.  No header.
     */
    public static void assertSVD(final SVD svd, final File svGT, final File vGT, final File pinvGT) {
        try {
            final List<String> gtSingularValuesStr = FileUtils.readLines(svGT);
            final double [] gtSingularValues = gtSingularValuesStr.stream().mapToDouble(Double::parseDouble).toArray();

            Assert.assertTrue(ArrayUtils.isSameLength(gtSingularValues, svd.getSingularValues()));

            for (int i = 0; i < gtSingularValues.length; i++){
                Assert.assertEquals(svd.getSingularValues()[i], gtSingularValues[i], 1e-5, "Singular value [" + i + "] did not match ground truth within tolerance.");
            }

            final List<String> lines = FileUtils.readLines(vGT);
            for (int i = 0; i < lines.size(); i++) {
                final String [] tmp = lines.get(i).split("\t");
                for (int j = 0; j < tmp.length ; j ++) {

                    //Note, it is okay if one is negative of the other.  In an SVD, there are multiple possible decompositions.  So long as the norm is the same, you are good.
                    Assert.assertEquals(Math.abs(svd.getV().getEntry(i, j)), Math.abs(Double.parseDouble(tmp[j].trim())), 1e-5, "Failure in (" + i + ", " + j + ")");
                }
            }

            // Test that the pinv was calculated correctly.  Using truncated file that only has 10 lines.
            final List<String> linesPinv = FileUtils.readLines(pinvGT);
            for (int i = 0; i < 10; i++) {
                final String [] tmp = linesPinv.get(i).split("\t");
                for (int j = 0; j < tmp.length ; j ++) {
                    Assert.assertEquals(Math.abs(svd.getPinv().getEntry(i, j)), Math.abs(Double.parseDouble(tmp[j].trim())), 1e-5, "Failure in (" + i + ", " + j + ")");
                }
            }

            SVDTestUtils.assertUnitaryMatrix(svd.getV());
            SVDTestUtils.assertUnitaryMatrix(svd.getU());
        } catch (final IOException ioe) {
            Assert.fail("Error in test data.", ioe);
        }
    }
}
