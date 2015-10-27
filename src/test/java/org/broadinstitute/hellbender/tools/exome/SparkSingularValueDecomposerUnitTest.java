package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Log;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SparkSingularValueDecomposer;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by lichtens on 10/28/15.
 */
public class SparkSingularValueDecomposerUnitTest extends BaseTest {
    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/tools/exome");
    private static final File CONTROL_PCOV_FULL_FILE = new File(TEST_FILE_DIR, "create-pon-control-full.pcov");
    private static final File CONTROL_PCOV_GT_SV = new File(TEST_FILE_DIR, "create-pon-control-full.pcov.gt.singular_values");
    private static final File CONTROL_PCOV_GT_V = new File(TEST_FILE_DIR, "create-pon-control-full.pcov.gt.v");
    private static final File CONTROL_PCOV_GT_PINV = new File(TEST_FILE_DIR, "create-pon-control-full.pcov.gt.10.pinv_truncated");

    @Test
    public void testBasicFunction(){
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        try {
            ReadCountCollection rcc = ReadCountCollectionUtils.parse(CONTROL_PCOV_FULL_FILE);
            final SVD svd = SparkSingularValueDecomposer.createSVD(ctx, rcc.counts());
            List<String> gtSingularValuesStr = FileUtils.readLines(CONTROL_PCOV_GT_SV);
            double [] gtSingularValues = gtSingularValuesStr.stream().mapToDouble(Double::parseDouble).toArray();

            Assert.assertTrue(ArrayUtils.isSameLength(gtSingularValues, svd.getSingularValues()));

            for (int i = 0; i < gtSingularValues.length; i++){
                Assert.assertEquals(svd.getSingularValues()[i], gtSingularValues[i], 1e-5, "Singular value [" + i + "] did not match ground truth within tolerance.");
            }

            List<String> lines = FileUtils.readLines(CONTROL_PCOV_GT_V);
            for (int i = 0; i < lines.size(); i++) {
                String [] tmp = lines.get(i).split("\t");
                for (int j = 0; j < tmp.length ; j ++) {

                    //Note, it is okay if one is negative of the other.  In an SVD, there are multiple possible decompositions.  So long as the norm is the same, you are good.
                    Assert.assertEquals(Math.abs(svd.getV().getEntry(i, j)), Math.abs(Double.parseDouble(tmp[j].trim())), 1e-5, "Failure in (" + i + ", " + j + ")");
                }
            }

            // Test that the pinv was calculated correctly.  Using truncated file that only has 10 lines.
            List<String> linesPinv = FileUtils.readLines(CONTROL_PCOV_GT_PINV);
            for (int i = 0; i < 10; i++) {
                String [] tmp = linesPinv.get(i).split("\t");
                for (int j = 0; j < tmp.length ; j ++) {
                    Assert.assertEquals(Math.abs(svd.getPinv().getEntry(i, j)), Math.abs(Double.parseDouble(tmp[j].trim())), 1e-5, "Failure in (" + i + ", " + j + ")");
                }
            }

        } catch (IOException ioe) {
            Assert.assertTrue(false, "Could not read test input file: " + CONTROL_PCOV_FULL_FILE);
        }

    }

}
