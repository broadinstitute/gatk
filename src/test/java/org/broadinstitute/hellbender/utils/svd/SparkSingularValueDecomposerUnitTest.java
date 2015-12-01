package org.broadinstitute.hellbender.utils.svd;

import htsjdk.samtools.util.Log;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;


public class SparkSingularValueDecomposerUnitTest extends BaseTest {
    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File CONTROL_PCOV_FULL_FILE = new File(TEST_FILE_DIR, "create-pon-control-full.pcov");
    private static final File CONTROL_PCOV_GT_SV = new File(TEST_FILE_DIR, "create-pon-control-full.pcov.gt.singular_values");
    private static final File CONTROL_PCOV_GT_V = new File(TEST_FILE_DIR, "create-pon-control-full.pcov.gt.v");
    private static final File CONTROL_PCOV_GT_PINV = new File(TEST_FILE_DIR, "create-pon-control-full.pcov.gt.10.pinv_truncated");

    @Test
    public void testBasicFunction(){
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        try {
            final ReadCountCollection rcc = ReadCountCollectionUtils.parse(CONTROL_PCOV_FULL_FILE);
            final SVD svd = SparkSingularValueDecomposer.createSVD(ctx, rcc.counts());
            SVDTestUtils.assertSVD(svd, CONTROL_PCOV_GT_SV, CONTROL_PCOV_GT_V, CONTROL_PCOV_GT_PINV);
        } catch (final IOException ioe) {
            Assert.fail("Error in test data.", ioe);
        }
    }
}
