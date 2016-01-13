package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.broadinstitute.hellbender.utils.hdf5.PoNTestUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;

public class TangentNormalizerUnitTest extends BaseTest {
    private final static String TEST_DIR = "src/test/resources/org/broadinstitute/hellbender/tools/exome/";
    private final static String TEST_PCOV_FILE = TEST_DIR + "create-pon-control-full.pcov";
    private final static String TEST_FULL_NORMALS_TN_FILE = TEST_DIR + "create-pon-all-targets.pon.normal_projection";
    private final static String TEST_SOME_NORMALS_TN_FILE = TEST_DIR + "create-pon-some-targets.pon.normal_projection";
    private final static String TEST_FULL_PON = TEST_DIR + "create-pon-all-targets.pon";
    private final static String TEST_SOME_PON = TEST_DIR + "create-pon-some-targets.pon";


    @Test
    public void testSparkTangentNormalizeSparkVsNoSpark() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(new File(TEST_PCOV_FILE), 20);

        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponHDF5File);

            final TangentNormalizationResult tnWithSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);
            final TangentNormalizationResult tnWithoutSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, null);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), tnWithoutSpark.getTangentNormalized().counts(), false);
            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getPreTangentNormalized().counts(), tnWithoutSpark.getPreTangentNormalized().counts(), false);
            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentBetaHats(), tnWithoutSpark.getTangentBetaHats(), false);
        }
    }

    @Test
    public void testSparkNormalizeNormalsFullPoN() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final File ponFile = new File(TEST_FULL_PON);
        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponHDF5File);
            final RealMatrix gtTNedNormals = PoNTestUtils.readTsvIntoMatrix(new File(TEST_FULL_NORMALS_TN_FILE));
            final TangentNormalizationResult tnWithSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), gtTNedNormals, false);
        }
    }
    @Test
    public void testSparkNormalizeNormalsSomePoN() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final File ponFile = new File(TEST_SOME_PON);
        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponHDF5File);
            final RealMatrix gtTNedNormals = PoNTestUtils.readTsvIntoMatrix(new File(TEST_SOME_NORMALS_TN_FILE));
            final TangentNormalizationResult tnWithSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), gtTNedNormals, false);
        }
    }
}
