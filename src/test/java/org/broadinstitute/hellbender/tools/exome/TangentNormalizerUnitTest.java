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
    private final static File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private final static File TEST_PCOV_FILE = new File(TEST_DIR, "create-pon-control-full.pcov");
    private final static File TEST_FULL_NORMALS_TN_FILE = new File(TEST_DIR, "create-pon-all-targets.pon.normal_projection");
    private final static File TEST_SOME_NORMALS_TN_FILE = new File(TEST_DIR, "create-pon-some-targets.pon.normal_projection");
    private final static File TEST_FULL_PON = new File(TEST_DIR, "create-pon-all-targets.pon");
    private final static File TEST_SOME_PON = new File(TEST_DIR, "create-pon-some-targets.pon");

    @Test
    public void testSparkTangentNormalizeSparkVsNoSpark() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(TEST_PCOV_FILE, 20);

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

        try (final HDF5File ponHDF5File = new HDF5File(TEST_FULL_PON)) {
            final PoN pon = new HDF5PoN(ponHDF5File);
            final RealMatrix gtTNedNormals = PoNTestUtils.readTsvIntoMatrix(TEST_FULL_NORMALS_TN_FILE);
            final TangentNormalizationResult tnWithSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), gtTNedNormals, false);
        }
    }
    @Test
    public void testSparkNormalizeNormalsSomePoN() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        try (final HDF5File ponHDF5File = new HDF5File(TEST_SOME_PON)) {
            final PoN pon = new HDF5PoN(ponHDF5File);
            final RealMatrix gtTNedNormals = PoNTestUtils.readTsvIntoMatrix(TEST_SOME_NORMALS_TN_FILE);
            final TangentNormalizationResult tnWithSpark = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), gtTNedNormals, false);
        }
    }
}
