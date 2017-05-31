package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.pon.PoNTestUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public final class PCATangentNormalizationUtilsUnitTest extends BaseTest {
    private static final File LARGE_CNV_TEST_FILE_DIR = new File(largeFileTestDir, "cnv");
    private static final File TEST_PCOV_FILE = new File(LARGE_CNV_TEST_FILE_DIR, "create-pon-control-full.pcov");
    private static final File TEST_FULL_NORMALS_TN_FILE = new File(LARGE_CNV_TEST_FILE_DIR, "create-pon-all-targets.pon.normal_projection");
    private static final File TEST_SOME_NORMALS_TN_FILE = new File(LARGE_CNV_TEST_FILE_DIR, "create-pon-some-targets.pon.normal_projection");
    private static final File TEST_FULL_PON = new File(LARGE_CNV_TEST_FILE_DIR, "create-pon-all-targets.pon");
    private static final File TEST_SOME_PON = new File(LARGE_CNV_TEST_FILE_DIR, "create-pon-some-targets.pon");

    @DataProvider(name="normalizeReadCountByTargetFactorsData")
    public Object[][] normalizeReadCountByTargetFactorsData() {
        final List<Object[]> result = new ArrayList<>(1);
        @SuppressWarnings("serial")
        final List<Target> targets = new ArrayList<Target>() {{
            add(new Target("A"));
            add(new Target("B"));
            add(new Target("C"));
        }};

        @SuppressWarnings("serial")
        final List<String> columnNames = new ArrayList<String>() {{
            add("1");
            add("2");
            add("3");
        }};
        result.add(new Object[] {
                new ReadCountCollection(targets, columnNames,
                        new Array2DRowRealMatrix(new double[][] {
                                new double[] { 1.1, 2.2, 3.3 },
                                new double[] { 0.1, 0.2, 0.3},
                                new double[] { 11.1, 22.2, 33.3 } }, false))
                , new double[] { 100.0, 200.0, 300.0 }});
        return result.toArray(new Object[1][]);
    }

    @Test(dataProvider="normalizeReadCountByTargetFactorsData")
    public void testNormalizeReadCountByTargetFactors(final ReadCountCollection readCount, final double[] targetFactors) {
        final RealMatrix before = readCount.counts().copy();
        PCATangentNormalizationUtils.factorNormalize(readCount.counts(), targetFactors);
        final RealMatrix after = readCount.counts();
        Assert.assertEquals(before.getColumnDimension(), after.getColumnDimension());
        Assert.assertEquals(before.getRowDimension(), after.getRowDimension());
        for (int i = 0; i < after.getColumnDimension(); i++) {
            for (int j = 0; j < after.getRowDimension(); j++) {
                final double beforeValue = before.getEntry(j, i);
                final double afterValue = after.getEntry(j, i);
                Assert.assertEquals(afterValue, beforeValue / targetFactors[j], 0.001);
            }
        }
    }

    @Test
    public void testSparkTangentNormalizeSparkVsNoSpark() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final File ponFile = PoNTestUtils.createDummyHDF5FilePoN(TEST_PCOV_FILE, 20);

        try (final HDF5File ponHDF5File = new HDF5File(ponFile)) {
            final PCACoveragePoN pon = new HDF5PCACoveragePoN(ponHDF5File);

            final PCATangentNormalizationResult tnWithSpark = pon.normalizeNormalsInPoN(ctx);
            final PCATangentNormalizationResult tnWithoutSpark = pon.normalizeNormalsInPoN();

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), tnWithoutSpark.getTangentNormalized().counts(), false);
            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getPreTangentNormalized().counts(), tnWithoutSpark.getPreTangentNormalized().counts(), false);
            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentBetaHats(), tnWithoutSpark.getTangentBetaHats(), false);
        }
    }

    @Test
    public void testSparkNormalizeNormalsFullPoN() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        try (final HDF5File ponHDF5File = new HDF5File(TEST_FULL_PON)) {
            final PCACoveragePoN pon = new HDF5PCACoveragePoN(ponHDF5File);
            final RealMatrix gtTNedNormals = PoNTestUtils.readTsvIntoMatrix(TEST_FULL_NORMALS_TN_FILE);
            final PCATangentNormalizationResult tnWithSpark = pon.normalizeNormalsInPoN(ctx);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), gtTNedNormals, false);
        }
    }
    @Test
    public void testSparkNormalizeNormalsSomePoN() {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        try (final HDF5File ponHDF5File = new HDF5File(TEST_SOME_PON)) {
            final PCACoveragePoN pon = new HDF5PCACoveragePoN(ponHDF5File);
            final RealMatrix gtTNedNormals = PoNTestUtils.readTsvIntoMatrix(TEST_SOME_NORMALS_TN_FILE);
            final PCATangentNormalizationResult tnWithSpark = pon.normalizeNormalsInPoN(ctx);

            PoNTestUtils.assertEqualsMatrix(tnWithSpark.getTangentNormalized().counts(), gtTNedNormals, false);
        }
    }
}
