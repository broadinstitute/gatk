package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.mockito.internal.util.collections.Sets;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class DepthMatrixTest extends GATKBaseTest {

    private static final String SAMPLE_1 = "sample1";
    private static final String SAMPLE_2 = "sample2";
    private static final String SAMPLE_3 = "sample3";

    private static final double RAW_DEPTH_1 = 29.;
    private static final double RAW_DEPTH_2 = 30.;
    private static final double RAW_DEPTH_3 = 33.;

    private static final double MEDIAN_1 = 32.;
    private static final double MEDIAN_2 = 35.;
    private static final double MEDIAN_3 = 40.;

    private static final double TEST_TOLERANCE = 1e-6;

    private static final String DEPTH_FILE_PATH = largeFileTestDir + "sv/aggregate_depth_test.rd.txt.gz";

    // Based on all_samples_depth_chr2_00000675 variant from AggregateDepthEvidence integration test
    private static final SimpleInterval TEST_INTERVAL = new SimpleInterval("chr2", 89295001, 89306001);

    // Based on all_samples_depth_chr1_00002cd3 variant from AggregateDepthEvidence integration test
    private static final SimpleInterval TEST_INTERVAL_LARGE = new SimpleInterval("chr1", 143503901, 143968588);

    // Based on all_samples_depth_chr1_00000141 variant from AggregateDepthEvidence integration test, includes gap
    private static final SimpleInterval TEST_INTERVAL_GAP = new SimpleInterval("chr1", 206001, 261667);

    private static DepthMatrix makeTestMatrix(final String contig, final int start, final int end, final int binSize) {
        Utils.validateArg(binSize > 0, "binSize must be greater than zero");
        final List<SimpleInterval> bins = new ArrayList<>();
        for (int i = start; i < end; i += binSize) {
            bins.add(new SimpleInterval(contig, i, i + binSize - 1));
        }
        final int numBins = bins.size();
        final double[] counts1 = new double[numBins];
        final double[] counts2 = new double[numBins];
        final double[] counts3 = new double[numBins];
        for (int i = 0; i < numBins; i++) {
            counts1[i] = RAW_DEPTH_1;
            counts2[i] = RAW_DEPTH_2;
            counts3[i] = RAW_DEPTH_3;
        }
        final Map<String, double[]> matrix = new HashMap<>();
        matrix.put(SAMPLE_1, counts1);
        matrix.put(SAMPLE_2, counts2);
        matrix.put(SAMPLE_3, counts3);
        return new DepthMatrix(bins, matrix);
    }

    private static Map<String, Double> makeSampleMedians() {
        final Map<String, Double> medians = new HashMap<>();
        medians.put(SAMPLE_1, MEDIAN_1);
        medians.put(SAMPLE_2, MEDIAN_2);
        medians.put(SAMPLE_3, MEDIAN_3);
        return medians;
    }

    private static Map<String, Double> makeSampleMedians(final List<String> sampleNames, final double median) {
        final Map<String, Double> medians = new HashMap<>();
        for (final String sampleName : sampleNames) {
            medians.put(sampleName, median);
        }
        return medians;
    }

    private static void assertMatricesEqual(final DepthMatrix mat1, final DepthMatrix mat2) {
        Assert.assertEquals(mat1.getNumBins(), mat2.getNumBins());
        Assert.assertEquals(mat1.getSampleSet(), mat2.getSampleSet());
        for (final String sample : mat1.getSampleSet()) {
            Assert.assertEquals(mat1.getSample(sample), mat2.getSample(sample));
        }
    }

    @Test
    public void testLoadTiny() {
        final int bins = 10;
        final int largeSizeCutoff = 100000; // Larger than interval size
        final int largeVariantRegionPoints = 15;
        final int largeVariantWindow = 1000;
        final FeatureDataSource<DepthEvidence> source = new FeatureDataSource<>(new File(DEPTH_FILE_PATH));
        final DepthMatrixLoader loader = new DepthMatrixLoader(source, bins, largeSizeCutoff, largeVariantRegionPoints, largeVariantWindow, SVTestUtils.hg38Dict);
        final Map<String, Double> sampleMedians = makeSampleMedians(((SVFeaturesHeader) source.getHeader()).getSampleNames(), 30.);
        final SimpleInterval interval = new SimpleInterval("chr2", 89295001, 89295002);
        final DepthMatrix matrix = loader.load(interval, sampleMedians);
        Assert.assertEquals(matrix.getNumBins(), 1);
        Assert.assertEquals(matrix.getSampleSet().size(), 156);
        Assert.assertEquals(matrix.getBins().get(0), new SimpleInterval("chr2", 89294920, 89295019));
        final double[] values = matrix.getSample("NA20845");
        SVTestUtils.assertFloatWithinTolerance(values[0], 0.5666666666666667, TEST_TOLERANCE);
    }

    @Test
    public void testLoad() {
        final int bins = 10;
        final int largeSizeCutoff = 100000; // Larger than interval size
        final int largeVariantRegionPoints = 15;
        final int largeVariantWindow = 1000;
        final FeatureDataSource<DepthEvidence> source = new FeatureDataSource<>(new File(DEPTH_FILE_PATH));
        final DepthMatrixLoader loader = new DepthMatrixLoader(source, bins, largeSizeCutoff, largeVariantRegionPoints, largeVariantWindow, SVTestUtils.hg38Dict);
        final Map<String, Double> sampleMedians = makeSampleMedians(((SVFeaturesHeader)source.getHeader()).getSampleNames(), 30.);
        final DepthMatrix matrix = loader.load(TEST_INTERVAL, sampleMedians);
        Assert.assertEquals(matrix.getNumBins(), 8);
        Assert.assertEquals(matrix.getSampleSet().size(), 156);
        Assert.assertEquals(matrix.getBins().get(0), new SimpleInterval("chr2", 89296420, 89297419));
        Assert.assertEquals(matrix.getBins().get(7), new SimpleInterval("chr2", 89303420, 89304419));
        final double[] values = matrix.getSample("NA20845");
        SVTestUtils.assertFloatWithinTolerance(values[0], 0.5676567656765676, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[1], 0.594059405940594, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[2], 0.47854785478547857, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[3], 0.42244224422442245, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[4], 0.7524752475247525, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[5], 0.43564356435643564, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[6], 0.43234323432343236, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[7], 0.3795379537953795, TEST_TOLERANCE);
    }

    @Test
    public void testLoadWithGap() {
        final int bins = 10;
        final int largeSizeCutoff = 100000;
        final int largeVariantRegionPoints = 15;
        final int largeVariantWindow = 1000;
        final FeatureDataSource<DepthEvidence> source = new FeatureDataSource<>(new File(DEPTH_FILE_PATH));
        final DepthMatrixLoader loader = new DepthMatrixLoader(source, bins, largeSizeCutoff, largeVariantRegionPoints, largeVariantWindow, SVTestUtils.hg38Dict);
        final Map<String, Double> sampleMedians = makeSampleMedians(((SVFeaturesHeader)source.getHeader()).getSampleNames(), 30.);
        final DepthMatrix matrix = loader.load(TEST_INTERVAL_GAP, sampleMedians);
        Assert.assertEquals(matrix.getNumBins(), 8);
        Assert.assertEquals(matrix.getSampleSet().size(), 156);
        Assert.assertEquals(matrix.getBins().get(0), new SimpleInterval("chr1", 211800, 217300));
        Assert.assertEquals(matrix.getBins().get(7), new SimpleInterval("chr1", 250300, 255800));
        final double[] values = matrix.getSample("NA20845");
        SVTestUtils.assertFloatWithinTolerance(values[0], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[1], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[2], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[3], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[4], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[5], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[6], 0.03327283726557774, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[7], 0.03327283726557774, TEST_TOLERANCE);
    }

    @Test
    public void testLoadLarge() {
        final int bins = 10;
        final int largeSizeCutoff = 100000; // Smaller than interval size
        final int largeVariantRegionPoints = 15;
        final int largeVariantWindow = 1000;
        final FeatureDataSource<DepthEvidence> source = new FeatureDataSource<>(new File(DEPTH_FILE_PATH));
        final DepthMatrixLoader loader = new DepthMatrixLoader(source, bins, largeSizeCutoff, largeVariantRegionPoints, largeVariantWindow, SVTestUtils.hg38Dict);
        final Map<String, Double> sampleMedians = makeSampleMedians(((SVFeaturesHeader)source.getHeader()).getSampleNames(), 30.);
        final DepthMatrix matrix = loader.load(TEST_INTERVAL_LARGE, sampleMedians);
        Assert.assertEquals(matrix.getNumBins(), 8);
        Assert.assertEquals(matrix.getSampleSet().size(), 156);
        Assert.assertEquals(matrix.getBins().get(0), new SimpleInterval("chr1", 143596838, 143597838));
        Assert.assertEquals(matrix.getBins().get(7), new SimpleInterval("chr1", 143813691, 143814691));
        final double[] values = matrix.getSample("NA20845");
        SVTestUtils.assertFloatWithinTolerance(values[0], 0.6363636363636364, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[1], 0.6666666666666666, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[2], 0.6666666666666666, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[3], 0.6363636363636364, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[4], 0.7272727272727273, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[5], 0.6060606060606061, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[6], 0.6666666666666666, TEST_TOLERANCE);
        SVTestUtils.assertFloatWithinTolerance(values[7], 0.6666666666666666, TEST_TOLERANCE);
    }

    @Test
    public void testSliceEmpty() {
        final DepthMatrix matrix = new DepthMatrix(Collections.singletonList(TEST_INTERVAL), Collections.emptyMap());
        final double[] slice = matrix.slice(3, Collections.emptySet());
        Assert.assertNotNull(slice);
        Assert.assertEquals(slice.length, 0);
    }

    @Test
    public void testSlice() {
        final DepthMatrix matrix = makeTestMatrix("chr1", 10000, 20000, 100);
        final double[] rawValues1 = matrix.getSample(SAMPLE_1);
        final double[] rawValues2 = matrix.getSample(SAMPLE_2);
        final double[] rawValues3 = matrix.getSample(SAMPLE_3);
        for (int i = 0; i < matrix.getBins().size(); i++) {
            rawValues1[i] = RAW_DEPTH_1 + i;
            rawValues2[i] = RAW_DEPTH_2 + i;
            rawValues3[i] = RAW_DEPTH_3 + i;
        }
        final double[] slice = matrix.slice(20, Sets.newSet(SAMPLE_1, SAMPLE_2));
        Assert.assertNotNull(slice);
        Assert.assertEquals(slice.length, 2);
        Assert.assertTrue(slice[0] == RAW_DEPTH_1 + 20 || slice[1] == RAW_DEPTH_1 + 20);
        Assert.assertTrue(slice[0] == RAW_DEPTH_2 + 20 || slice[1] == RAW_DEPTH_2 + 20);
    }

    @Test
    public void testCompressMatrix() {
        final DepthMatrix rawMatrix = makeTestMatrix("chr1", 1000, 4999, 100);
        final double[] rawValues1 = rawMatrix.getSample(SAMPLE_1);
        final double[] rawValues2 = rawMatrix.getSample(SAMPLE_2);
        final double[] rawValues3 = rawMatrix.getSample(SAMPLE_3);
        for (int i = 0; i < rawMatrix.getBins().size(); i++) {
            rawValues1[i] = RAW_DEPTH_1 + i;
            rawValues2[i] = RAW_DEPTH_2 + i;
            rawValues3[i] = RAW_DEPTH_3 + i;
        }
        final DepthMatrixLoader.CompressionResult compressionResult = new DepthMatrixLoader.CompressionResult(1.5, new SimpleInterval("chr1", 1150, 4750));
        final DepthMatrix compressedMatrix = DepthMatrixLoader.compressMatrix(rawMatrix, compressionResult);
        Assert.assertEquals(compressedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(compressedMatrix.getBins().size(), 26);
        final SimpleInterval firstBin = compressedMatrix.getBins().get(0);
        final SimpleInterval thirdBin = compressedMatrix.getBins().get(2);
        final SimpleInterval lastBin = compressedMatrix.getBins().get(compressedMatrix.getBins().size() - 1);
        Assert.assertEquals(firstBin.getStart(), 1000);
        Assert.assertEquals(firstBin.getEnd(), 1099);
        Assert.assertEquals(thirdBin.getStart(), 1300);
        Assert.assertEquals(thirdBin.getEnd(), 1399);
        Assert.assertEquals(lastBin.getStart(), 4700);
        Assert.assertEquals(lastBin.getEnd(), 4799);
        final double[] values1 = compressedMatrix.getSample(SAMPLE_1);
        Assert.assertEquals(values1.length, compressedMatrix.getBins().size());
        Assert.assertEquals(values1[0], 29.);
        Assert.assertEquals(values1[1], 30.);
        Assert.assertEquals(values1[2], 32.);
        Assert.assertEquals(values1[3], 33.);
        Assert.assertEquals(values1[4], 35.);
        Assert.assertEquals(values1[25], 66.);
    }


    @Test
    public void testCompressMatrixNull() {
        final DepthMatrix rawMatrix = new DepthMatrix(Collections.singletonList(TEST_INTERVAL), Collections.emptyMap());
        final DepthMatrix compressedEmptyMatrix = DepthMatrixLoader.compressMatrix(rawMatrix, new DepthMatrixLoader.CompressionResult(1.5, new SimpleInterval("chr1", 1150, 4750)));
        assertMatricesEqual(rawMatrix, compressedEmptyMatrix);

        final DepthMatrix compressedEmptyMatrix2 = DepthMatrixLoader.compressMatrix(rawMatrix, new DepthMatrixLoader.CompressionResult(1., new SimpleInterval("chr1", 1150, 4750)));
        assertMatricesEqual(rawMatrix, compressedEmptyMatrix2);

        final DepthMatrix compressedEmptyMatrix3 = DepthMatrixLoader.compressMatrix(rawMatrix, null);
        assertMatricesEqual(rawMatrix, compressedEmptyMatrix3);
    }

    @DataProvider(name = "testCompressMediansData")
    public Object[][] testCompressMediansData() {
        return new Object[][]{
                {1.},
                {0.5},
                {1.5},
                {2.0},
                {11.5}
        };
    }

    @Test(dataProvider = "testCompressMediansData")
    public void testCompressMedians(final double compression) {
        final Map<String, Double> sampleMedians = makeSampleMedians();
        final DepthMatrixLoader.CompressionResult compressionResult = new DepthMatrixLoader.CompressionResult(compression, new SimpleInterval("chr1", 1150, 4750));
        final Map<String, Double> compressedMedians = DepthMatrixLoader.compressMedians(sampleMedians, compressionResult);
        Assert.assertEquals(compressedMedians.size(), sampleMedians.size());
        Assert.assertEquals(compressedMedians.keySet(), sampleMedians.keySet());
        final double compressionFactor = compression <= 1 ? 1. : compression;
        Assert.assertEquals(compressedMedians.get(SAMPLE_1), sampleMedians.get(SAMPLE_1) * compressionFactor);
        Assert.assertEquals(compressedMedians.get(SAMPLE_2), sampleMedians.get(SAMPLE_2) * compressionFactor);
        Assert.assertEquals(compressedMedians.get(SAMPLE_3), sampleMedians.get(SAMPLE_3) * compressionFactor);
    }

    @Test(dataProvider = "testCompressMediansData")
    public void testCompressMediansZero(final double compression) {
        final Map<String, Double> sampleMedians = makeSampleMedians();
        sampleMedians.put(SAMPLE_1, 0.0);
        final DepthMatrixLoader.CompressionResult compressionResult = new DepthMatrixLoader.CompressionResult(compression, new SimpleInterval("chr1", 1150, 4750));
        final Map<String, Double> compressedMedians = DepthMatrixLoader.compressMedians(sampleMedians, compressionResult);
        Assert.assertEquals(compressedMedians.size(), sampleMedians.size());
        Assert.assertEquals(compressedMedians.keySet(), sampleMedians.keySet());
        if (compression <= 1) {
            Assert.assertEquals(compressedMedians.get(SAMPLE_1), sampleMedians.get(SAMPLE_1));
            Assert.assertEquals(compressedMedians.get(SAMPLE_2), sampleMedians.get(SAMPLE_2));
            Assert.assertEquals(compressedMedians.get(SAMPLE_3), sampleMedians.get(SAMPLE_3));
        } else {
            Assert.assertEquals(compressedMedians.get(SAMPLE_1), compression);
            Assert.assertEquals(compressedMedians.get(SAMPLE_2), sampleMedians.get(SAMPLE_2) * compression);
            Assert.assertEquals(compressedMedians.get(SAMPLE_3), sampleMedians.get(SAMPLE_3) * compression);
        }
    }

    @DataProvider(name = "testCalculateCompressionData")
    public Object[][] testCalculateCompressionData() {
        return new Object[][]{
                {10001, 20000, 10, 1.0, 10001, 20000},
                {10001, 20000, 5, 2.0, 10001, 20000},
                {10001, 20000, 1, 10.0, 10001, 20000},
                {10001, 20000, 7, 1.1428571428571428, 11001, 19000},
                {10001, 20000, 20, 1.0, 10001, 20000},
                {9000, 21000, 5, 2.0, 10001, 20000},
                {11000, 19000, 5, 1.2, 12001, 18000},
                {11000, 19000, 8, 1.0, 11001, 19000},
                {14000, 16000, 4, 1.0, 14001, 16000},
                {14000, 16000, 1, 2.0, 14001, 16000},
                {14001, 15000, 2, 1.0, 10001, 20000},
        };
    }

    @Test(dataProvider = "testCalculateCompressionData")
    public void testCalculateCompression(final int start, final int end, final int bins,
                                         final double expectedCompression, final int expectedStart, final int expectedEnd) {
        final DepthMatrix matrix = makeTestMatrix("chr1", 10001, 20000, 1000);
        final SimpleInterval interval = new SimpleInterval("chr1", start, end);
        final DepthMatrixLoader.CompressionResult compressionResult = DepthMatrixLoader.calculateCompression(matrix.getBins(), interval, bins);
        Assert.assertNotNull(compressionResult);
        SVTestUtils.assertFloatWithinTolerance(compressionResult.compression(), expectedCompression, TEST_TOLERANCE);
        final SimpleInterval expectedInterval = new SimpleInterval("chr1", expectedStart, expectedEnd);
        Assert.assertEquals(compressionResult.adjustedRegion(), expectedInterval);
    }

    @DataProvider(name = "testTrimMatrixBeforeCompressionData")
    public Object[][] testTrimMatrixBeforeCompressionData() {
        return new Object[][]{
                {1, 1001, 4999, 0, 0},
                {1, 1001, 5000, 0, 0},
                {1, 1001, 5001, 0, 0},
                {1, 1000, 5001, 0, 0},
                {1, 1101, 4999, 1, 0},
                {1, 1001, 4901, 0, 1},
                {1, 1101, 4901, 1, 1},
                {0.5, 1101, 4901, 1, 1},
                {2, 1101, 4901, 1, 1},
                {1, 1102, 4900, 1, 1},
                {1, 1201, 4800, 2, 2},
        };
    }

    @Test(dataProvider = "testTrimMatrixBeforeCompressionData")
    public void testTrimMatrixBeforeCompression(final double compression, final int adjustedStart, final int adjustedEnd, final int expectedFirstOffset, final int expectedLastOffset) {
        final DepthMatrix rawMatrix = makeTestMatrix("chr1", 1001, 4999, 100);
        final double[] rawValues1 = rawMatrix.getSample(SAMPLE_1);
        final double[] rawValues2 = rawMatrix.getSample(SAMPLE_2);
        final double[] rawValues3 = rawMatrix.getSample(SAMPLE_3);
        for (int i = 0; i < rawMatrix.getBins().size(); i++) {
            rawValues1[i] = RAW_DEPTH_1 + i;
            rawValues2[i] = RAW_DEPTH_2 + i;
            rawValues3[i] = RAW_DEPTH_3 + i;
        }
        final DepthMatrix trimmedMatrix = DepthMatrixLoader.trimMatrixBeforeCompression(
                rawMatrix,
                new DepthMatrixLoader.CompressionResult(compression, new SimpleInterval("chr1", adjustedStart, adjustedEnd))
        );
        Assert.assertEquals(trimmedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(trimmedMatrix.getBins(), rawMatrix.getBins().subList(expectedFirstOffset, rawMatrix.getBins().size() - expectedLastOffset));
        final double[] values1 = trimmedMatrix.getSample(SAMPLE_1);
        final double[] values2 = trimmedMatrix.getSample(SAMPLE_2);
        final double[] values3 = trimmedMatrix.getSample(SAMPLE_3);
        final int numBins = trimmedMatrix.getBins().size();
        Assert.assertEquals(values1.length, numBins);
        Assert.assertEquals(values2.length, numBins);
        Assert.assertEquals(values3.length, numBins);
        for (int i = 0; i < numBins; i++) {
            Assert.assertEquals(values1[i], RAW_DEPTH_1 + i + expectedFirstOffset);
            Assert.assertEquals(values2[i], RAW_DEPTH_2 + i + expectedFirstOffset);
            Assert.assertEquals(values3[i], RAW_DEPTH_3 + i + expectedFirstOffset);
        }
    }

    @Test
    public void testTrimMatrixBeforeCompressionNull() {
        final DepthMatrix rawMatrix = makeTestMatrix("chr1", 1001, 4999, 100);
        final DepthMatrix trimmedMatrix = DepthMatrixLoader.trimMatrixBeforeCompression(rawMatrix, null);
        Assert.assertEquals(trimmedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(trimmedMatrix.getBins(), rawMatrix.getBins());
    }

    @Test
    public void testCalculateCompressionNull() {
        final SimpleInterval interval = new SimpleInterval("chr1", 10001, 20000);
        final DepthMatrixLoader.CompressionResult compressionResult = DepthMatrixLoader.calculateCompression(Collections.emptyList(), interval, 10);
        Assert.assertNull(compressionResult);
    }

    @Test
    public void testTrimOuterBins() {
        final DepthMatrix rawMatrix = makeTestMatrix("chr1", 1000, 4999, ((5000 - 1000) / DepthMatrixLoader.MIN_BINS_FOR_TRIMMING));
        Assert.assertTrue(rawMatrix.getBins().size() == DepthMatrixLoader.MIN_BINS_FOR_TRIMMING, "Test assertion failed, adjust bin size");
        final double[] rawValues1 = rawMatrix.getSample(SAMPLE_1);
        final double[] rawValues2 = rawMatrix.getSample(SAMPLE_2);
        final double[] rawValues3 = rawMatrix.getSample(SAMPLE_3);
        for (int i = 0; i < rawMatrix.getBins().size(); i++) {
            rawValues1[i] = RAW_DEPTH_1 + i;
            rawValues2[i] = RAW_DEPTH_2 + i;
            rawValues3[i] = RAW_DEPTH_3 + i;
        }
        final DepthMatrix trimmedMatrix = DepthMatrixLoader.trimOuterBins(rawMatrix);
        Assert.assertEquals(trimmedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(trimmedMatrix.getBins(), rawMatrix.getBins().subList(1, rawMatrix.getBins().size() - 1));
        final double[] values1 = trimmedMatrix.getSample(SAMPLE_1);
        final double[] values2 = trimmedMatrix.getSample(SAMPLE_2);
        final double[] values3 = trimmedMatrix.getSample(SAMPLE_3);
        final int numBins = trimmedMatrix.getBins().size();
        Assert.assertEquals(values1.length, numBins);
        Assert.assertEquals(values2.length, numBins);
        Assert.assertEquals(values3.length, numBins);
        for (int i = 0; i < numBins; i++) {
            Assert.assertEquals(values1[i], RAW_DEPTH_1 + i + 1);
            Assert.assertEquals(values2[i], RAW_DEPTH_2 + i + 1);
            Assert.assertEquals(values3[i], RAW_DEPTH_3 + i + 1);
        }
    }

    @Test
    public void testTrimOuterBinsNoTrim() {
        final DepthMatrix rawMatrix = makeTestMatrix("chr1", 1000, 4999, ((5000 - 1000) / (DepthMatrixLoader.MIN_BINS_FOR_TRIMMING - 1)));
        Assert.assertTrue(rawMatrix.getBins().size() < DepthMatrixLoader.MIN_BINS_FOR_TRIMMING, "Test assertion failed, adjust bin size");
        final double[] rawValues1 = rawMatrix.getSample(SAMPLE_1);
        final double[] rawValues2 = rawMatrix.getSample(SAMPLE_2);
        final double[] rawValues3 = rawMatrix.getSample(SAMPLE_3);
        for (int i = 0; i < rawMatrix.getBins().size(); i++) {
            rawValues1[i] = RAW_DEPTH_1 + i;
            rawValues2[i] = RAW_DEPTH_2 + i;
            rawValues3[i] = RAW_DEPTH_3 + i;
        }
        final DepthMatrix trimmedMatrix = DepthMatrixLoader.trimOuterBins(rawMatrix);
        Assert.assertEquals(trimmedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(trimmedMatrix.getBins(), rawMatrix.getBins());
        final double[] values1 = trimmedMatrix.getSample(SAMPLE_1);
        final double[] values2 = trimmedMatrix.getSample(SAMPLE_2);
        final double[] values3 = trimmedMatrix.getSample(SAMPLE_3);
        final int numBins = trimmedMatrix.getBins().size();
        Assert.assertEquals(values1.length, numBins);
        Assert.assertEquals(values2.length, numBins);
        Assert.assertEquals(values3.length, numBins);
        for (int i = 0; i < numBins; i++) {
            Assert.assertEquals(values1[i], RAW_DEPTH_1 + i);
            Assert.assertEquals(values2[i], RAW_DEPTH_2 + i);
            Assert.assertEquals(values3[i], RAW_DEPTH_3 + i);
        }
    }

    @Test
    public void testSetZeroesToOnes() {
        final DepthMatrix matrix = makeTestMatrix("chr1", 1000, 5000, 200);
        final double[] rawValues1 = matrix.getSample(SAMPLE_1);
        final double[] rawValues2 = matrix.getSample(SAMPLE_2);
        final double[] rawValues3 = matrix.getSample(SAMPLE_3);
        rawValues1[1] = 0;
        rawValues2[2] = 0;
        rawValues3[3] = 0;
        final Set<String> samples = new HashSet<>(matrix.getSampleSet());
        final List<SimpleInterval> intervals = new ArrayList<>(matrix.getBins());
        DepthMatrixLoader.setZeroesToOnes(matrix);
        Assert.assertEquals(matrix.getSampleSet(), samples);
        Assert.assertEquals(matrix.getBins(), intervals);
        final double[] values1 = matrix.getSample(SAMPLE_1);
        final double[] values2 = matrix.getSample(SAMPLE_2);
        final double[] values3 = matrix.getSample(SAMPLE_3);
        final int numBins = matrix.getBins().size();
        Assert.assertEquals(values1.length, numBins);
        Assert.assertEquals(values2.length, numBins);
        Assert.assertEquals(values3.length, numBins);
        Assert.assertEquals(values1[1], 1);
        Assert.assertEquals(values2[2], 1);
        Assert.assertEquals(values3[3], 1);
        for (int i = 0; i < numBins; i++) {
            if (i != 1) {
                Assert.assertEquals(values1[i], RAW_DEPTH_1);
            }
            if (i != 2) {
                Assert.assertEquals(values2[i], RAW_DEPTH_2);
            }
            if (i != 3) {
                Assert.assertEquals(values3[i], RAW_DEPTH_3);
            }
        }
    }

    @Test
    public void testNormalizeMatrix() {
        final DepthMatrix rawMatrix = makeTestMatrix("chr1", 1000, 5000, 200);
        final Map<String, Double> sampleMedians = makeSampleMedians();
        final DepthMatrix normalizedMatrix = DepthMatrixLoader.normalizeMatrix(rawMatrix, sampleMedians);
        Assert.assertEquals(normalizedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(normalizedMatrix.getBins(), rawMatrix.getBins());
        final double[] values1 = normalizedMatrix.getSample(SAMPLE_1);
        final double[] values2 = normalizedMatrix.getSample(SAMPLE_2);
        final double[] values3 = normalizedMatrix.getSample(SAMPLE_3);
        final int numBins = normalizedMatrix.getBins().size();
        Assert.assertEquals(values1.length, numBins);
        Assert.assertEquals(values2.length, numBins);
        Assert.assertEquals(values3.length, numBins);
        for (int i = 0; i < numBins; i++) {
            Assert.assertEquals(values1[i], RAW_DEPTH_1 / MEDIAN_1);
            Assert.assertEquals(values2[i], RAW_DEPTH_2 / MEDIAN_2);
            Assert.assertEquals(values3[i], RAW_DEPTH_3 / MEDIAN_3);
        }
    }

    @Test
    public void testNormalizeMatrixNoSamples() {
        final DepthMatrix rawMatrix = new DepthMatrix(Collections.singletonList(TEST_INTERVAL), Collections.emptyMap());
        final Map<String, Double> sampleMedians = makeSampleMedians();
        final DepthMatrix normalizedMatrix = DepthMatrixLoader.normalizeMatrix(rawMatrix, sampleMedians);
        Assert.assertEquals(normalizedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(normalizedMatrix.getBins(), rawMatrix.getBins());
    }

    @Test
    public void testNormalizeMatrixNoBins() {
        final DepthMatrix rawMatrix = new DepthMatrix(Collections.emptyList(), Collections.singletonMap(SAMPLE_1, new double[]{}));
        final Map<String, Double> sampleMedians = makeSampleMedians();
        final DepthMatrix normalizedMatrix = DepthMatrixLoader.normalizeMatrix(rawMatrix, sampleMedians);
        Assert.assertEquals(normalizedMatrix.getSampleSet(), rawMatrix.getSampleSet());
        Assert.assertEquals(normalizedMatrix.getBins(), rawMatrix.getBins());
        Assert.assertEquals(normalizedMatrix.getSample(SAMPLE_1).length, 0);
    }

}