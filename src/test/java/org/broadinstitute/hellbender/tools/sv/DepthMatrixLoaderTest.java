package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.*;

public class DepthMatrixLoaderTest extends GATKBaseTest {

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