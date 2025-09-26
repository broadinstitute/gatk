package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.DepthMatrix;
import org.broadinstitute.hellbender.tools.sv.DepthMatrixTest;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class DepthEvidenceTestTest {

    private static final String SAMPLE_PREFIX = "sample";
    private static final double DEFAULT_COVERAGE = 1.0;
    private static final double DEFAULT_VAR_OFFSET = 0.1;
    private static final double DEFAULT_AMPLITUDE = 0.1;
    private static final double DEFAULT_PERIOD = 0.001;
    private static final double FLOAT_TOLERANCE = 1e-8;

    private static List<String> generateSamples(final int n) {
        final List<String> samples = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            samples.add(SAMPLE_PREFIX + i);
        }
        return samples;
    }

    private static DepthMatrix generateMatrix(final String contig, final int start, final int end, final int binSize,
            final List<String> samples, final Set<String> carriers, final double carrierOffset, final double jitter) {
        final List<SimpleInterval> intervals = DepthMatrixTest.makeBins(contig, start, end, binSize);
        final Map<String, double[]> counts = new HashMap<>();
        int i = 0;
        for (final String sample : samples) {
            final double[] array = new double[intervals.size()];
            final double sampleFrac = (++i / (double) samples.size());
            double cov = DEFAULT_COVERAGE + ((sampleFrac - 0.5) * jitter);
            if (carriers.contains(sample)) {
                cov += carrierOffset;
            }
            for (int j = 0; j < array.length; j++) {
                array[j] = cov + DEFAULT_AMPLITUDE * Math.sin(Math.PI * DEFAULT_PERIOD * (start + binSize * j) + sampleFrac);
            }
            counts.put(sample, array);
        }
        return new DepthMatrix(intervals, counts);
    }

    @DataProvider(name = "testData")
    public Object[][] testData() {
        return new Object[][]{
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 0, 0, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Double.NaN, Double.NaN, Double.NaN, true},
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 1, 0, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Double.NaN, Double.NaN, Double.NaN, true},
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 1, 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Double.NaN, Double.NaN, Double.NaN, false},
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 2, 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Double.NaN, Double.NaN, 0.25, false},
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 30, 5, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 5.60222062839788E-6, 1.2407337293796061E-5, 0.25000000000000044, false},
                {"chr1", 10001, 20001, 1000, 0.3, 0.8, 3, 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 2.034760087224475E-4, 3.934273404714304E-4, 0.2500000000000002, false},
                {"chr1", 10001, 20001, 1000, 0.3, 0.8, 3, 2, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, Double.NaN, Double.NaN, 0.2500000000000001, false},
                {"chr1", 10001, 20001, 1000, 0.4, 0.999, 10, 2, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 1.7062657260868974E-5, 3.664593219676604E-5, 0.30000000000000027, false},
                {"chr1", 10001, 20001, 1000, 0.4, 0.8, 30, 5, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 1.170372187608848E-5, 1.6292180690546942E-5, 0.30000000000000027, false},
                {"chr1", 11001, 21001, 1000, 0.5, 0.8, 30, 5, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 1.9390832026622284E-5, 2.423101722470733E-5, 0.34999999999999987, false},
                {"chr1", 100001, 200001, 1000, 0.6, 0.8, 30, 5, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 2.7939964670098405E-5, 3.26516698243573E-5, 0.40000000000000013, false},
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 30, 25, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 5.60222062839788E-6, 1.2380024605906925E-5, 0.25, false},
                {"chr1", 10001, 20001, 100, 0.3, 0.8, 30, 5, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, 0.9280674441568251, 0.982075541216225, -0.05000000000000049, false},
                {"chr1", 10001, 20001, 100, 0.3, 0.999, 5, 3, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, 0.9138783402445532, 0.9862873851507764, -0.050000000000000266, false},
                {"chr1", 10001, 20001, 100, 0.1, 0.8, 30, 5, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, 2.5254565399279727E-4, 0.4738946487731285, 0.0499999999999996, false},
        };
    }

    @Test(dataProvider = "testData")
    public void test(final String contig, final int start, final int end, final int binSize, final double jitter, final double powerThreshold, final int numSamples, final int numCarriers, final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                     final double expectedP, final double expectedSecondMaxP, final double expectedMedianSeparation, final boolean expectNull) {
        final DepthEvidenceTest test = new DepthEvidenceTest(powerThreshold);
        final List<String> samples = generateSamples(numSamples);
        final Set<String> carriers = new HashSet<>(samples.subList(0, numCarriers));
        final SVCallRecord record;
        final double coverageOffset;
        if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            record = SVTestUtils.makeDeletionRecordWithCoordsAndCarriers(contig, start, end, samples, carriers);
            coverageOffset = -DEFAULT_VAR_OFFSET;
        } else if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            record = SVTestUtils.makeDuplicationRecordWithCoordsAndCarriers(contig, start, end, samples, carriers);
            coverageOffset = DEFAULT_VAR_OFFSET;
        } else {
            throw new TestException("Unsupported SV type " + svtype);
        }
        final DepthMatrix matrix = generateMatrix(contig, start, end, binSize, samples, carriers, coverageOffset, jitter);
        final DepthEvidenceTest.DepthTestResult result = test.test(record, matrix);
        if (expectNull) {
            Assert.assertNull(result);
        } else {
            SVTestUtils.assertFloatWithinTolerance(result.pValue(), expectedP, FLOAT_TOLERANCE);
            SVTestUtils.assertFloatWithinTolerance(result.secondMaxP(), expectedSecondMaxP, FLOAT_TOLERANCE);
            SVTestUtils.assertFloatWithinTolerance(result.medianSeparation(), expectedMedianSeparation, FLOAT_TOLERANCE);
        }
    }

}