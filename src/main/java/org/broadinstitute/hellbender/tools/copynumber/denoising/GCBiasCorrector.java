package org.broadinstitute.hellbender.tools.copynumber.denoising;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Learn multiplicative correction factors as a function of GC content using a simple regression.
 * Our regression curve is obtained by filling GC-content bins of width 0.01 with the coverages of
 * genomic intervals corresponding to each GC and taking the median to get a robust estimate of the curve.  
 * In order to smooth out bins with few data (i.e. extreme GC values that occur rarely) 
 * we then convolve these medians with an exponential kernel.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class GCBiasCorrector {
    private static final Logger logger = LogManager.getLogger(GCBiasCorrector.class);
    
    // GC bins are 0%, 1% . . . 100%
    private static final int NUMBER_OF_GC_BINS = 101;

    // scale (in units of GC content from 0 to 1) over which GC bias correlation decays
    // i.e. the bias at GC content = 0.3 and at 0.2 are correlated ~exp(-0.1 / correlationLength)
    // this is effectively a smoothness parameter used for regularizing the GC-bias estimate for
    // GC content bins that contain few genomic intervals)
    private static final double correlationLength = 0.02;
    private static final double correlationDecayRatePerBin = 1.0 / (correlationLength * NUMBER_OF_GC_BINS);

    // multiply by these to get a GC-bias correction as a function of GC
    private final double[] gcCorrectionFactors;

    // Apache Commons median doesn't work on empty arrays; this value is a placeholder to avoid exceptions
    private static final double DUMMY_VALUE_NEVER_USED = 1.0;

    /**
     * Learn multiplicative correction factors as a function of GC content.  Basically, learn a
     * regression curve of coverage vs. GC in order to divide by that curve later.
     *
     * @param readCounts           integer read counts for a single sample
     * @param intervalGCContent    GC content (assumed to be from 0.0 to 1.0) of genomic intervals for {@code readCounts}
     */
    private GCBiasCorrector(final RealVector readCounts,
                            final double[] intervalGCContent) {
        Utils.nonNull(readCounts);
        Utils.nonNull(intervalGCContent);
        Utils.validateArg(intervalGCContent.length > 0, "Number of intervals must be positive.");
        Utils.validateArg(intervalGCContent.length == readCounts.getDimension(),
                "Number of intervals in read-counts matrix and GC-content array do not match.");

        final List<List<Double>> readCountsByGC = new ArrayList<>(NUMBER_OF_GC_BINS);
        IntStream.range(0, NUMBER_OF_GC_BINS).forEach(n -> readCountsByGC.add(new ArrayList<>()));
        IntStream.range(0, intervalGCContent.length).forEach(n -> readCountsByGC.get(gcContentToBinIndex(intervalGCContent[n])).add(readCounts.getEntry(n)));
        gcCorrectionFactors = calculateCorrectionFactors(readCountsByGC);
    }

    /**
     * Corrects GC bias of read counts in place.
     * @param readCounts            samples x intervals, assumed to be integer counts
     * @param intervalGCContent     GC content (assumed to be from 0.0 to 1.0) of genomic intervals for {@code readCounts}
     */
    public static void correctGCBias(final RealMatrix readCounts,
                                     final double[] intervalGCContent) {
        Utils.nonNull(readCounts);
        Utils.nonNull(intervalGCContent);
        ParamUtils.isPositive(intervalGCContent.length, "Number of intervals must be positive.");
        Utils.validateArg(readCounts.getColumnDimension() == intervalGCContent.length,
                "Number of intervals in read-counts matrix and GC-content array do not match.");
        final double[] totalCoveragePerSample = IntStream.range(0, readCounts.getRowDimension())
                .mapToDouble(r -> readCounts.getRowVector(r).getL1Norm())
                .toArray();
        // each row (sample) has its own GC-bias curve, hence its own GC-bias corrector
        final List<GCBiasCorrector> gcBiasCorrectors = IntStream.range(0, readCounts.getRowDimension())
                .mapToObj(n -> new GCBiasCorrector(readCounts.getRowVector(n), intervalGCContent))
                .collect(Collectors.toList());

        // correct the input counts in-place
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int interval, double readCount) {
                return gcBiasCorrectors.get(row).correctedCoverage(readCount, intervalGCContent[interval]);
            }
        });

        // we would like the average correction factor to be 1.0 in the sense that average coverage before and after
        // correction should be equal
        final double[] sampleNormalizationFactors = IntStream.range(0, readCounts.getRowDimension())
                .mapToDouble(r -> totalCoveragePerSample[r] / readCounts.getRowVector(r).getL1Norm())
                .toArray();
        readCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int interval, double coverage) {
                return coverage * sampleNormalizationFactors[row];
            }
        });
    }

    /**
     * As described above, calculate medians of each GC bin and convolve with an exponential kernel.
     *
     * @param readCountsByGC integer read counts for each GC bin
     * @return multiplicative correction factors for each GC bin
     */
    private double[] calculateCorrectionFactors(final List<List<Double>> readCountsByGC) {
        final RealVector medians = new ArrayRealVector(readCountsByGC.stream()
                .mapToDouble(GCBiasCorrector::medianOrDefault)
                .toArray());
        return IntStream.range(0, NUMBER_OF_GC_BINS)
                .mapToDouble(bin -> {
                    final RealVector weights = new ArrayRealVector(IntStream.range(0, NUMBER_OF_GC_BINS)
                            .mapToDouble(n -> readCountsByGC.get(n).size() * Math.exp(-Math.abs(bin - n) * correlationDecayRatePerBin))
                            .toArray());
                    return weights.dotProduct(medians) / weights.getL1Norm();})
                .map(x -> 1 / x)
                .toArray();
    }

    private double correctedCoverage(final double readCount,
                                     final double gcContent) {
        return gcCorrectionFactors[gcContentToBinIndex(gcContent)] * readCount;
    }

    // return a median of coverages or dummy default value if no coverage exists at this GC bin
    // this default is never used because empty bins get zero weight in {@code calculateCorrectionFactors}
    private static double medianOrDefault(final List<Double> list) {
        return list.size() > 0 ? new Median().evaluate(list.stream().mapToDouble(d -> d).toArray()) : DUMMY_VALUE_NEVER_USED;
    }

    private static int gcContentToBinIndex(final double gcContent) {
        return (int) Math.round(gcContent * (NUMBER_OF_GC_BINS - 1));
    }
}
