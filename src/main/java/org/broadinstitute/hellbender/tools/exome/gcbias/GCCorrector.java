package org.broadinstitute.hellbender.tools.exome.gcbias;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Learn multiplicative correction factors as a function of GC content from coverage vs. GC data.  Basically, learn a
 * regression curve of coverage vs. GC, and divide by that curve to get GC-corrected coverage.
 *
 * Our regression curve is obtained by filling GC content bins of width 0.01 with the coverages of targets corresponding to each GC
 * and taking the median to get a robust estimate of the curve.  In order to smooth out bins with few data (i.e. extreme
 * GC values that occur rarely) we then convolve these medians with an exponential kernel.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
class GCCorrector {
    // GC bins are 0%, 1% . . . 100%
    private final static int NUMBER_OF_GC_BINS = 101;

    // scale (in units of GCcontent from 0 to 1) over which gc bias correlation decays
    // i.e. the bias at GC content = 0.3 and at 0.2 are correlated ~exp(-0.1/correlationLength)
    // this is effectively a smoothness parameter used for regularizing the GC bias estimate for
    // GC content bins that have few targets)
    private static final double correlationLength = 0.02;
    private static final double correlationDecayRatePerBin = 1.0 / (correlationLength * NUMBER_OF_GC_BINS);

    // multiply by these to get a GC correction as a function of GC
    private final double[] gcCorrectionFactors;

    //Apache commons median doesn't work on empty arrays; this value is a placeholder to avoid exceptions
    private static final double DUMMY_VALUE_NEVER_USED = 1.0;

    /**
     * Learn multiplicative correction factors as a function of GC from coverage vs. GC data.  Basically, learn a
     * regression curve of coverage vs. GC in order to divide by that curve later.
     *
     * @param gcContents GC content (from 0.0 to 1.0) of targets in {@code coverage}
     * @param coverage raw of proportional coverage
     */
    public GCCorrector(final double[] gcContents, final RealVector coverage) {
        Utils.nonNull(gcContents);
        Utils.nonNull(coverage);
        Utils.validateArg(gcContents.length > 0, "must have at lest one datum");
        Utils.validateArg(gcContents.length == coverage.getDimension(), "must have one gc value per coverage.");

        final List<List<Double>> coveragesByGC = new ArrayList<>(NUMBER_OF_GC_BINS);
        IntStream.range(0, NUMBER_OF_GC_BINS).forEach(n -> coveragesByGC.add(new ArrayList<>()));
        IntStream.range(0, gcContents.length).forEach(n -> coveragesByGC.get(gcContentToBinIndex(gcContents[n])).add(coverage.getEntry(n)));
        gcCorrectionFactors = calculateCorrectionFactors(coveragesByGC);
    }

    /**
     * As described above, calculate medians of each GC bin and convolve with an exponential kernel.
     *
     * @param coveragesByGC list of coverages for each GC bin
     * @return multiplicative correction factors for each GC bin
     */
    private double[] calculateCorrectionFactors(final List<List<Double>> coveragesByGC) {
        final RealVector medians = new ArrayRealVector(coveragesByGC.stream().mapToDouble(GCCorrector::medianOrDefault).toArray());
        return IntStream.range(0, NUMBER_OF_GC_BINS).mapToDouble(bin -> {
            final RealVector weights = new ArrayRealVector(IntStream.range(0, NUMBER_OF_GC_BINS)
                    .mapToDouble(n -> coveragesByGC.get(n).size() * Math.exp(-Math.abs(bin - n) * correlationDecayRatePerBin)).toArray());
            return weights.dotProduct(medians) / weights.getL1Norm();
        }).map(x -> 1/x).toArray();
    }

    /**
     *
     * @param inputCounts raw coverage before GC correction
     * @param gcContentByTarget      array of gc contents, one per target of the input
     * @return              GC-corrected coverage
     */
    public static ReadCountCollection correctCoverage(final ReadCountCollection inputCounts, final double[] gcContentByTarget) {
        // each column (sample) has its own GC bias curve, hence its own GC corrector
        final List<GCCorrector> gcCorrectors = IntStream.range(0, inputCounts.columnNames().size())
                .mapToObj(n -> new GCCorrector(gcContentByTarget, inputCounts.counts().getColumnVector(n))).collect(Collectors.toList());

        // gc correct a copy of the input counts in-place
        final RealMatrix correctedCounts = inputCounts.counts().copy();
        correctedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int target, int column, double coverage) {
                return gcCorrectors.get(column).correctedCoverage(coverage, gcContentByTarget[target]);
            }
        });

        // we would like the average correction factor to be 1.0 in the sense that average coverage before and after
        // correction should be equal
        final double[] columnNormalizationFactors = IntStream.range(0, inputCounts.columnNames().size())
                .mapToDouble(c -> inputCounts.counts().getColumnVector(c).getL1Norm() / correctedCounts.getColumnVector(c).getL1Norm()).toArray();
        correctedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int target, int column, double coverage) {
                return coverage * columnNormalizationFactors[column];
            }
        });

        return new ReadCountCollection(inputCounts.targets(), inputCounts.columnNames(), correctedCounts);
    }

    private double correctedCoverage(final double coverage, final double gcContent) {
        return gcCorrectionFactors[gcContentToBinIndex(gcContent)] * coverage;
    }

    // return a median of coverages or dummy default value if no coverage exists at this gc bin
    // this default is never used because empty bins get zero weight in {@code calculateCorrectionFactors}
    private static double medianOrDefault(final List<Double> list) {
        return list.size() > 0 ? new Median().evaluate(list.stream().mapToDouble(d->d).toArray()) : DUMMY_VALUE_NEVER_USED;
    }

    private static int gcContentToBinIndex(final double gcContent) {
        return (int) Math.round(gcContent * (NUMBER_OF_GC_BINS-1));
    }
}
