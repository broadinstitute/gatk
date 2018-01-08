package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simulates {@link CopyRatioSegmentedData} given parameter values for use in test classes.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioSimulatedData {
    private static final int MIN_INTERVALS_PER_SEGMENT = 3;
    private static final double LOG2_COPY_RATIO_MIN = CopyRatioModeller.LOG2_COPY_RATIO_MIN;
    private static final double LOG2_COPY_RATIO_MAX = CopyRatioModeller.LOG2_COPY_RATIO_MAX;
    private static final double SEGMENT_MEAN_MIN = -10.;
    private static final double SEGMENT_MEAN_MAX = 5.;

    private static PoissonDistribution makePoisson(final RandomGenerator rng, final double mean) {
        return new PoissonDistribution(rng, mean, PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
    }

    private final CopyRatioSegmentedData data;
    private final CopyRatioState trueState;

    CopyRatioSimulatedData(final SampleLocatableMetadata metadata,
                           final double variance,
                           final double outlierProbability,
                           final int numSegments,
                           final double averageIntervalsPerSegment,
                           final RandomGenerator rng) {
        final List<Double> segmentMeans = new ArrayList<>(numSegments);
        final List<Boolean> outlierIndicators = new ArrayList<>();
        final List<CopyRatio> copyRatios = new ArrayList<>();
        final List<SimpleInterval> segments = new ArrayList<>();

        final double standardDeviation = Math.sqrt(variance);

        final PoissonDistribution segmentLengthGenerator = makePoisson(rng, averageIntervalsPerSegment);
        final UniformRealDistribution segmentMeanGenerator = new UniformRealDistribution(rng, SEGMENT_MEAN_MIN, SEGMENT_MEAN_MAX);
        final UniformRealDistribution outlierGenerator = new UniformRealDistribution(rng, LOG2_COPY_RATIO_MIN, LOG2_COPY_RATIO_MAX);

        //put each segment on its own chromosome and sort in sequence-dictionary order
        final List<String> chromosomes = IntStream.range(0, numSegments)
                .mapToObj(i -> metadata.getSequenceDictionary().getSequence(i).getSequenceName())
                .collect(Collectors.toList());

        for (final String chromosome : chromosomes) {
            // calculate the range of interval indices for this segment
            final int numIntervalsInSegment = Math.max(MIN_INTERVALS_PER_SEGMENT, segmentLengthGenerator.sample());

            final double segmentMean = segmentMeanGenerator.sample();
            segmentMeans.add(segmentMean);

            //we will put all the intervals in this segment/chromosome at loci 1, 2, 3 etc
            segments.add(new SimpleInterval(chromosome, 1, numIntervalsInSegment));
            for (int interval = 1; interval < numIntervalsInSegment + 1; interval++) {

                final double log2CopyRatio;
                if (rng.nextDouble() < outlierProbability) {
                    outlierIndicators.add(true);
                    log2CopyRatio = outlierGenerator.sample();
                } else {
                    outlierIndicators.add(false);
                    log2CopyRatio = segmentMean + rng.nextGaussian() * standardDeviation;
                }

                copyRatios.add(new CopyRatio(new SimpleInterval(chromosome, interval, interval), log2CopyRatio));
            }
        }

        data = new CopyRatioSegmentedData(
                new CopyRatioCollection(metadata, copyRatios),
                new SimpleIntervalCollection(metadata, segments));
        trueState = new CopyRatioState(variance, outlierProbability, new CopyRatioState.SegmentMeans(segmentMeans), new CopyRatioState.OutlierIndicators(outlierIndicators));
    }

    CopyRatioSegmentedData getData() {
        return data;
    }

    CopyRatioCollection getCopyRatios() {
        return data.getCopyRatios();
    }

    CopyRatioState getTrueState() { return trueState; }
}
