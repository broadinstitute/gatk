package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simulates {@link AlleleFractionSegmentedData} given parameter values for use in test classes.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionSimulatedData {
    private static final int MIN_HETS_PER_SEGMENT = 3;

    private static PoissonDistribution makePoisson(final RandomGenerator rng, final double mean) {
        return new PoissonDistribution(rng, mean, PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
    }

    private final AlleleFractionSegmentedData data;
    private final AlleleFractionState trueState;

    AlleleFractionSimulatedData(final SampleLocatableMetadata metadata,
                                final AlleleFractionGlobalParameters globalParameters,
                                final int numSegments,
                                final double averageHetsPerSegment,
                                final double averageDepth,
                                final RandomGenerator rng) {
        final AlleleFractionState.MinorFractions minorFractions = new AlleleFractionState.MinorFractions(numSegments);
        final List<AllelicCount> allelicCounts = new ArrayList<>();
        final List<SimpleInterval> segments = new ArrayList<>();

        final PoissonDistribution segmentLengthGenerator = makePoisson(rng, averageHetsPerSegment);
        final PoissonDistribution readDepthGenerator = makePoisson(rng, averageDepth);
        final UniformRealDistribution minorFractionGenerator = new UniformRealDistribution(rng, 0.0, 0.5);

        final double meanBias = globalParameters.getMeanBias();
        final double biasVariance = globalParameters.getBiasVariance();
        final double outlierProbability = globalParameters.getOutlierProbability();

        //translate to ApacheCommons' parametrization of the gamma distribution
        final double gammaShape = meanBias * meanBias / biasVariance;
        final double gammaScale = biasVariance / meanBias;
        final GammaDistribution biasGenerator = new GammaDistribution(rng, gammaShape, gammaScale);

        //put each segment on its own chromosome and sort in sequence-dictionary order
        final List<String> chromosomes = IntStream.range(0, numSegments)
                .mapToObj(i -> metadata.getSequenceDictionary().getSequence(i).getSequenceName())
                .collect(Collectors.toList());

        for (final String chromosome : chromosomes) {
            // calculate the range of het indices for this segment
            final int numHetsInSegment = Math.max(MIN_HETS_PER_SEGMENT, segmentLengthGenerator.sample());

            final double minorFraction = minorFractionGenerator.sample();
            minorFractions.add(minorFraction);

            //we will put all the hets in this segment/chromosome at loci 1, 2, 3 etc
            segments.add(new SimpleInterval(chromosome, 1, numHetsInSegment));
            for (int het = 1; het < numHetsInSegment + 1; het++) {
                final double bias = biasGenerator.sample();

                //flip a coin to decide alt minor (alt fraction = minor fraction) or ref minor (alt fraction = 1 - minor fraction)
                final boolean isAltMinor = rng.nextDouble() < 0.5;
                final double altFraction =  isAltMinor ? minorFraction : 1 - minorFraction;

                //the probability of an alt read is the alt fraction modified by the bias or, in the case of an outlier, random
                final double pAlt;
                if (rng.nextDouble() < outlierProbability) {
                    pAlt = rng.nextDouble();
                } else {
                    pAlt = altFraction / (altFraction + (1 - altFraction) * bias);
                }

                final int numReads = readDepthGenerator.sample();
                final int numAltReads = new BinomialDistribution(rng, numReads, pAlt).sample();
                final int numRefReads = numReads - numAltReads;
                allelicCounts.add(new AllelicCount(new SimpleInterval(chromosome, het, het), numRefReads, numAltReads));
            }
        }

        data = new AlleleFractionSegmentedData(
                new AllelicCountCollection(metadata, allelicCounts),
                new SimpleIntervalCollection(metadata, segments));
        trueState = new AlleleFractionState(meanBias, biasVariance, outlierProbability, minorFractions);
    }

    AlleleFractionSegmentedData getData() {
        return data;
    }

    AllelicCountCollection getAllelicCounts() {
        return data.getAllelicCounts();
    }

    AlleleFractionState getTrueState() { return trueState; }
}
