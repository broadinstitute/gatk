package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;
import java.util.HashMap;
import java.util.OptionalInt;
import java.util.stream.IntStream;

/**
 * We store a memo to avoid repeated computation of statistical power to detect a variant.
 * The key of the memo is a pair of numbers: number of reads and estimated allele fraction
 */
public class TumorPowerCalculator {
    private final double errorProbability;
    private final double tumorLODThreshold;
    private final double contamination;
    private final boolean enableSmoothing;

    private final HashMap<PowerCacheKey, Double> cache = new HashMap<>();

    public TumorPowerCalculator(double errorProbability, double constantLodThreshold, double contamination) {
        this(errorProbability, constantLodThreshold, contamination, true);
    }

    public TumorPowerCalculator(double errorProbability, double tumorLODThreshold, double contamination, boolean enableSmoothing) {
        this.errorProbability = errorProbability;
        this.tumorLODThreshold = tumorLODThreshold;
        this.contamination = contamination;
        this.enableSmoothing = enableSmoothing;
    }

    /**
     * A helper class that acts as the key to the memo of pre-computed power
     *
     * TODO: Not ideal to use double as a key. Refactor such that we use as keys numAlts and numReads, which are integers. Then calculate numAlts/numReads when we need allele fraction.
     *
     */
    private static class PowerCacheKey extends Pair<Integer, Double> {
        private final Double alleleFraction;
        private final Integer numReads;

        public PowerCacheKey(final int numReads, final double alleleFraction) {
            super(numReads, alleleFraction);
            this.alleleFraction = alleleFraction;
            this.numReads = numReads;
        }

        private boolean closeEnough(final double x, final double y, final double epsilon){
            return(Math.abs(x - y) < epsilon);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            PowerCacheKey that = (PowerCacheKey) o;
            return closeEnough(alleleFraction, that.alleleFraction, 0.001) && numReads != that.numReads;
        }

        @Override
        public int hashCode() {
            final long temp = alleleFraction != +0.0d ? Double.doubleToLongBits(alleleFraction) : 0L;
            return 31 * numReads + (int) (temp ^ (temp >>> 32));
        }
    }

    /**
     *
     * @param numReads                total number of reads, REF and ALT combined, in + or - strand
     * @param alleleFraction          the true allele fraction estimated as the combined allele fraction from + and - reads
     * @return                        probability of correctly calling the variant (i.e. power) given the above estimated allele fraction and number of reads.
     *                                we compute power separately for each strand (+ and -)
     */
    public double cachedPowerCalculation(final int numReads, final double alleleFraction) {
        final PowerCacheKey key = new PowerCacheKey(numReads, alleleFraction);
        // we first look up if power for given number of read and allele fraction has already been computed and stored in the cache.
        // if not we compute it and store it in teh cache.
        Double power = cache.get(key);
        if (power == null) {
            power = calculatePower(numReads, alleleFraction);
            cache.put(key, power);
        }
        return power;
    }

    /* helper function for calculateTumorLod */
    private double calculateLogLikelihood(final int numReads, final int numAlts, final double alleleFraction) {
        return((numReads-numAlts) * Math.log10( alleleFraction * errorProbability + (1 - alleleFraction)*(1 - errorProbability) ) +
                numAlts * Math.log10(alleleFraction * (1 - errorProbability) + (1 - alleleFraction) * errorProbability));
    }

    private double calculateTumorLod(final int numReads, final int numAlts) {
        final double alleleFraction = (double) numAlts / (double) numReads;
        final double altLikelihod = calculateLogLikelihood(numReads, numAlts, alleleFraction);
        final double refLikelihood = calculateLogLikelihood(numReads, numAlts, contamination);
        return altLikelihod - refLikelihood;
    }

    private double calculatePower(final int numReads, final double alleleFraction) {
        if (numReads==0) return 0;

        // TODO: add the factor of 1/3
        final double probAltRead = alleleFraction*(1 - errorProbability) + (1/3)*(1 - alleleFraction) * errorProbability;
        final BinomialDistribution binom = new BinomialDistribution(numReads, probAltRead);

        //TODO: use IndexRange::mapToDouble when latest gatk-public snapshot is used
        final double[] binomialProbabilities = IntStream.range(0, numReads + 1).mapToDouble(binom::probability).toArray();

        // find the smallest number of ALT reads k such that tumorLOD(k) > tumorLODThreshold
        final OptionalInt smallestKAboveLogThreshold = IntStream.range(0, numReads + 1)
                .filter(k -> calculateTumorLod(numReads, k) > tumorLODThreshold)
                .findFirst();

        if (! smallestKAboveLogThreshold.isPresent()){
            return 0;
        } else if (smallestKAboveLogThreshold.getAsInt() <= 0){
            throw new IllegalStateException("smallest k that meets the tumor LOD threshold is less than or equal to 0");
        }

        double power = Arrays.stream(binomialProbabilities, smallestKAboveLogThreshold.getAsInt(), binomialProbabilities.length).sum();

        // here we correct for the fact that the exact lod threshold is likely somewhere between
        // the k and k-1 bin, so we prorate the power from that bin
        if ( enableSmoothing ){
            final double tumorLODAtK = calculateTumorLod(numReads, smallestKAboveLogThreshold.getAsInt());
            final double tumorLODAtKMinusOne = calculateTumorLod(numReads, smallestKAboveLogThreshold.getAsInt()-1);
            final double weight = 1 - (tumorLODThreshold - tumorLODAtKMinusOne ) / (tumorLODAtK - tumorLODAtKMinusOne);
            power += weight * binomialProbabilities[smallestKAboveLogThreshold.getAsInt() - 1];
        }

        return power;
    }
}