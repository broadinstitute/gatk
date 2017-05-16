package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * This class constructs a collection of integer copy number transition probability caches, allowing
 * streamlined query of the transition matrix for each contig, sex genotype, and distance.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionProbabilityCacheCollection implements Serializable {

    private static final long serialVersionUID = -7708106959687292204L;

    private final Map<SexGenotypeContigPairKey, IntegerCopyNumberTransitionProbabilityCache>
            sexGenotypeAndContigToTransitionMatrixCacheMap;

    /**
     * This is an arbitrary integer and is used as a placeholder for {@link #maxCopyNumber} if the transition
     * matrix collection is not padded; see {@link #getMaxCopyNumber()}.
     */
    private static final int UNDEFINED_MAX_COPY_NUMBER = -1;

    private final boolean padded;
    private final int maxCopyNumber;
    private final Set<String> sexGenotypesSet;
    private final Set<String> contigSet;

    /**
     * Public constructor from a {@link IntegerCopyNumberTransitionMatrixCollection}
     *
     * @param transitionMatrixCollection collection of transition matrices
     * @param padStates whether or not pad all transition matrices to the maximum copy number in the collection
     */
    public IntegerCopyNumberTransitionProbabilityCacheCollection(
            @Nonnull final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection,
            final boolean padStates) {
        Utils.nonNull(transitionMatrixCollection, "The transition matrix cache collection must be non-null");
        sexGenotypeAndContigToTransitionMatrixCacheMap = constructCacheMapFromMatrixCollection(padStates
                ? transitionMatrixCollection.cloneWithPaddedStates()
                : transitionMatrixCollection);
        padded = padStates;
        maxCopyNumber = padded
                ? transitionMatrixCollection.getMaxCopyNumberInCollection()
                : UNDEFINED_MAX_COPY_NUMBER;
        sexGenotypesSet = Collections.unmodifiableSet(transitionMatrixCollection.getSexGenotypes());
        contigSet = Collections.unmodifiableSet(transitionMatrixCollection.getContigs());
    }

    private static Map<SexGenotypeContigPairKey, IntegerCopyNumberTransitionProbabilityCache> constructCacheMapFromMatrixCollection(
            @Nonnull final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection) {
        /* assert that the transition matrix collection is complete */
        transitionMatrixCollection.assertCompleteness();
        final Map<SexGenotypeContigPairKey, IntegerCopyNumberTransitionProbabilityCache> cacheMap = new HashMap<>();
        for (final String sexGenotype : transitionMatrixCollection.getSexGenotypes()) {
            for (final String contig : transitionMatrixCollection.getContigs()) {
                cacheMap.put(SexGenotypeContigPairKey.of(sexGenotype, contig),
                        new IntegerCopyNumberTransitionProbabilityCache(transitionMatrixCollection.get(sexGenotype, contig)));
            }
        }
        return cacheMap;
    }

    /**
     * Calculates the log transition probability between two integer copy number states for a given
     * distance, sex genotype, and contig
     *
     * @param distance distance between the two targets (must be positive)
     * @param sexGenotype a string identifier for the sex genotype
     * @param contig a string identifier for the contig
     * @param to the destination integer copy number state
     * @param from the departure integer copy number state
     * @return a double value
     */
    public double logTransitionProbability(final int distance, final String sexGenotype, final String contig,
                                           final IntegerCopyNumberState to, final IntegerCopyNumberState from) {
        ParamUtils.isPositive(distance, "The distance between two two targets must be positive");
        Utils.nonNull(to, "The destination state must be non-null");
        Utils.nonNull(from, "The departure state must be non-null");
        return sexGenotypeAndContigToTransitionMatrixCacheMap.get(checkKey(sexGenotype, contig))
                .logTransitionProbability(distance, to, from);
    }

    public double logStationaryProbability(final String sexGenotype, final String contig,
                                           final IntegerCopyNumberState state) {
        return sexGenotypeAndContigToTransitionMatrixCacheMap.get(checkKey(sexGenotype, contig))
                .logStationaryProbability(state);
    }

    public void clearCaches() {
        sexGenotypeAndContigToTransitionMatrixCacheMap.values()
                .forEach(IntegerCopyNumberTransitionProbabilityCache::clearCache);
    }

    public int getMaxCopyNumber(final String sexGenotype, final String contig) {
        return sexGenotypeAndContigToTransitionMatrixCacheMap.get(checkKey(sexGenotype, contig))
                .getMaxCopyNumber();
    }

    public int getMaxCopyNumber() {
        if (!padded) {
            throw new IllegalStateException("The global maximum copy number is only well-defined if the collection" +
                    " is padded with extra states upon construction");
        } else {
            return maxCopyNumber;
        }
    }

    public boolean isPadded() {
        return padded;
    }

    public Set<String> getSexGenotypes() {
        return sexGenotypesSet;
    }

    public Set<String> getContigs() {
        return contigSet;
    }

    private SexGenotypeContigPairKey checkKey(final String sexGenotype, final String contig) {
        final SexGenotypeContigPairKey key = SexGenotypeContigPairKey.of(sexGenotype, contig);
        if (!sexGenotypeAndContigToTransitionMatrixCacheMap.containsKey(key)) {
            throw new IllegalArgumentException(String.format("Can not calculate the transition matrix" +
                    " for sex genotype \"%s\" and contig \"%s\"", sexGenotype, contig));
        }
        return key;
    }
}
