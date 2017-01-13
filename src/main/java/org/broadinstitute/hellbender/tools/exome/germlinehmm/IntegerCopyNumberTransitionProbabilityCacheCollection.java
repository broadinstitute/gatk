package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.Reader;
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

    private final boolean padded;
    private final int maxCopyNumber;
    private final Set<String> sexGenotypesSet;
    private final Set<String> contigsSet;

    /**
     * Public constructor from a {@link IntegerCopyNumberTransitionMatrixCollection}
     *
     * @param transitionMatrixCollection collection of transition matrices
     */
    public IntegerCopyNumberTransitionProbabilityCacheCollection(
            @Nonnull final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection,
            final boolean padStates) {
        sexGenotypeAndContigToTransitionMatrixCacheMap = new HashMap<>();
        constructCacheMapFromMatrixCollection(padStates ? transitionMatrixCollection.cloneWithPaddedStates()
                : transitionMatrixCollection);
        padded = padStates;
        if (padded) {
            maxCopyNumber = transitionMatrixCollection.getMaxCopyNumberInCollection();
        } else {
            maxCopyNumber = 0;
        }
        sexGenotypesSet = Collections.unmodifiableSet(transitionMatrixCollection.getSexGenotypes());
        contigsSet = Collections.unmodifiableSet(transitionMatrixCollection.getContigs());
    }

    /**
     * Public constructor from a transition matrix collection reader
     *
     * @param reader
     * @param parentPath
     * @param padStates
     */
    public IntegerCopyNumberTransitionProbabilityCacheCollection(@Nonnull final Reader reader,
                                                                 @Nonnull final String parentPath,
                                                                 final boolean padStates) {
        final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection =
                IntegerCopyNumberTransitionMatrixCollection.read(reader, parentPath);
        sexGenotypeAndContigToTransitionMatrixCacheMap = new HashMap<>();
        constructCacheMapFromMatrixCollection(padStates ? transitionMatrixCollection.cloneWithPaddedStates()
                : transitionMatrixCollection);
        padded = padStates;
        if (padded) {
            maxCopyNumber = transitionMatrixCollection.getMaxCopyNumberInCollection();
        } else {
            maxCopyNumber = 0;
        }
        sexGenotypesSet = Collections.unmodifiableSet(transitionMatrixCollection.getSexGenotypes());
        contigsSet = Collections.unmodifiableSet(transitionMatrixCollection.getContigs());
    }

    /**
     * Public constructor from a transition matrix collection file
     *
     * @param inputFile
     * @param padStates
     */
    public IntegerCopyNumberTransitionProbabilityCacheCollection(@Nonnull final File inputFile,
                                                                 final boolean padStates) {
        final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection =
                IntegerCopyNumberTransitionMatrixCollection.read(inputFile);
        sexGenotypeAndContigToTransitionMatrixCacheMap = new HashMap<>();
        constructCacheMapFromMatrixCollection(padStates ? transitionMatrixCollection.cloneWithPaddedStates()
                : transitionMatrixCollection);
        padded = padStates;
        if (padded) {
            maxCopyNumber = transitionMatrixCollection.getMaxCopyNumberInCollection();
        } else {
            maxCopyNumber = 0;
        }
        sexGenotypesSet = Collections.unmodifiableSet(transitionMatrixCollection.getSexGenotypes());
        contigsSet = Collections.unmodifiableSet(transitionMatrixCollection.getContigs());
    }

    private void constructCacheMapFromMatrixCollection(
            @Nonnull final IntegerCopyNumberTransitionMatrixCollection transitionMatrixCollection) {
        /* assert that the transition matrix collection is complete */
        transitionMatrixCollection.assertCompleteness();
        for (final String sexGenotype : transitionMatrixCollection.getSexGenotypes()) {
            for (final String contig : transitionMatrixCollection.getContigs()) {
                sexGenotypeAndContigToTransitionMatrixCacheMap.put(SexGenotypeContigPairKey.of(sexGenotype, contig),
                        new IntegerCopyNumberTransitionProbabilityCache(transitionMatrixCollection.get(sexGenotype, contig)));
            }
        }
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

    public void cacheTransitionMatrix(final int distance, final String sexGenotype, final String contig) {
        ParamUtils.isPositive(distance, "The distance between two two targets must be positive");
        sexGenotypeAndContigToTransitionMatrixCacheMap.get(checkKey(sexGenotype, contig))
                .cacheLogTransitionMatrix(distance);
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
        return contigsSet;
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
