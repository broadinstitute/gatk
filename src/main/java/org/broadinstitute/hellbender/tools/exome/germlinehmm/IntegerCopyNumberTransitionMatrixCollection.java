package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * This class holds a collection of {@link IntegerCopyNumberTransitionMatrixData} for different sex genotypes
 * and contigs.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionMatrixCollection {

    /**
     * The map {@link SexGenotypeContigPairKey} => integer copy number prior data
     */
    private final Map<SexGenotypeContigPairKey, IntegerCopyNumberTransitionMatrixData>
            sexGenotypeAndContigToTransitionMatrixMap;

    private int maxCopyNumberInCollection;

    /**
     *
     */
    public IntegerCopyNumberTransitionMatrixCollection() {
        sexGenotypeAndContigToTransitionMatrixMap = new HashMap<>();
        maxCopyNumberInCollection = -1;
    }

    /**
     *
     * @param sexGenotype
     * @param contig
     * @param priorData
     */
    public void add(@Nonnull final String sexGenotype,
                    @Nonnull final String contig,
                    @Nonnull final IntegerCopyNumberTransitionMatrixData priorData) {
        final SexGenotypeContigPairKey key = SexGenotypeContigPairKey.of(sexGenotype, contig);
        if (sexGenotypeAndContigToTransitionMatrixMap.keySet().contains(key)) {
            throw new IllegalArgumentException(String.format("The integer copy number transition matrix for sex" +
                    " genotype \"%s\" and contig \"%s\" is already defined. Check for duplicate entries in the" +
                    " table file.", sexGenotype, contig));
        }
        sexGenotypeAndContigToTransitionMatrixMap.put(key, priorData);
        maxCopyNumberInCollection = FastMath.max(maxCopyNumberInCollection, priorData.getMaxCopyNumber());
    }

    /**
     *
     * @param sexGenotype
     * @param contig
     * @return
     */
    public IntegerCopyNumberTransitionMatrixData get(@Nonnull final String sexGenotype,
                                                     @Nonnull final String contig) {
        return sexGenotypeAndContigToTransitionMatrixMap.get(checkKey(sexGenotype, contig));
    }

    /**
     * Returns a set of sex genotypes that are annotated at least once
     *
     * @return
     */
    public Set<String> getSexGenotypes() {
        return sexGenotypeAndContigToTransitionMatrixMap.keySet().stream()
                .map(SexGenotypeContigPairKey::getSexGenotype)
                .collect(Collectors.toSet());
    }

    /**
     * Returns a set of contigs that are annotated at least once
     *
     * @return
     */
    public Set<String> getContigs() {
        return sexGenotypeAndContigToTransitionMatrixMap.keySet().stream()
                .map(SexGenotypeContigPairKey::getContig)
                .collect(Collectors.toSet());
    }

    public static IntegerCopyNumberTransitionMatrixCollection read(final File tableFile) {
        return IntegerCopyNumberTransitionMatrixCollectionReader.read(tableFile);
    }

    public static IntegerCopyNumberTransitionMatrixCollection read(final Reader reader,
                                                                   final String parentPath) {
        return IntegerCopyNumberTransitionMatrixCollectionReader.read(reader, parentPath);
    }

    /**
     * Asserts that all of provided sex genotypes and contigs have assigned copy number transition
     * matrices.
     *
     * @param sexGenotypes a set of sex genotype identifier strings
     * @param contigs      a set of contig strings
     * @throws UserException.BadInput if a required data entry is not available
     */
    public void assertCompleteness(@Nonnull final Set<String> sexGenotypes,
                                   @Nonnull final Set<String> contigs) {
        for (final String sexGenotype : sexGenotypes) {
            for (final String contig : contigs) {
                checkKey(sexGenotype, contig);
            }
        }
    }

    /**
     * Asserts that all pairwise choices of (complete sex genotype, contig) are annotated
     */
    public void assertCompleteness() {
        assertCompleteness(getSexGenotypes(), getContigs());
    }

    public int getMaxCopyNumberInCollection() {
        return maxCopyNumberInCollection;
    }

    /**
     *
     * @param sexGenotype
     * @param contig
     * @return
     */
    private SexGenotypeContigPairKey checkKey(@Nonnull final String sexGenotype,
                                              @Nonnull final String contig) {
        final SexGenotypeContigPairKey key = SexGenotypeContigPairKey.of(sexGenotype, contig);
        if (!sexGenotypeAndContigToTransitionMatrixMap.containsKey(key)) {
            throw new UserException.BadInput(String.format("The collection does not include" +
                            " the integer copy number transition matrix for sex genotype \"%s\" and contig \"%s\"",
                    sexGenotype, contig));
        }
        return key;
    }

    /**
     * Finds the highest copy number for all sex genotypes and contigs, returns a clone with all transition
     * matrices padded to the highest copy number.
     *
     * @return a new instance of {@link IntegerCopyNumberTransitionMatrixCollection}
     */
    public IntegerCopyNumberTransitionMatrixCollection cloneWithPaddedStates() {
        final IntegerCopyNumberTransitionMatrixCollection newCollection =
                new IntegerCopyNumberTransitionMatrixCollection();
        if (!sexGenotypeAndContigToTransitionMatrixMap.isEmpty()) {
            final int highestCopyNumber = sexGenotypeAndContigToTransitionMatrixMap.values().stream()
                    .mapToInt(IntegerCopyNumberTransitionMatrixData::getMaxCopyNumber)
                    .max().getAsInt();
            sexGenotypeAndContigToTransitionMatrixMap.keySet().stream()
                    .forEach(key -> {
                        final IntegerCopyNumberTransitionMatrixData dat = sexGenotypeAndContigToTransitionMatrixMap.get(key);
                        newCollection.add(key.getSexGenotype(), key.getContig(),
                                new IntegerCopyNumberTransitionMatrixData(dat.getTransitionMatrix(),
                                        highestCopyNumber - dat.getMaxCopyNumber()));
                    });
        }
        return newCollection;
    }
}