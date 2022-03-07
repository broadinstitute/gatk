package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * This class holds caches of {@link GenotypeAlleleCounts} for multiple fixed ploidy, allele count pairs,
 * allowing for fast random access of genotypes.  Note that the increment method of GenotypeAlleleCounts is always fast,
 * so the caches here are only necessary when incremental traversal over genotypes in the canonical order is not possible.
 *
 *
 * This class is thread-safe since modifying the caches is synchronized.
 */
public final class GenotypesCache {

    /**
     * Maximum possible number of cached {@link GenotypeAlleleCounts} for each fixed ploidy
     */
    public static final int MAX_CACHE_SIZE = 5000;

    /**
     * Cache of GenotypeAlleleCounts objects by ploidy.  Format is caches[p][n] = nth genotype of ploidy p in canonical order,
     * with p up to the current maximum ploidy and n up to the maximum number of cached genotypes per table.
     */
    private static List<List<GenotypeAlleleCounts>> caches = new ArrayList<>();

    private GenotypesCache(){ }

    /**
     * Returns the GenotypeAlleleCounts associated to a particular ploidy and genotype index.
     *
     *  If the requested index is larger than {@link GenotypesCache#MAX_CACHE_SIZE},
     *  this method will construct the result iteratively from the largest cached object.  Thus if you are iterating
     *  through all genotype-allele-counts you should do sequentially using the iterator method to avoid a large efficiency drop.
     *
     * @param ploidy the ploidy
     * @param genotypeIndex  the genotype index in the canonical order
     * @return never {@code null}.
     */
    public static GenotypeAlleleCounts get(final int ploidy, final int genotypeIndex) {
        ensureCapacity(genotypeIndex, ploidy);
        Utils.validateArg(ploidy >= 0, "ploidy may not be negative");
        Utils.validateArg(genotypeIndex >= 0, "genotype index may not be negative");
        final List<GenotypeAlleleCounts> cache = caches.get(ploidy);
        if (genotypeIndex < cache.size()) {
            return cache.get(genotypeIndex);
        } else {
            final GenotypeAlleleCounts result = cache.get(cache.size() - 1).copy();
            result.increase(genotypeIndex + 1 - cache.size());
            return result;
        }
    }

    /**
     * Extends the genotype allele counts cache for a certain ploidy up to a given size
     *
     * This method is synchronized since it modifies the shared cache.
     */
    private static synchronized void extendCache(final int ploidy, final int newSize) {
        final List<GenotypeAlleleCounts> cache = caches.get(ploidy);

        if (cache.isEmpty()) {
            cache.add(GenotypeAlleleCounts.first(ploidy));
        }

        while (cache.size() < newSize) {
            cache.add(cache.get(cache.size() - 1).next());
        }
    }

    /**
     * Update cache if necessary
     */
    private static void ensureCapacity(final int genotypeIndex, final int ploidy) {
        // add empty lists of genotypes until we have initialized all ploidies up to and including this one
        while (ploidy >= caches.size()) {
            caches.add(new ArrayList<>());
        }

        final List<GenotypeAlleleCounts> cache = caches.get(ploidy);

        if (cache.size() <= genotypeIndex && cache.size() < MAX_CACHE_SIZE) {
            final int newSize = Math.min(Math.max(cache.size() * 2 + 1, genotypeIndex), MAX_CACHE_SIZE);
            extendCache(ploidy, newSize);
        }
    }
}