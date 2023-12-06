package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.IntFunction;

/**
 * A class which holds a cache of GenotypingEngines for handling different ploidies
 * It is lazily initialized so GenotypingEngines are only created on as needed.
 * @param <T> the GenotypingEngine implementation this provides
 */
public final class MultiPloidyGenotyperCache<T extends GenotypingEngine> {
    private final Map<Integer, T> ploidyToGenotyperMap;
    private final IntFunction<T> ploidyToGenotyper;
    private final OverlapDetector<SimpleCount> ploidyRegions;
    private final int defaultPloidy;

    public MultiPloidyGenotyperCache(IntFunction<T> ploidyToGenotyper, int defaultPloidy, OverlapDetector<SimpleCount> ploidyRegions){
        this.defaultPloidy = defaultPloidy;
        this.ploidyRegions = ploidyRegions;
        this.ploidyToGenotyperMap = new LinkedHashMap<>();
        this.ploidyToGenotyper = ploidyToGenotyper;
    }

    /**
     * Determines the appropriate ploidy to use at given site for different genotyping engines.
     * @param region Current region of interest
     * @return if the region overlaps with one or more of the alternate ploidy regions, return the highest alternate ploidy
     * otherwise return the default ploidy
     */
   private int getPloidyToUseAtThisSite(Locatable region) {
       return ploidyRegions.getOverlaps(region)
               .stream()
               .mapToInt(SimpleCount::getCount)
               .max()
               .orElse(defaultPloidy);
   }

    /**
     *
     * @return the default genotyper that is configured for the default ploidy
     */
    public T getDefaultGenotypingEngine(){
        return ploidyToGenotyperMap.get(defaultPloidy);
    }

    /**
     * get an appropriate GenotypingEngine for the given region.
     * @param region the region to be genotyped
     * @return a GenotypingEngine configured for the highest ploidy in the region given
     */
    public T getGenotypingEngine(final Locatable region) {
        final int currentPloidy = getPloidyToUseAtThisSite(region);
        return ploidyToGenotyperMap.computeIfAbsent(currentPloidy, ploidyToGenotyper::apply);
    }
}
