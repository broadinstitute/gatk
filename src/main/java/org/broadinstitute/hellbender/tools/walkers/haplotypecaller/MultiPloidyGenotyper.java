package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.function.IntFunction;

public class MultiPloidyGenotyper<T extends GenotypingEngine> {
    public static final int DEFAULT_PLOIDY_SENTINEL = -1;
    private final Map<Integer, T> ploidyToGenotyperMap;
    private final OverlapDetector<SimpleCount> ploidyRegions;
    private final int defaultPloidy;

    public MultiPloidyGenotyper(IntFunction<T> ploidyToGenotyper, Set<Integer> allCustomPloidies, int defaultPloidy, OverlapDetector<SimpleCount> ploidyRegions){
        this.defaultPloidy = defaultPloidy;
        this.ploidyRegions = ploidyRegions;

        ploidyToGenotyperMap = new LinkedHashMap<>();
        ploidyToGenotyperMap.put(DEFAULT_PLOIDY_SENTINEL, ploidyToGenotyper.apply(defaultPloidy));
        // Create other custom genotyping engines if user provided ploidyRegions
        for (final int ploidy : allCustomPloidies) {
            this.ploidyToGenotyperMap.put(ploidy, ploidyToGenotyper.apply(ploidy));
        }
    }

    /**
     * Determines the appropriate ploidy to use at given site for different genotyping engines.
     * @param region Current region of interest
     * @return Ploidy value to use here given user inputs, or -1 if fall back to default
     */
   private int getPloidyToUseAtThisSite(Locatable region) {
       return ploidyRegions.getOverlaps(region)
               .stream()
               .mapToInt(SimpleCount::getCount)
               .max()
               .orElse(DEFAULT_PLOIDY_SENTINEL);
   }

    public T getDefaultGenotypingEngine(){
        return ploidyToGenotyperMap.get(DEFAULT_PLOIDY_SENTINEL);
    }

    public T getGenotypingEngine(final Locatable region) {
        final int currentPloidy = getPloidyToUseAtThisSite(region);
        return ploidyToGenotyperMap.get(currentPloidy);
    }
}
