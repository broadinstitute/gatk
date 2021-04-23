package org.broadinstitute.hellbender.tools.walkers.coverage;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * A class helper for storing running intervalPartition data.
 *
 * The primary purpose of this class is to associate various {@link DepthOfCoverageStats} objects with eachother and make
 * sure that they all get mutated together per locus.
 *
 * Generally speaking, one of these objects should be constructed for each interval or aggregation of intervals that
 * needs to be reported on as a discrete unit together. initialize() should be called so that the samples for each of
 * the underlying DepthOfCoverageStats objects can be initialized with the samples for its traversal and update() should
 * be called on this object once for each locus this object corresponds to.
 */
public class DepthOfCoveragePartitionedDataStore {
    private Map<DoCOutputType.Partition,DepthOfCoverageStats> coverageProfiles;

    /**
     * @param typesToUse Partitions by which to subdivide the coverage counts
     * @param start histogram bin starting coverage
     * @param stop histogram bin ending coverage
     * @param nBins number of bins to create between start and stop, these are used to generate a histogram.
     */
    public DepthOfCoveragePartitionedDataStore(final Collection<DoCOutputType.Partition> typesToUse, final int start, final int stop, final int nBins,
                                               final boolean includeDeletions, final boolean omitLocusTable, final Map<DoCOutputType.Partition, List<String>> globalIdentifierMap) {
        coverageProfiles = new TreeMap<>();
        for ( DoCOutputType.Partition type : typesToUse ) {
            coverageProfiles.put(type,new DepthOfCoverageStats(CoverageUtils.calculateCoverageHistogramBinEndpoints(start,stop,nBins), includeDeletions, omitLocusTable));
        }
        initializeStats(omitLocusTable, globalIdentifierMap);
    }

    // Calls out to the underlying DepthOfCoverageStats objects to ensure that they had all of their initilization steps called on them
    public void initializeStats(final boolean omitLocusTable, final Map<DoCOutputType.Partition, List<String>> globalIdentifierMap) {
        for ( DoCOutputType.Partition t : coverageProfiles.keySet() ) {
            // Make sure the identifiers have been added to the underlying DepthOfCoverageStats object
            for ( String sample : globalIdentifierMap.get(t) ) {
                coverageProfiles.get(t).initializeSample(sample);
            }
            if ( ! omitLocusTable ) {
                coverageProfiles.get(t).initializeLocusCounts();
            }
        }
    }

    // Adds each of the DepthOfCoverageStats objects information summarizing a new locus
    public void addLocusData(final Map<DoCOutputType.Partition, Map<String,int[]>> countsByIdentifierByType) {
        for ( DoCOutputType.Partition t : coverageProfiles.keySet() ) {
            coverageProfiles.get(t).update(countsByIdentifierByType.get(t));
        }
    }

    // Returns the underlying DepthOfCoverageStats object for a given partition
    public DepthOfCoverageStats getCoverageByAggregationType( final DoCOutputType.Partition t) {
        return coverageProfiles.get(t);
    }
}
