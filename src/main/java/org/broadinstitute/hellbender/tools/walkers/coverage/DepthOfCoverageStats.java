package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.utils.BaseUtils;

import java.util.*;

/**
 * A class for storing summarized coverage statistics for DepthOfCoverage.
 */
public class DepthOfCoverageStats {

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD DATA
    ////////////////////////////////////////////////////////////////////////////////////

    private Map<String,long[]> granularHistogramBySample; // holds the counts per each bin
    private Map<String,Long> totalCoverages; // holds total coverage per sample
    private int[] binLeftEndpoints; // describes the left endpoint for each bin
    private long[][] locusCoverageCounts; // holds counts of number of bases with >=X samples at >=Y coverage
    private boolean tabulateLocusCounts = false;
    private long nLoci; // number of loci seen
    private long totalDepthOfCoverage;
    private boolean includeDeletions;

    ////////////////////////////////////////////////////////////////////////////////////
    // TEMPORARY DATA ( not worth re-instantiating )
    ////////////////////////////////////////////////////////////////////////////////////

    private int[] locusHistogram; // holds a histogram for each locus; reset after each update() call
    private int totalLocusDepth; // holds the total depth of coverage for each locus; reset after each update() call

    ////////////////////////////////////////////////////////////////////////////////////
    // INITIALIZATION METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    /**
     *
     * @param leftEndpoints Array of the left endpoints for desired coverage bins to track. See {@link}
     * @param includeDeletions whether to count deletions at a site as coverage in the total
     * @param dontComputeLocusTable if true then skip calculations for per-locus output tables
     */
    public DepthOfCoverageStats(int[] leftEndpoints, boolean includeDeletions, boolean dontComputeLocusTable) {
        this.binLeftEndpoints = leftEndpoints;
        this.granularHistogramBySample = new HashMap<>();
        this.totalCoverages = new HashMap<>();
        this.nLoci = 0;
        this.totalLocusDepth = 0;
        this.totalDepthOfCoverage = 0;
        this.includeDeletions = includeDeletions;
        if ( ! dontComputeLocusTable ) {
            initializeLocusCounts();
        }
    }

    // Adds to granularHistogramBySample a per-locus histogram for the given sample to track
    void initializeSample(String sample) {
        if ( granularHistogramBySample.containsKey(sample) ) {
            return;
        }
        // Note: we add 1 here because bin left endpoints has values for ever bin except for the first one (which covers everything below the defined lower bound)
        granularHistogramBySample.put(sample, new long[this.binLeftEndpoints.length+1]);
        totalCoverages.put(sample,0L);
    }

    // Create the per-locus coverage counting arrays and histogram
    void initializeLocusCounts() {
        // Note: we add 1 here because bin left endpoints has values for ever bin except for the first one (which covers everything below the defined lower bound)
        locusCoverageCounts = new long[granularHistogramBySample.size()][binLeftEndpoints.length+1];
        locusHistogram = new int[binLeftEndpoints.length+1];
        tabulateLocusCounts = true;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // UPDATE METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    private void updateDepths(Map<String, Integer> depthBySample) {
        int coverageThresholdIndex;
        for ( String sample : granularHistogramBySample.keySet() ) {
            if ( depthBySample.containsKey(sample) ) {
                coverageThresholdIndex = updateSample(sample,depthBySample.get(sample));
                totalLocusDepth += depthBySample.get(sample);
            } else {
                coverageThresholdIndex = updateSample(sample,0);
            }

            if ( tabulateLocusCounts ) {
                for ( int i = 0; i <= coverageThresholdIndex; i ++ ) {
                    locusHistogram[i]++;
                }
            }
        }
        updateLocusCounts(locusHistogram);

        nLoci++;
        totalDepthOfCoverage += totalLocusDepth;
        totalLocusDepth = 0;
    }

    public void update(Map<String,int[]> countsBySample) {
        if ( countsBySample == null ) {
            this.updateDepths(new HashMap<>(1));
            return;
        }
        final HashMap<String,Integer> depthBySample = new HashMap<String,Integer>();
        for ( Map.Entry<String, int[]> entry : countsBySample.entrySet() ) {
            int total = 0;
            int[] counts = entry.getValue();
            for ( byte base : BaseUtils.BASES_EXTENDED ) {
                if ( includeDeletions || base != BaseUtils.Base.D.base ) { // note basesAreEqual assigns TRUE to (N,D) as both have simple index -1
                    total += counts[BaseUtils.extendedBaseToBaseIndex(base)];
                }
            }
            depthBySample.put(entry.getKey(),total);
        }
        
        this.updateDepths(depthBySample);
    }

    // returns the index of the last bin that the given depth exceeds the left endpoint of
    private int updateSample(final String sample, final int depth) {
        totalCoverages.put(sample,totalCoverages.get(sample)+depth);

        long[] granularBins = granularHistogramBySample.get(sample);
        for ( int b = 0; b < binLeftEndpoints.length; b ++ ) {
            if ( depth < binLeftEndpoints[b] ) {
                granularBins[b]++;
                return b;
            }
        }

        granularBins[binLeftEndpoints.length]++; // greater than all left-endpoints
        return binLeftEndpoints.length;
    }

    public void merge(DepthOfCoverageStats newStats) {
        this.mergeSamples(newStats);
        if ( this.tabulateLocusCounts && newStats.tabulateLocusCounts ) {
            this.mergeLocusCounts(newStats.getLocusCounts());
        }
        nLoci += newStats.getTotalLoci();
        totalDepthOfCoverage += newStats.getTotalCoverage();
    }

    private void mergeSamples(DepthOfCoverageStats otherStats) {
        Map<String,long[]> otherHistogram = otherStats.getHistograms();
        Map<String,Double> otherMeans = otherStats.getMeans();
        for ( String s : granularHistogramBySample.keySet() ) {
            long[] internalCounts = granularHistogramBySample.get(s);
            long[] externalCounts = otherHistogram.get(s);
            for ( int b = 0; b < internalCounts.length; b++ ) {
                internalCounts[b] += externalCounts[b];
            }

            this.totalCoverages.put(s, this.totalCoverages.get(s) + otherStats.totalCoverages.get(s));
        }
    }

    private void mergeLocusCounts( long[][] otherCounts ) {
        for ( int a = 0; a < locusCoverageCounts.length; a ++ ) {
            for ( int b = 0; b < locusCoverageCounts[0].length; b ++ ) {
                locusCoverageCounts[a][b] += otherCounts[a][b];
            }
        }
    }

    /*
     * Update locus counts -- takes an array in which the number of samples
     * with depth ABOVE [i] is held. So if the bin left endpoints were 2, 5, 10
     * then we'd have an array that represented:
     * [# samples with depth 0 - inf], [# samples with depth 2 - inf],
     * [# samples with depth 5 - inf], [# samples with depth 10-inf];
     *
     * this is
     * @argument cumulativeSamplesByDepthBin - see above
     */
    private void updateLocusCounts(int[] cumulativeSamplesByDepthBin) {
        if ( tabulateLocusCounts ) {
            for ( int bin = 0; bin < cumulativeSamplesByDepthBin.length; bin ++ ) {
                int numSamples = cumulativeSamplesByDepthBin[bin];
                for ( int i = 0; i < numSamples; i ++ ) {
                    locusCoverageCounts[i][bin]++;
                }

                cumulativeSamplesByDepthBin[bin] = 0; // reset counts in advance of next update()
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // ACCESSOR METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    // Return the per-sample histograms of depths
    public Map<String,long[]> getHistograms() {
        return granularHistogramBySample;
    }

    // Return array with the counts of number of bases with >=X samples at >=Y coverage
    public long[][] getLocusCounts() {
        return locusCoverageCounts;
    }

    // Return the key used to interpret the histogram bin edges
    public int[] getEndpoints() {
        return binLeftEndpoints;
    }

    public Map<String,Double> getMeans() {
        HashMap<String,Double> means = new HashMap<String,Double>();
        for ( String s : granularHistogramBySample.keySet() ) {
            means.put(s,( (double)totalCoverages.get(s))/( (double) nLoci ));
        }

        return means;
    }

    // Returns the total coverage (sum of coverage at each base) per sample
    public Map<String,Long> getTotals() {
        return totalCoverages;
    }

    // Returns the number of loci counted
    public long getTotalLoci() {
        return nLoci;
    }

    /**
     * @return  An unordered set of all the samples covered by this stats object
     */
    public Set<String> getAllSamples() {
        return granularHistogramBySample.keySet();
    }

    // Returns the mean of coverage for all loci seen
    public double getTotalMeanCoverage() {
        return ( (double) totalDepthOfCoverage )/ ( (double) nLoci );
    }

    // Returns the total coverage across all samples (sum of coverage at each base)
    public long getTotalCoverage() {
        return totalDepthOfCoverage;
    }

    /**
     * Return a list of the cumulative depth counts as a perportion of total number of loci
     * @param sample
     * @return
     */
    public double[] getCoverageProportions(String sample) {
        long[] hist = granularHistogramBySample.get(sample);
        double[] distribution = new double[hist.length];
        long count = 0;
        for ( int i = hist.length-1; i >= 0; i -- ) {
            count += hist[i];
            distribution[i] = ( (double) count) / nLoci;
        }

        return distribution;
    }

    public int value2bin(int value) {
        for ( int index = 0; index < binLeftEndpoints.length; index++ ) {
            if ( binLeftEndpoints[index] > value ) {
                return index;
            }
        }

        return binLeftEndpoints.length-1;
    }

    /**
     * Returns a count of the number of loci summarized in this object
     */
    public long getNumberOfLociCovered() {
        return nLoci;
    }

}