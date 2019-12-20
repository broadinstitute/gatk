/*
* Copyright 2012-2016 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.utils.BaseUtils;

import java.util.*;

/**
 * A class for storing summarized
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
     * @param includeDeletions
     * @param dontComputeLocusTable
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

    // Adds to granularHistogramBySample a per-locus histogram for
    void initializeSample(String sample) {
        if ( granularHistogramBySample.containsKey(sample) ) {
            return;
        }
        granularHistogramBySample.put(sample, new long[this.binLeftEndpoints.length+1]);
        totalCoverages.put(sample,0L);
    }

    // Create the per-locus
    void initializeLocusCounts() {
        locusCoverageCounts = new long[granularHistogramBySample.size()][binLeftEndpoints.length+1];
        locusHistogram = new int[binLeftEndpoints.length+1];
        tabulateLocusCounts = true;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // UPDATE METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    private void updateDepths(Map<String, Integer> depthBySample) {
        int b;
        for ( String sample : granularHistogramBySample.keySet() ) {
            if ( depthBySample.containsKey(sample) ) {
                b = updateSample(sample,depthBySample.get(sample));
                totalLocusDepth += depthBySample.get(sample);
            } else {
                b = updateSample(sample,0);
            }

            if ( tabulateLocusCounts ) {
                for ( int i = 0; i <= b; i ++ ) {
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
        // todo -- do we want to do anything special regarding base count or deletion statistics?
        HashMap<String,Integer> depthBySample = new HashMap<String,Integer>();
        // todo -- needs fixing with advent of new baseutils functionality using ENUMS and handling N,D
        for ( String s : countsBySample.keySet() ) {
            int total = 0;
            int[] counts = countsBySample.get(s);
            for ( byte base : BaseUtils.BASES_EXTENDED ) {
                if ( includeDeletions || ! ( base == BaseUtils.Base.D.base) ) { // note basesAreEqual assigns TRUE to (N,D) as both have simple index -1
                    total += counts[BaseUtils.extendedBaseToBaseIndex(base)];
                }
            }
            depthBySample.put(s,total);
        }
        
        this.updateDepths(depthBySample);
    }

    private int updateSample(String sample, int depth) {
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

    public Map<String,long[]> getHistograms() {
        return granularHistogramBySample;
    }

    public long[][] getLocusCounts() {
        return locusCoverageCounts;
    }

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

    public Map<String,Long> getTotals() {
        return totalCoverages;
    }

    public long getTotalLoci() {
        return nLoci;
    }

    /**
     * @return  An unordered set of all the samples covered by this stats object
     */
    public Set<String> getAllSamples() {
        return granularHistogramBySample.keySet();
    }

    public double getTotalMeanCoverage() {
        return ( (double) totalDepthOfCoverage )/ ( (double) nLoci );
    }

    public long getTotalCoverage() {
        return totalDepthOfCoverage;
    }

    /**
     * Return a list of the counts
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