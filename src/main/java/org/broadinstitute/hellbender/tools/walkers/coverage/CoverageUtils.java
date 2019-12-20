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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Utility classes to be used for DepthOfCoverage tool
 */
public class CoverageUtils {

    public enum CountPileupType {
        /**
         * Count all reads independently (even if from the same fragment).
         */
        COUNT_READS,
        /**
         * Count all fragments (even if the reads that compose the fragment are not consistent at that base).
         */
        COUNT_FRAGMENTS,
        /**
         * Count all fragments (but only if the reads that compose the fragment are consistent at that base).
         */
        COUNT_FRAGMENTS_REQUIRE_SAME_BASE
    }

    public static String getTypeID(SAMReadGroupRecord rg, DoCOutputType.Partition type ) {
        switch (type) {
            case sample:
                return rg.getSample();
            case readgroup:
                return rg.getSample()+"_rg_"+rg.getReadGroupId();
            case library:
                return rg.getLibrary();
            case center:
                return rg.getSequencingCenter();
            case platform:
                return rg.getPlatform();
            case sample_by_center:
                return rg.getSample()+"_cn_"+rg.getSequencingCenter();
            case sample_by_platform:
                return rg.getSample()+"_pl_"+rg.getPlatform();
            case sample_by_platform_by_center:
                return rg.getSample()+"_pl_"+rg.getPlatform()+"_cn_"+rg.getSequencingCenter();
            default:
                throw new GATKException(String.format("Invalid aggrigation type %s", type));
        }
    }

    /**
     * Returns the counts of bases from reads with maxMapQ > MAPQ > minMapQ and maxBaseQ > base quality > minBaseQ in the context
     * as an array of ints, indexed by the index fields of BaseUtils. These counts are computed seperately for each combination of
     * {@link DoCOutputType.Partition} and then by sampleID.
     *
     * @param context Alignment context to compute base counts over
     * @param minMapQ nimmun mapping qualiyt
     * @param maxMapQ
     * @param minBaseQ
     * @param maxBaseQ
     * @param countType
     * @param types
     * @param header
     * @return
     */
    public static Map<DoCOutputType.Partition,Map<String,int[]>>
                    getBaseCountsByPartition(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ, CountPileupType countType, Collection<DoCOutputType.Partition> types, SAMFileHeader header) {

        Map<DoCOutputType.Partition,Map<String,int[]>> countsByIDByType = new HashMap<>();
        Map<SAMReadGroupRecord,int[]> countsByRG = getBaseCountsByReadGroup(context, minMapQ, maxMapQ, minBaseQ, maxBaseQ, countType, header);
        for (DoCOutputType.Partition t : types ) {
            // iterate through the read group counts and build the type associations
            for ( Map.Entry<SAMReadGroupRecord,int[]> readGroupCountEntry : countsByRG.entrySet() ) {
                String typeID = getTypeID(readGroupCountEntry.getKey(),t);

                if ( ! countsByIDByType.keySet().contains(t) ) {
                    countsByIDByType.put(t,new HashMap<>());
                }

                if ( ! countsByIDByType.get(t).keySet().contains(typeID) ) {
                    countsByIDByType.get(t).put(typeID,readGroupCountEntry.getValue().clone());
                } else {
                    addCounts(countsByIDByType.get(t).get(typeID),readGroupCountEntry.getValue());
                }
            }
        }


        return countsByIDByType;
    }

    private static void addCounts(int[] updateMe, int[] leaveMeAlone ) {
        for ( int index = 0; index < leaveMeAlone.length; index++ ) {
            updateMe[index] += leaveMeAlone[index];
        }
    }

    /**
     * Takes an AlignmentContext object and extracts all the reads that pass the provided filters, and then returns a
     * count breakdown for each base (plus D and N bases) present at the site.
     *
     * NOTE: this currently doesn't support counts by fragments as was the case in gatk3
     */
    private static Map<SAMReadGroupRecord,int[]> getBaseCountsByReadGroup(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ, CountPileupType countType, SAMFileHeader header) {
        Map<SAMReadGroupRecord, int[]> countsByRG = new HashMap<>();

        Map<String, int[]> countsByRGName = new HashMap<>();
        Map<String, SAMReadGroupRecord> RGByName = new HashMap<String, SAMReadGroupRecord>();

        List<PileupElement> countPileup = new LinkedList<>();

        switch (countType) {

            case COUNT_READS:
                for (PileupElement e : context.getBasePileup()) {
                    if (filterElementFromCount(e, minMapQ, maxMapQ, minBaseQ, maxBaseQ)) {
                        countPileup.add(e);
                    }
                }
                break;

            // TODO see reconcile FragmentUtils.create() and its various idiosyncrasies to reenable this feature
            case COUNT_FRAGMENTS: // ignore base identities and put in FIRST base that passes filters:
                throw new GATKException("Fragment based counting is currently unsupported");

            case COUNT_FRAGMENTS_REQUIRE_SAME_BASE:
                throw new GATKException("Fragment based counting is currently unsupported");

            default:
                throw new UserException("Must use valid CountPileupType");
        }

        for (PileupElement e : countPileup) {
            SAMReadGroupRecord readGroup = ReadUtils.getSAMReadGroupRecord(e.getRead(), header);
            Utils.nonNull(readGroup, String.format("Read %s was missing read group information", e.getRead()));

            // uniqueReadGroupID is unique across the library, read group ID, and the sample
            String uniqueReadGroupId = readGroup.getSample() + "_" + readGroup.getReadGroupId() + "_" + readGroup.getLibrary() + "_" + readGroup.getPlatformUnit();
            int[] counts = countsByRGName.get(uniqueReadGroupId);
            if (counts == null) {
                counts = new int[6];
                countsByRGName.put(uniqueReadGroupId, counts);
                RGByName.put(uniqueReadGroupId, readGroup);
            }

            updateCounts(counts, e);
        }

        for (String readGroupId : RGByName.keySet()) {
            countsByRG.put(RGByName.get(readGroupId), countsByRGName.get(readGroupId));
        }

        return countsByRG;
    }

    // Applies the provided mapping and base quality filters to the provided read
    private static boolean filterElementFromCount(PileupElement e, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ) {
        return (e.getMappingQual() >= minMapQ && e.getMappingQual() <= maxMapQ && ( e.getQual() >= minBaseQ && e.getQual() <= maxBaseQ || e.isDeletion() ));
    }

    private static void updateCounts(int[] counts, PileupElement e) {
        if ( e.isDeletion() ) {
            counts[BaseUtils.Base.D.ordinal()]++;
        } else if ( BaseUtils.basesAreEqual(BaseUtils.Base.N.base, e.getBase()) ) {
            counts[BaseUtils.Base.N.ordinal()]++;
        } else {
            try {
                counts[BaseUtils.simpleBaseToBaseIndex(e.getBase())]++;
            } catch (ArrayIndexOutOfBoundsException exc) {
                throw new GATKException("Expected a simple base, but actually received"+(char)e.getBase());
            }
        }
    }

    /**
     * This method is used to construct the bins we use to compute coverage histograms.
     *
     * Returns the coverage "leftEndpoints" which correspond to the lowest value inclusive to be included in the given
     * histogram bins. This method will attempt to distribute the bins between the specified values based a logarithmic
     * distribution?
     *
     * @param lower lower bound of the second bin in the histogram (the first bin includes everything below this value)
     * @param upper lower bound of the last bin in the histogram (the last bin will correspond to everything over this value)
     * @param nBins the number of bins to try and fit between lower and upper (must be < upper - lower)
     * @return an array corresponding to all of the inclusive left edges for a histogram with the specified parameters.
     */
    public static int[] calculateCoverageHistogramBinEndpoints(int lower, int upper, int nBins) {
        if ( nBins > upper - lower || lower < 1 ) {
            throw new UserException.BadInput("the start must be at least 1 and the number of bins may not exceed stop - start");
        }

        int[] binLeftEndpoints = new int[nBins+1];
        binLeftEndpoints[0] = lower;

        int length = upper - lower;
        double scale = Math.log10((double) length)/nBins;

        for ( int b = 1; b < nBins ; b++ ) {
            int leftEnd = lower + (int) Math.floor(Math.pow(10.0,(b-1.0)*scale));
            // todo -- simplify to length^(scale/bins); make non-constant to put bin ends in more "useful"
            // todo -- positions on the number line
            while ( leftEnd <= binLeftEndpoints[b-1] ) {
                leftEnd++;
            }

            binLeftEndpoints[b] = leftEnd;
        }

        binLeftEndpoints[binLeftEndpoints.length-1] = upper;

        return binLeftEndpoints;
    }

    /*
     * @updateTargetTable
     * The idea is to have counts for how many *targets* have at least K samples with
     * median coverage of at least X.
     * To that end:
     * Iterate over the samples the DOCS object, determine how many there are with
     * median coverage > leftEnds[0]; how many with median coverage > leftEnds[1]
     * and so on. Then this target has at least N, N-1, N-2, ... 1, 0 samples covered
     * to leftEnds[0] and at least M,M-1,M-2,...1,0 samples covered to leftEnds[1]
     * and so on.
     */
    public static void updateTargetTable(int[][] table, DepthOfCoverageStats stats) {
        int[] cutoffs = stats.getEndpoints();
        int[] countsOfMediansAboveCutoffs = new int[cutoffs.length+1]; // 0 bin to catch everything

        for ( String s : stats.getAllSamples() ) {
            int medianBin = getQuantile(stats.getHistograms().get(s),0.5);
            for ( int i = 0; i <= medianBin; i ++) {
                countsOfMediansAboveCutoffs[i]++;
            }
        }

        for ( int medianBin = 0; medianBin < countsOfMediansAboveCutoffs.length; medianBin++) {
            for ( ; countsOfMediansAboveCutoffs[medianBin] > 0; countsOfMediansAboveCutoffs[medianBin]-- ) {
                table[countsOfMediansAboveCutoffs[medianBin]-1][medianBin]++;
                // the -1 is due to counts being 1-based and offsets being 0-based
            }
        }
    }

    /**
     * Returns the percentage of the overall histogram dataset that is beyond a certain bin.
     *
     * @param histogram histogram to test
     * @param bin bin to evaluate
     * @return value from 0-100 corresponding to the portion of the total count over the given bin
     */
    public static double getPctBasesAbove(long[] histogram, int bin) {
        long below = 0;
        long above = 0;
        for ( int index = 0; index < histogram.length; index++) {
            if ( index < bin ) {
                below+=histogram[index];
            } else {
                above+=histogram[index];
            }
        }

        return 100*( (double) above )/( above + below );
    }

    /**
     * Returns the bin that corresponds to provided quantile value.
     *
     * @param histogram histogram to compute quantile for
     * @param prop Quantile proportion
     * @return index of histogram bin containing at lest total*prop values in lower indexed bins
     */
    public static int getQuantile(long[] histogram, double prop) {
        Utils.validate(prop >= 0 && prop <= 1, "Quantile proportion must fall between 0.0 and 1.0 inclusive");
        long total = MathUtils.sum(histogram);

        long counts = 0;
        int bin = -1;
        while ( counts < prop*total ) {
            counts += histogram[bin+1];
            bin++;
        }

        return bin == -1 ? 0 : bin;
    }

}
