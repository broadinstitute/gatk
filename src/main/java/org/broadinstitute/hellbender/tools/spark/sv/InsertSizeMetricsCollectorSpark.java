package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.SortedMap;

/**
 * Worker class to collect insert size metrics.
 */
public final class InsertSizeMetricsCollectorSpark implements Serializable{
    private static final long serialVersionUID = 1L;

    private final InsertSizeMetrics metrics = new InsertSizeMetrics();
    private Histogram<Integer> htslibHistogram;

    /**
     * Constructor and do collection work.
     * @param filteredReads          reads that pass filters
     * @param samRgRecords           a list of read groups
     * @param HistogramMADTolerance  MAD tolerance when producing histogram plot
     */
    // TODO: * @param accumLevels            accumulation level {SAMPLE, LIBRARY, READGROUP}
    public InsertSizeMetricsCollectorSpark(final JavaRDD<GATKRead> filteredReads,
                                           //final Set<MetricAccumulationLevel> accumLevels,
                                           final List<SAMReadGroupRecord> samRgRecords,
                                           final double HistogramMADTolerance) {

        final JavaRDD<Integer> lengthRDD = filteredReads.map(read -> Math.abs(read.getFragmentLength()) );
        lengthRDD.cache();

        metrics.READ_PAIRS = filteredReads.count();

        collectLevelInfo();

        collectMedianAndMAD(lengthRDD);

        collectSymmetricalBinWidth(lengthRDD);

        trimHistogramAndCollectMeanAndSD(lengthRDD, HistogramMADTolerance);
    }

    public final InsertSizeMetrics getMetrics(){
        return metrics;
    }

    public final Histogram<Integer> getHitogram(){
        return htslibHistogram;
    }

    // simply return an double array from an RDD of integers for using Apache DescriptiveStatistics
    private static final double[] getRDDArray(final JavaRDD<Integer> rdd){
        final List<Integer> tempList = rdd.collect();
        final double lengthArr[] = new double[tempList.size()];
        for(int i=0; i<lengthArr.length; ++i){
            lengthArr[i] = tempList.get(i);
        }
        return lengthArr;
    }

    private void collectMedianAndMAD(final JavaRDD<Integer> lengthRDD){

        double lengthArr[] = getRDDArray(lengthRDD);

        final DescriptiveStatistics percent = new DescriptiveStatistics(lengthArr);

        metrics.MIN_INSERT_SIZE     = (int) percent.getMin();
        metrics.MAX_INSERT_SIZE     = (int) percent.getMax();

        metrics.MEDIAN_INSERT_SIZE  = percent.getPercentile(50.0);

        final int median = (int) metrics.MEDIAN_INSERT_SIZE;
        double lengthDiffArr[] = getRDDArray( lengthRDD.map(length -> Math.abs(length - median)) );
        final DescriptiveStatistics percentOnDiff = new DescriptiveStatistics(lengthDiffArr);
        metrics.MEDIAN_ABSOLUTE_DEVIATION = percentOnDiff.getPercentile(50.0);
    }

    private void collectSymmetricalBinWidth(final JavaRDD<Integer> lengthRDD){

        final long bin_widths[] = computeRanges(lengthRDD.countByValue(),
                                                (int) metrics.MEDIAN_INSERT_SIZE,
                                                (double) lengthRDD.count());

        // extend left and right symmetrically half of these width values from the median,
        //   the indicated percentages of valid reads will have fragment length in that range.
        metrics.WIDTH_OF_10_PERCENT = (int) bin_widths[0];
        metrics.WIDTH_OF_20_PERCENT = (int) bin_widths[1];
        metrics.WIDTH_OF_30_PERCENT = (int) bin_widths[2];
        metrics.WIDTH_OF_40_PERCENT = (int) bin_widths[3];
        metrics.WIDTH_OF_50_PERCENT = (int) bin_widths[4];
        metrics.WIDTH_OF_60_PERCENT = (int) bin_widths[5];
        metrics.WIDTH_OF_70_PERCENT = (int) bin_widths[6];
        metrics.WIDTH_OF_80_PERCENT = (int) bin_widths[7];
        metrics.WIDTH_OF_90_PERCENT = (int) bin_widths[8];
        metrics.WIDTH_OF_99_PERCENT = (int) bin_widths[9];
    }

    /**
     * This is where actual work for collecting the symmetrical bin boundaries around the median is done.
     * @param unsortedHist unsorted histogram containing mapping form insert size to count
     * @param start        starting value for searching
     * @param totalCount   total count of valid fragments
     * @return             returns the bin widths (right_edge - left_edge + 1) collectively as an array
     */
    @VisibleForTesting
    static long[] computeRanges(final Map<Integer, Long> unsortedHist, final int start, final double totalCount){

        final SortedMap<Integer, Long> hist = new TreeMap<>(unsortedHist);

        double sum = 0.0;  // for calculating coverage, stored as sum to avoid frequent casting

        int left  = start; // left and right boundaries of histogram bins
        int right = left;  //      start from median, and gradually open up

        long bin_widths[] = new long[10];   // for storing distance between left and right boundaries of histogram bins
                                            // dimension is 10 because metrics requires 10 histogram bin width values.
        int i = 0;
        int j = 0;                          // represent lowest and highest indices of bin_widths that needs to be updated

        while(i<10){                        // until all width values are computed
            Long sz1 = hist.get(left);
            Long sz2 = (left!=right)? hist.get(right) : 0L;
            if(null!=sz1) { sum += sz1; }   // since left and right are incremented/decremented by 1,
            if(null!=sz2) { sum += sz2; }   //    they may end up not in hist.keySet()

            j = (int) (10.*sum/totalCount); // if coverage increased by enough value, update necessary ones
            for(int k=i; k<j; ++k){
                bin_widths[k] = right - left + 1;
            }
            i = j;                          // and update pointers

            --left;
            ++right;
        }

        return bin_widths;
    }

    /**
     * Trim down width of histogram by tolerance.
     * Mean and SD are severely impacted by outliers, so they are computed here after trimming.
     * @param HistogramMADTolerance how many MADs from the median will be tolerated when down-trimming
     */
    private void trimHistogramAndCollectMeanAndSD(final JavaRDD<Integer> lengthRDD, final double HistogramMADTolerance){

        // range of read length that is going to be tolerated when computing histogram, mean and SD, but will not affect other metrics.
        final int lengthToleranceLowerEnd = (int) metrics.MEDIAN_INSERT_SIZE - (int) (HistogramMADTolerance * metrics.MEDIAN_ABSOLUTE_DEVIATION);
        final int lengthToleranceHigherEnd = (int) metrics.MEDIAN_INSERT_SIZE + (int) (HistogramMADTolerance * metrics.MEDIAN_ABSOLUTE_DEVIATION);
        final JavaRDD<Integer> trimmedLengthRDD = lengthRDD.filter(length -> length >= lengthToleranceLowerEnd)
                                                           .filter(length -> length <= lengthToleranceHigherEnd);
        final double lengthArr[] = getRDDArray(trimmedLengthRDD);

        final DescriptiveStatistics percent = new DescriptiveStatistics(lengthArr);

        metrics.MEAN_INSERT_SIZE   = percent.getMean();
        metrics.STANDARD_DEVIATION = percent.getStandardDeviation();

        final SortedMap<Integer, Long> hist = new TreeMap<>(lengthRDD.countByValue());
        // TODO: fix second argument used in ctor, and potentially change how Histogram is constructed from TreeMap
        htslibHistogram = new Histogram<>("insert_size", "All_Reads.fr_count");
        for(final int size: hist.keySet()){
            htslibHistogram.prefillBins(size);
            htslibHistogram.increment(size, hist.get(size));
        }
    }

    /**
     * Collect bulk information: pair orientation, sample name, library name and read groups.
     */
    private void collectLevelInfo(){
        // TODO: change these, unchecked yet
        metrics.PAIR_ORIENTATION = SamPairUtil.PairOrientation.FR;
        metrics.SAMPLE = null;
        metrics.LIBRARY = null;
        metrics.READ_GROUP = null;
    }

    @VisibleForTesting
    void produceMetricsFile(final MetricsFile<InsertSizeMetrics, Integer> metricsFile) {
        metricsFile.addHistogram(htslibHistogram);
        metricsFile.addMetric(metrics);
    }
}