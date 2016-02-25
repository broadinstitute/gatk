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

/**
 * Worker class to collect insert size metrics.
 * Currently not using the Apache statistics libraries because it requires array structure, which must be converted/constructed
 *    from the map structure returned by countByValue()
 */
public final class InsertSizeMetricsCollectorSpark implements Serializable{
    private static final long serialVersionUID = 1L;

    private final JavaRDD<Integer> lengthRDD;
    private final InsertSizeMetrics metric = new InsertSizeMetrics();
    final Histogram<Integer> htslibHistogram;

    /**
     * Constructor and perform map & reduce.
     * @param filteredReads          reads that pass filters
     * @param samRgRecords           a list of read groups
     * @param HistogramMADTolerance  MAD tolerance when producing histogram plot
     */
    // TODO: * @param accumLevels            accumulation level {SAMPLE, LIBRARY, READGROUP}
    public InsertSizeMetricsCollectorSpark(final JavaRDD<GATKRead> filteredReads,
                                           final long pairsCount,
                                           //final Set<MetricAccumulationLevel> accumLevels,
                                           final List<SAMReadGroupRecord> samRgRecords,
                                           final double HistogramMADTolerance) {

        lengthRDD = filteredReads.map(read -> Math.abs(read.getFragmentLength()) );

        double lengthArr[] = returnRDDArray(lengthRDD);

        DescriptiveStatistics percent = new DescriptiveStatistics(lengthArr);

        metric.READ_PAIRS          = pairsCount;
        metric.MIN_INSERT_SIZE     = (int) percent.getMin();
        metric.MAX_INSERT_SIZE     = (int) percent.getMax();

        metric.WIDTH_OF_10_PERCENT = (int) Math.round(percent.getPercentile(55.0) - percent.getPercentile(45.0));
        metric.WIDTH_OF_20_PERCENT = (int) Math.round(percent.getPercentile(60.0) - percent.getPercentile(40.0));
        metric.WIDTH_OF_30_PERCENT = (int) Math.round(percent.getPercentile(65.0) - percent.getPercentile(35.0));
        metric.WIDTH_OF_40_PERCENT = (int) Math.round(percent.getPercentile(70.0) - percent.getPercentile(30.0));
        metric.WIDTH_OF_50_PERCENT = (int) Math.round(percent.getPercentile(75.0) - percent.getPercentile(25.0));
        metric.WIDTH_OF_60_PERCENT = (int) Math.round(percent.getPercentile(80.0) - percent.getPercentile(20.0));
        metric.WIDTH_OF_70_PERCENT = (int) Math.round(percent.getPercentile(85.0) - percent.getPercentile(15.0));
        metric.WIDTH_OF_80_PERCENT = (int) Math.round(percent.getPercentile(90.0) - percent.getPercentile(10.0));
        metric.WIDTH_OF_90_PERCENT = (int) Math.round(percent.getPercentile(95.0) - percent.getPercentile(5.0));
        metric.WIDTH_OF_99_PERCENT = (int) Math.round(percent.getPercentile(99.5) - percent.getPercentile(0.5));

        metric.MEDIAN_INSERT_SIZE  = percent.getPercentile(50.0);

        final int median = (int) metric.MEDIAN_INSERT_SIZE;
        final JavaRDD<Integer> absLengthDiffRDD = lengthRDD.map(length -> Math.abs(length - median));
        lengthArr = returnRDDArray(absLengthDiffRDD);      // reuse
        percent = new DescriptiveStatistics(lengthArr);    // reuse
        metric.MEDIAN_ABSOLUTE_DEVIATION = percent.getPercentile(50.0);

        // trim down by tolerance
        // mean and SD are severely impacted by outliers, so they are computed here
        // range of read length that is going to be tolerated when computing histogram, mean and SD, but will not affect other metrics.
        final int toleranceRange[] = {median - (int) (HistogramMADTolerance*metric.MEDIAN_ABSOLUTE_DEVIATION),
                                      median + (int) (HistogramMADTolerance*metric.MEDIAN_ABSOLUTE_DEVIATION)};
        final JavaRDD<Integer> trimedLengthRDD = lengthRDD.filter(length -> length>=toleranceRange[0] && length<=toleranceRange[1]);
        lengthArr = returnRDDArray(trimedLengthRDD);

        percent = new DescriptiveStatistics(lengthArr);

        metric.MEAN_INSERT_SIZE   = percent.getMean();
        metric.STANDARD_DEVIATION = percent.getStandardDeviation();

        final Map<Integer, Long> temp = trimedLengthRDD.countByValue();
        final TreeMap<Integer, Long> hist = new TreeMap<>(temp);
        // TODO: fix second argument used in ctor, and potentially change how Histogram is constructed from TreeMap
        htslibHistogram = new Histogram<>("insert_size", "All_Reads.fr_count");
        for(final Integer sz: hist.keySet()){
            htslibHistogram.prefillBins(sz);
            htslibHistogram.increment(sz, hist.get(sz));
        }

        // TODO: change these, unchecked yet
        metric.PAIR_ORIENTATION = SamPairUtil.PairOrientation.FR;
        metric.SAMPLE = null;
        metric.LIBRARY = null;
        metric.READ_GROUP = null;
    }

    // simply return an doubel array from an RDD of integers for using Apache DescriptiveStatistics
    private final double[] returnRDDArray(final JavaRDD<Integer> rdd){
        final List<Integer> tempList = rdd.collect();
        final double lengthArr[] = new double[(int) tempList.size()];
        for(int i=0; i<lengthArr.length; ++i){
            lengthArr[i] = tempList.get(i);
        }
        return lengthArr;
    }

    public final InsertSizeMetrics getMetrics(){
        return metric;
    }

    @VisibleForTesting
    void produceMetricsFile(final MetricsFile<InsertSizeMetrics, Integer> metricsFile) {
        metricsFile.addHistogram(htslibHistogram);
        metricsFile.addMetric(metric);
    }

    /**
     * The reason why this function is named "WithLeftBias" is explained below.
     * Consider this case (simplified): median insert size is 100, and there are exactly 10% valid reads with this
     *   insert size. There are two equally closely valued "bins"--insert size 99 and insert size 101, with 10% and 20%
     *   of valid reads having that particular insert size respectively. Then one is faced with a dilemma about which
     *   one to choose: to use the 99 bin or the 100 bin, and the choice has consequences.
     * Now the name "WithLeftBias" clearly suggests that whenever there is a tie like this, following implementation
     *   always chooses the left bin.
     */
    /*private void collectPercentagesWithLeftBias(InsertSizeMetrics metric) {

        double bins[] = new double[10];
        double percentages[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99};

        Integer keyArray[] = new Integer[hist.keySet().size()];
        int ptr = 0;
        if(hist.get(metric.MEAN_INSERT_SIZE) != null){
            keyArray[ptr] = (int) metric.MEDIAN_INSERT_SIZE;
        }else{

            final Integer left = hist.lowerKey((int)metric.MEDIAN_INSERT_SIZE);
            final Integer right = hist.higherKey((int)metric.MEDIAN_INSERT_SIZE);
            keyArray[ptr] = left;
            keyArray[ptr+1] = right;
        }

        metric.WIDTH_OF_10_PERCENT = (int) bins[0];
        metric.WIDTH_OF_20_PERCENT = (int) bins[1];
        metric.WIDTH_OF_30_PERCENT = (int) bins[2];
        metric.WIDTH_OF_40_PERCENT = (int) bins[3];
        metric.WIDTH_OF_50_PERCENT = (int) bins[4];
        metric.WIDTH_OF_60_PERCENT = (int) bins[5];
        metric.WIDTH_OF_70_PERCENT = (int) bins[6];
        metric.WIDTH_OF_80_PERCENT = (int) bins[7];
        metric.WIDTH_OF_90_PERCENT = (int) bins[8];
        metric.WIDTH_OF_99_PERCENT = (int) bins[9];
    }*/
}