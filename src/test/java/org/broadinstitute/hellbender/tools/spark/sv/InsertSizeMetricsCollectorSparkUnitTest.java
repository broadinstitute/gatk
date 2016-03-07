package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.testng.Assert;
import com.google.common.collect.Sets;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;

/**
 * Testing several of the utility functions defined and used in InsertSizeMetricsCollectorSpark.
 */
public final class InsertSizeMetricsCollectorSparkUnitTest {

    /**
     * Testing "long[] InsertSizeMetricsCollectorSpark::computeRanges(final Map<Integer, Long>, final int, final double)",
     *     using actual valid reads (13 of them) information from "insert_size_metrics_test.bam", where inferred length are
     *     {36 36 36 38 38 40 41 41 41 41 44 44 45} with median == 40.
     */
    @Test(groups = "sv")
    public void testRanges() throws Exception {

        final Map<Integer, Long> map = new TreeMap<>();
        map.put(36, 3L);
        map.put(38, 2L);
        map.put(40, 1L);
        map.put(41, 4L);
        map.put(44, 2L);
        map.put(45, 1L);

        Histogram<Integer> hist = new Histogram<>("dummy", "test");
        for(final int size: map.keySet()){
            hist.prefillBins(size);
            hist.increment(size, map.get(size));
        }

        long[] calculatedBinWidths = InsertSizeMetricsCollectorSpark.computeRanges(hist, 41, 13);
        long[] expectedBinWidths = {1L, 1L, 1L, 7L, 7L, 7L, 9L, 11L, 11L, 11L};
        Assert.assertEquals(calculatedBinWidths, expectedBinWidths);
    }

    /**
     * Tests "ArrayList<Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> aggregateHistograms(final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>,
                                                                                                                            final Set<MetricAccumulationLevel>))
     * Test case: 1 sample, 2 libraries
     *            3 read groups, 2 belonging to 1 library, 1 belonging to the other
     * @throws Exception
     */
    @Test(groups = "sv")
    public void testAggregator() throws Exception{


        final InsertSizeMetricsCollectorSpark.GroupMetaInfo readGroup1A = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("testSample",
                                                                                                                            "testLibrary1",
                                                                                                                            "testReadGroup1A",
                                                                                                                            MetricAccumulationLevel.READ_GROUP);
        final InsertSizeMetricsCollectorSpark.GroupMetaInfo readGroup1B = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("testSample",
                                                                                                                            "testLibrary1",
                                                                                                                            "testReadGroup1B",
                                                                                                                            MetricAccumulationLevel.READ_GROUP);
        final InsertSizeMetricsCollectorSpark.GroupMetaInfo readGroup2 = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("testSample",
                                                                                                                            "testLibrary2",
                                                                                                                            "testReadGroup2",
                                                                                                                            MetricAccumulationLevel.READ_GROUP);

        final InsertSizeMetricsCollectorSpark.GroupMetaInfo library1 = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("testSample",
                                                                                                                         "testLibrary1",
                                                                                                                         null,
                                                                                                                         MetricAccumulationLevel.LIBRARY);
        final InsertSizeMetricsCollectorSpark.GroupMetaInfo library2 = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("testSample",
                                                                                                                        "testLibrary2",
                                                                                                                        null,
                                                                                                                        MetricAccumulationLevel.LIBRARY);
        final InsertSizeMetricsCollectorSpark.GroupMetaInfo sample1 = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("testSample",
                                                                                                                        null,
                                                                                                                        null,
                                                                                                                        MetricAccumulationLevel.SAMPLE);
        final InsertSizeMetricsCollectorSpark.GroupMetaInfo allReads = new InsertSizeMetricsCollectorSpark.GroupMetaInfo(null,
                                                                                                                         null,
                                                                                                                         null,
                                                                                                                         MetricAccumulationLevel.ALL_READS);

        final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> histogramsOfReadGroup1A = new HashMap<>();
        final SortedMap<Integer, Long> testHist1A = new TreeMap<>();
        testHist1A.put(1,2L);
        testHist1A.put(3,4L);
        histogramsOfReadGroup1A.put(SamPairUtil.PairOrientation.FR, testHist1A);

        final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> histogramsOfReadGroup1B = new HashMap<>();
        final SortedMap<Integer, Long> testHist1B = new TreeMap<>();
        testHist1B.put(50,60L);
        histogramsOfReadGroup1B.put(SamPairUtil.PairOrientation.FR, testHist1B);

        final Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> histogramsOfReadGroup2 = new HashMap<>();
        final SortedMap<Integer, Long> testHist2 = new TreeMap<>();
        testHist2.put(100,200L);
        testHist2.put(1,  500L);
        histogramsOfReadGroup2.put(SamPairUtil.PairOrientation.FR, testHist2);

        Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> histogramsOfAllTestReadGroups = new HashMap<>();
        histogramsOfAllTestReadGroups.put(readGroup1A, histogramsOfReadGroup1A);
        histogramsOfAllTestReadGroups.put(readGroup1B, histogramsOfReadGroup1B);
        histogramsOfAllTestReadGroups.put(readGroup2 , histogramsOfReadGroup2 );

        final Set<MetricAccumulationLevel> accumLevels = Sets.newHashSet(MetricAccumulationLevel.READ_GROUP, MetricAccumulationLevel.LIBRARY, MetricAccumulationLevel.SAMPLE, MetricAccumulationLevel.ALL_READS);
        final ArrayList<Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>> listOfResults = InsertSizeMetricsCollectorSpark.aggregateHistograms(histogramsOfAllTestReadGroups, accumLevels);

        int sz = 0;
        for(final Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> entry : listOfResults) { sz += entry.size(); }
        Assert.assertEquals(sz, 7);

        Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> histogramsOfLibraries = listOfResults.get(1);
        Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> histogramsOfSamples   = listOfResults.get(2);
        Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> histogramsOfAllReads  = listOfResults.get(3);

        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(1), Long.valueOf(502L));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(3), Long.valueOf(4L));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(50), Long.valueOf(60L));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(100), Long.valueOf(200L));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(123), null);

        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(1),
                            histogramsOfSamples.get(sample1).get(SamPairUtil.PairOrientation.FR).get(1));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(3),
                            histogramsOfSamples.get(sample1).get(SamPairUtil.PairOrientation.FR).get(3));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(50),
                            histogramsOfSamples.get(sample1).get(SamPairUtil.PairOrientation.FR).get(50));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(100),
                            histogramsOfSamples.get(sample1).get(SamPairUtil.PairOrientation.FR).get(100));
        Assert.assertEquals(histogramsOfAllReads.get(allReads).get(SamPairUtil.PairOrientation.FR).get(123),
                            histogramsOfSamples.get(sample1).get(SamPairUtil.PairOrientation.FR).get(123));

        Assert.assertEquals(histogramsOfLibraries.get(library1).get(SamPairUtil.PairOrientation.FR).get(1), Long.valueOf(2L));
        Assert.assertEquals(histogramsOfLibraries.get(library1).get(SamPairUtil.PairOrientation.FR).get(3), Long.valueOf(4L));
        Assert.assertEquals(histogramsOfLibraries.get(library1).get(SamPairUtil.PairOrientation.FR).get(50), Long.valueOf(60L));
        Assert.assertEquals(histogramsOfLibraries.get(library1).get(SamPairUtil.PairOrientation.FR).get(100), null);

        Assert.assertEquals(histogramsOfLibraries.get(library2).get(SamPairUtil.PairOrientation.FR).get(1), Long.valueOf(500L));
        Assert.assertEquals(histogramsOfLibraries.get(library2).get(SamPairUtil.PairOrientation.FR).get(3), null);
        Assert.assertEquals(histogramsOfLibraries.get(library2).get(SamPairUtil.PairOrientation.FR).get(50), null);
        Assert.assertEquals(histogramsOfLibraries.get(library2).get(SamPairUtil.PairOrientation.FR).get(100), Long.valueOf(200L));
    }

    /**
     * Tests void convertSortedMapToHTSHistogram(Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>>>,
                                                 final Map<GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>>,
                                                 final double)
     * @throws Exception
     */
    @Test(groups = "sv")
    public void testSortedMapToHTSJDKHistogramConverter() throws Exception{

        final InsertSizeMetricsCollectorSpark.GroupMetaInfo testGroup = new InsertSizeMetricsCollectorSpark.GroupMetaInfo("sample", "library", "readGroup", MetricAccumulationLevel.READ_GROUP);
        final SortedMap<Integer, Long> hist = new TreeMap<>();
        hist.put(36, 3L);
        hist.put(38, 2L);
        hist.put(40, 1L);
        hist.put(41, 4L);
        hist.put(44, 2L);
        hist.put(45, 1L);
        Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>>> testHist = new HashMap<>();
        Map<SamPairUtil.PairOrientation, SortedMap<Integer, Long>> FRMap = new HashMap<>();
        FRMap.put(SamPairUtil.PairOrientation.FR, hist);
        testHist.put(testGroup, FRMap);

        Map<InsertSizeMetricsCollectorSpark.GroupMetaInfo, Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>>> testHTSJDKHistogram = new HashMap<>();
        Histogram<Integer> htsjdkHistogram = new Histogram<>("dummy", "test");
        InsertSizeMetrics metrics = new InsertSizeMetrics();
        Map<SamPairUtil.PairOrientation, Tuple2<Histogram<Integer>, InsertSizeMetrics>> FRHist = new HashMap<>();
        FRHist.put(SamPairUtil.PairOrientation.FR, new Tuple2<>(htsjdkHistogram, metrics));

        testHTSJDKHistogram.put(testGroup, FRHist);

        InsertSizeMetricsCollectorSpark.convertSortedMapToHTSHistogram(testHTSJDKHistogram, testHist, 10.0);

        final Histogram<Integer> filledHistogram = testHTSJDKHistogram.get(testGroup).get(SamPairUtil.PairOrientation.FR)._1();
        final InsertSizeMetrics filledMetrics = testHTSJDKHistogram.get(testGroup).get(SamPairUtil.PairOrientation.FR)._2();

        Assert.assertEquals(filledHistogram.getMax(), 45.0);
        Assert.assertEquals(filledHistogram.getMin(), 36.0);
        Assert.assertEquals(filledHistogram.getCount(), 13.0);
        Assert.assertEquals(filledHistogram.getMean(), 40.1, 0.05);
        Assert.assertEquals(filledHistogram.getMeanBinSize(), 2.17, 0.05);
        Assert.assertEquals(filledHistogram.getStandardDeviation(), 3.1, 0.05);
        Assert.assertEquals(filledHistogram.getMedian(), 41.0);
        Assert.assertEquals(filledHistogram.getMedianBinSize(), 2.0);
        Assert.assertEquals(filledHistogram.getMedianAbsoluteDeviation(), 3.0);
        Assert.assertEquals(filledHistogram.getSum(), 521.0);
        Assert.assertEquals(filledHistogram.getMode(), 41.0);
        Assert.assertEquals(filledHistogram.getPercentile(0.10), 36.0);
        Assert.assertEquals(filledHistogram.getPercentile(0.50), 41.0);
        Assert.assertEquals(filledHistogram.getPercentile(0.90), 44.0);

        Assert.assertEquals(filledMetrics.READ_GROUP, "readGroup");
        Assert.assertEquals(filledMetrics.LIBRARY, "library");
        Assert.assertEquals(filledMetrics.SAMPLE, "sample");

        Assert.assertEquals(filledMetrics.READ_PAIRS, 13);
        Assert.assertEquals(filledMetrics.MIN_INSERT_SIZE, 36);
        Assert.assertEquals(filledMetrics.MAX_INSERT_SIZE, 45);
        Assert.assertEquals(filledMetrics.MEAN_INSERT_SIZE, 40.1, 0.05);
        Assert.assertEquals(filledMetrics.STANDARD_DEVIATION, 3.1, 0.05);
        Assert.assertEquals(filledMetrics.MEDIAN_INSERT_SIZE, 41.0);
        Assert.assertEquals(filledMetrics.MEDIAN_ABSOLUTE_DEVIATION, 3.0);

        Assert.assertEquals(filledMetrics.WIDTH_OF_10_PERCENT, 1);
        Assert.assertEquals(filledMetrics.WIDTH_OF_20_PERCENT, 1);
        Assert.assertEquals(filledMetrics.WIDTH_OF_30_PERCENT, 1);
        Assert.assertEquals(filledMetrics.WIDTH_OF_40_PERCENT, 7);
        Assert.assertEquals(filledMetrics.WIDTH_OF_50_PERCENT, 7);
        Assert.assertEquals(filledMetrics.WIDTH_OF_60_PERCENT, 7);
        Assert.assertEquals(filledMetrics.WIDTH_OF_70_PERCENT, 9);
        Assert.assertEquals(filledMetrics.WIDTH_OF_80_PERCENT, 11);
        Assert.assertEquals(filledMetrics.WIDTH_OF_90_PERCENT, 11);
        Assert.assertEquals(filledMetrics.WIDTH_OF_99_PERCENT, 11);
    }
}
