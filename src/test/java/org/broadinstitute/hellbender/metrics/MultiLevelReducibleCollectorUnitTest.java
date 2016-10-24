package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.Serializable;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeSet;

/**
 * Test all of the functionality of MultiLevelReducibleCollector. This includes duplicates of
 * all of the basic functionality of multi-level distribution based on the MultiLevelCollector
 * tests, with additional test code to verify the combine functionality used for combing/reducing
 * multiple collectors into a single collector.
 */
public final class MultiLevelReducibleCollectorUnitTest {

    public static File TESTFILE = new File("src/test/resources/org/broadinstitute/hellbender/metrics/test.sam");

    public String noneOrStr(final String str) {
        return str == null ? "" : str;
    }

    static final class TestArg {
        public final SAMRecord samRecord;
        public final ReferenceSequence refSeq;

        public TestArg(final SAMRecord samRecord, final ReferenceSequence refSeq) {
            this.samRecord = samRecord;
            this.refSeq    = refSeq;
        }
    }

    /** We will just Tally up the number of times records were added to this metric and change FINISHED
     * to true when FINISHED is called
     */
    static final class TotalNumberMetric extends MultiLevelMetrics implements Serializable {
        private static final long serialVersionUID = 1L;

        /** The number of these encountered **/
        public Integer TALLY = 0;
        public boolean FINISHED = false;
    }

    private class RecordCountPerUnitCollector implements PerUnitMetricCollector<TotalNumberMetric, Integer, TestArg>{
        private static final long serialVersionUID = 1L;
        private final TotalNumberMetric metric;

        RecordCountMultiLevelCollector outerCollector;

        public RecordCountPerUnitCollector(
                final RecordCountMultiLevelCollector outerCollector, final String sample, final String library, final String readGroup) {
            this.outerCollector  = outerCollector;
            metric = new TotalNumberMetric();
            metric.SAMPLE     = sample;
            metric.LIBRARY    = library;
            metric.READ_GROUP = readGroup;
            outerCollector.unitsToMetrics.put(noneOrStr(sample) + "_" + noneOrStr(library) + "_" + noneOrStr(readGroup), metric);
        }

        @Override
        public void acceptRecord(final TestArg args) {
            outerCollector.numProcessed += 1;
            metric.TALLY += 1;
            if(metric.SAMPLE != null) {
                Assert.assertEquals(metric.SAMPLE, args.samRecord.getReadGroup().getSample());
            }
            if(metric.LIBRARY != null) {
                Assert.assertEquals(metric.LIBRARY, args.samRecord.getReadGroup().getLibrary());
            }

            if(metric.READ_GROUP != null) {
                Assert.assertEquals(metric.READ_GROUP, args.samRecord.getReadGroup().getPlatformUnit());
            }
        }

        @Override
        public void finish() {
            metric.FINISHED = true;
        }

        @Override
        public void addMetricsToFile(final MetricsFile<TotalNumberMetric, Integer> totalNumberMetricIntegerMetricsFile) {
            totalNumberMetricIntegerMetricsFile.addMetric(metric);
        }

        public RecordCountPerUnitCollector combine(RecordCountPerUnitCollector source) {
            Assert.assertEquals(this.metric.FINISHED, true);
            Assert.assertEquals(source.metric.FINISHED, true);
            Assert.assertEquals(this.metric.SAMPLE, source.metric.SAMPLE);
            Assert.assertEquals(this.metric.LIBRARY, source.metric.LIBRARY);
            Assert.assertEquals(this.metric.READ_GROUP, source.metric.READ_GROUP);

            metric.TALLY += source.metric.TALLY;
            return this;
        }
    }

    class RecordCountMultiLevelCollector extends MultiLevelReducibleCollector<TotalNumberMetric, Integer, TestArg, RecordCountPerUnitCollector> {
        private static final long serialVersionUID = 1L;

        public RecordCountMultiLevelCollector(
                final Set<MetricAccumulationLevel> accumulationLevels,
                final List<SAMReadGroupRecord> samRgRecords) {
            setup(accumulationLevels, samRgRecords);
        }

        //The number of times records were accepted by a RecordCountPerUnitCollectors (note since the same
        //samRecord might be aggregated by multiple PerUnit collectors, this may be greater than the number of
        //records in the file
        private int numProcessed = 0;

        public int getNumProcessed() {
            return numProcessed;
        }

        private Map<String, TotalNumberMetric> unitsToMetrics = new LinkedHashMap<>();

        public Map<String, TotalNumberMetric> getUnitsToMetrics() {
            return unitsToMetrics;
        }

        public void setUnitsToMetrics(Map<String, TotalNumberMetric> inMap) {
            unitsToMetrics = inMap;
        }

        @Override
        protected TestArg makeArg(final SAMRecord samRec, final ReferenceSequence refSeq) {
            return new TestArg(samRec, refSeq);
        }

        @Override
        protected RecordCountPerUnitCollector makeChildCollector(
                final String sample, final String library, final String readGroup) {
            return new RecordCountPerUnitCollector(this, sample, library, readGroup);
        }

        /*
        * Normally this would not be overridden by subclasses since the interesting metrics are at
        * the per-unit collector level, but this test collects aggregate values for test purposes
        * so we need to combine those manually.
         */
        public void combine(RecordCountMultiLevelCollector source) {
            // first, combine the per-unit metrics by delegating to the default combine method
            super.combine(source);
            // combine the test-specific stuff
            this.numProcessed = this.getNumProcessed() + source.getNumProcessed();
            Map<String, TotalNumberMetric> combinedUnitsToMetrics = new LinkedHashMap<>(this.getUnitsToMetrics());
            combinedUnitsToMetrics.putAll(source.getUnitsToMetrics());
            combinedUnitsToMetrics.putAll(this.getUnitsToMetrics());
            this.setUnitsToMetrics(combinedUnitsToMetrics);
        }

        @Override
        public RecordCountPerUnitCollector combineUnit(RecordCountPerUnitCollector c1, RecordCountPerUnitCollector c2) {
           return c1.combine(c2);
        }
    }

    public static final Map<MetricAccumulationLevel, Map<String, Integer>> accumulationLevelToPerUnitReads =
            new LinkedHashMap<>();

    static {
        HashMap<String, Integer> curMap = new LinkedHashMap<>();
        curMap.put("__", 19);
        accumulationLevelToPerUnitReads.put(MetricAccumulationLevel.ALL_READS, curMap);

        curMap = new LinkedHashMap<>();
        curMap.put("Ma__", 10);
        curMap.put("Pa__", 9);
        accumulationLevelToPerUnitReads.put(MetricAccumulationLevel.SAMPLE, curMap);

        curMap = new LinkedHashMap<>();
        curMap.put("Ma_whatever_", 10);
        curMap.put("Pa_lib1_",     4);
        curMap.put("Pa_lib2_",     5);
        accumulationLevelToPerUnitReads.put(MetricAccumulationLevel.LIBRARY, curMap);


        curMap = new LinkedHashMap<>();
        curMap.put("Ma_whatever_me",     10);
        curMap.put("Pa_lib1_myself", 4);
        curMap.put("Pa_lib2_i",      3);
        curMap.put("Pa_lib2_i2",     2);
        accumulationLevelToPerUnitReads.put(MetricAccumulationLevel.READ_GROUP, curMap);
    }

    @DataProvider(name = "variedAccumulationLevels")
    public Object[][] variedAccumulationLevels() {
        return new Object[][] {
                {makeSet(MetricAccumulationLevel.ALL_READS)},
                {makeSet(MetricAccumulationLevel.ALL_READS,    MetricAccumulationLevel.SAMPLE)},
                {makeSet(MetricAccumulationLevel.SAMPLE,       MetricAccumulationLevel.LIBRARY)},
                {makeSet(MetricAccumulationLevel.READ_GROUP,   MetricAccumulationLevel.LIBRARY)},
                {makeSet(MetricAccumulationLevel.SAMPLE,       MetricAccumulationLevel.LIBRARY, MetricAccumulationLevel.READ_GROUP)},
                {makeSet(MetricAccumulationLevel.SAMPLE,       MetricAccumulationLevel.LIBRARY, MetricAccumulationLevel.READ_GROUP, MetricAccumulationLevel.ALL_READS)},
        };
    }

    @Test(dataProvider = "variedAccumulationLevels")
    public void multilevelCollectorTest(final Set<MetricAccumulationLevel> accumulationLevels) {
        final SamReader in = SamReaderFactory.makeDefault().open(TESTFILE);

        final RecordCountMultiLevelCollector collector1 = new RecordCountMultiLevelCollector(
                accumulationLevels, in.getFileHeader().getReadGroups());
        final RecordCountMultiLevelCollector collector2 = new RecordCountMultiLevelCollector(
                accumulationLevels, in.getFileHeader().getReadGroups());

        //distribute the reads across the two collectors
        int count = 1;
        for (final SAMRecord rec : in) {
            if (count % 2 == 0) {
                collector1.acceptRecord(rec, null);
            }
            else {
                collector2.acceptRecord(rec, null);
            }
            count++;
        }
        collector1.finish();
        collector2.finish();

        // combine the results into collector1
        collector1.combine(collector2);

        int totalProcessed = 0;
        int totalMetrics = 0;
        for(final MetricAccumulationLevel level : accumulationLevels) {
            final Map<String, Integer> keyToMetrics = accumulationLevelToPerUnitReads.get(level);
            for(final Map.Entry<String, Integer> entry : keyToMetrics.entrySet()) {
                final TotalNumberMetric metric = collector1.getUnitsToMetrics().get(entry.getKey());
                Assert.assertEquals(entry.getValue(), metric.TALLY);
                Assert.assertTrue(metric.FINISHED);
                totalProcessed += metric.TALLY;
                totalMetrics   += 1;
            }
        }

        Assert.assertEquals(collector1.getUnitsToMetrics().size(), totalMetrics);
        Assert.assertEquals(totalProcessed, collector1.getNumProcessed());
        CloserUtil.close(in);
    }
}
