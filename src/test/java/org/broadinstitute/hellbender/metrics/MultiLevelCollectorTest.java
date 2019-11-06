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

public final class MultiLevelCollectorTest {

    public static File TESTFILE = new File("src/test/resources/org/broadinstitute/hellbender/metrics/test.sam");

    public String noneOrStr(final String str) {
        final String out;
        if(str == null) {
            out = "";
        } else {
            out = str;
        }
        return out;
    }

    class TestArg {
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
    class TotalNumberMetric extends MultiLevelMetrics implements Serializable {
        private static final long serialVersionUID = 1L;

        /** The number of these encountered **/
        public Integer TALLY = 0;
        public boolean FINISHED = false;
    }

    class RecordCountMultiLevelCollector extends MultiLevelCollector<TotalNumberMetric, Integer, TestArg> {
        private static final long serialVersionUID = 1L;

        public RecordCountMultiLevelCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords) {
            setup(accumulationLevels, samRgRecords);
        }

        //The number of times records were accepted by a RecordCountPerUnitCollectors (note since the same
        //samRecord might be aggregated by multiple PerUnit collectors, this may be greater than the number of
        //records in the file
        private int numProcessed = 0;

        public int getNumProcessed() {
            return numProcessed;
        }

        private final Map<String, TotalNumberMetric> unitsToMetrics = new LinkedHashMap<>();

        public Map<String, TotalNumberMetric> getUnitsToMetrics() {
            return unitsToMetrics;
        }

        @Override
        protected TestArg makeArg(final SAMRecord samRec, final ReferenceSequence refSeq) {
            return new TestArg(samRec, refSeq);
        }

        @Override
        protected PerUnitMetricCollector<TotalNumberMetric, Integer, TestArg> makeChildCollector(final String sample, final String library, final String readGroup) {
            return new RecordCountPerUnitCollector(sample, library, readGroup);
        }

        private class RecordCountPerUnitCollector implements PerUnitMetricCollector<TotalNumberMetric, Integer, TestArg>{
            private static final long serialVersionUID = 1L;
            private final TotalNumberMetric metric;

            public RecordCountPerUnitCollector(final String sample, final String library, final String readGroup) {
                metric = new TotalNumberMetric();
                metric.SAMPLE     = sample;
                metric.LIBRARY    = library;
                metric.READ_GROUP = readGroup;
                unitsToMetrics.put(noneOrStr(sample) + "_" + noneOrStr(library) + "_" + noneOrStr(readGroup), metric);
            }

            @Override
            public void acceptRecord(final TestArg args) {
                numProcessed += 1;
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
        }
    }

    public static final Map<MetricAccumulationLevel, Map<String, Integer>> accumulationLevelToPerUnitReads = new LinkedHashMap<>();
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
        final RecordCountMultiLevelCollector collector = new RecordCountMultiLevelCollector(accumulationLevels, in.getFileHeader().getReadGroups());

        for (final SAMRecord rec : in) {
            collector.acceptRecord(rec, null);
        }

        collector.finish();

        int totalProcessed = 0;
        int totalMetrics = 0;
        for(final MetricAccumulationLevel level : accumulationLevels) {
            final Map<String, Integer> keyToMetrics = accumulationLevelToPerUnitReads.get(level);
            for(final Map.Entry<String, Integer> entry : keyToMetrics.entrySet()) {
                final TotalNumberMetric metric = collector.getUnitsToMetrics().get(entry.getKey());
                Assert.assertEquals(entry.getValue(), metric.TALLY);
                Assert.assertTrue(metric.FINISHED);
                totalProcessed += metric.TALLY;
                totalMetrics   += 1;
            }
        }

        Assert.assertEquals(collector.getUnitsToMetrics().size(), totalMetrics);
        Assert.assertEquals(totalProcessed, collector.getNumProcessed());
        CloserUtil.close(in);
    }
}
