package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class ReadsSparkSourceUnitTest extends BaseTest {
    @DataProvider(name = "loadReads")
    public Object[][] loadReads() {
        String dir = "src/test/resources/org/broadinstitute/hellbender/tools/BQSR/";
        return new Object[][]{
                {dir + "HiSeq.1mb.1RG.2k_lines.alternate.bam"},
                {dir + "expected.HiSeq.1mb.1RG.2k_lines.bqsr.DIQ.alternate.bam"},
        };
    }

    @Test(dataProvider = "loadReads", groups = "spark")
    public void readsSparkSourceTest(String bam) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddSerialReads = getSerialReads(ctx, bam);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(bam);

        List<GATKRead> serialReads = rddSerialReads.collect();
        List<GATKRead> parallelReads = rddParallelReads.collect();
        Assert.assertEquals(serialReads.size(), parallelReads.size());
    }

    /**
     * Loads Reads using samReaderFactory, then calling ctx.parallelize.
     * @param bam file to load
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getSerialReads(final JavaSparkContext ctx, final String bam) {
        final SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

        ReadsDataSource bam2 = new ReadsDataSource(new File(bam), samReaderFactory);
        bam2.setIntervalsForTraversal(intervals);
        List<GATKRead> records = Lists.newArrayList();
        for ( GATKRead read : bam2 ) {
            records.add(read);
        }
        return ctx.parallelize(records);
    }
}