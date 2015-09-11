package org.broadinstitute.hellbender.engine.spark.datasources;


import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;


public class ReadsSparkSinkUnitTest extends BaseTest {
    @DataProvider(name = "loadReads")
    public Object[][] loadReads() {
        String dir = "src/test/resources/org/broadinstitute/hellbender/tools/BQSR/";
        return new Object[][]{
                {dir + "HiSeq.1mb.1RG.2k_lines.alternate.bam", "ReadsSparkSinkUnitTest1.bam"},
                {dir + "expected.HiSeq.1mb.1RG.2k_lines.bqsr.DIQ.alternate.bam", "ReadsSparkSinkUnitTest2.bam"},
        };
    }

    @Test(dataProvider = "loadReads", groups = "spark")
    public void readsSinkTest(String inputBam, String outputFile) throws IOException {
        new File(outputFile).deleteOnExit();
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputBam);
        SAMFileHeader header = ReadsSparkSource.getHeader(ctx, inputBam);

        ReadsSparkSink.writeReads(ctx, outputFile, rddParallelReads, header, true);

        JavaRDD<GATKRead> rddParallelReads2 = readSource.getParallelReads(outputFile);
        final List<GATKRead> writtenReads = rddParallelReads2.collect();

        final ReadCoordinateComparator comparator = new ReadCoordinateComparator(header);
        // Assert that the reads are sorted.
        final int size = writtenReads.size();
        for (int i = 0; i < size-1; ++i) {
            final GATKRead smaller = writtenReads.get(i);
            final GATKRead larger = writtenReads.get(i + 1);
            Assert.assertTrue(comparator.compare(smaller, larger) < 0);
        }
        Assert.assertEquals(rddParallelReads.count(), rddParallelReads2.count());
    }

}