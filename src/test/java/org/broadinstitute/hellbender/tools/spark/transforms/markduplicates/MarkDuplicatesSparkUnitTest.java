package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.engine.spark.SparkTestUtils;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.tools.picard.sam.markduplicates.MarkDuplicatesIntegrationTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;

public class MarkDuplicatesSparkUnitTest extends BaseTest {
    @DataProvider(name = "md")
    public Object[][] loadReads() {
        String dir = MarkDuplicatesIntegrationTest.TEST_DATA_DIR.getAbsolutePath();
        return new Object[][]{
                {dir + "/example.chr1.1-1K.unmarkedDups.noDups.bam", 20, 0},
                {dir + "/example.chr1.1-1K.unmarkedDups.bam", 77, 6}, // was 90... unmapped reads?
                {dir + "/example.chr1.1-1K.markedDups.bam", 77, 6},
        };
    }

    @Test(dataProvider = "md")
    public void markDupesTest(final String input, final long totalExpected, final long dupsExpected) throws IOException {
        JavaSparkContext ctx = SparkTestUtils.getTestContext();

        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        JavaRDD<GATKRead> reads = readSource.getParallelReads(input);
        Assert.assertEquals(reads.count(), totalExpected);

        SAMFileHeader header = ReadsSparkSource.getHeader(ctx, input);
        OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();
        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
        JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.mark(reads, header, finder);

        Assert.assertEquals(markedReads.count(), totalExpected);
        JavaRDD<GATKRead> dupes = markedReads.filter(GATKRead::isDuplicate);

        Assert.assertEquals(dupes.count(), dupsExpected);

        ctx.stop();
    }

}