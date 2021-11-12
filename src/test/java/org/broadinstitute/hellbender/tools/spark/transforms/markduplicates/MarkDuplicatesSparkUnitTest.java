//package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;
//
//import htsjdk.samtools.SAMFileHeader;
//import org.apache.spark.api.java.JavaRDD;
//import org.apache.spark.api.java.JavaSparkContext;
//import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
//import org.broadinstitute.hellbender.engine.GATKPath;
//import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
//import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
//import org.broadinstitute.hellbender.utils.read.GATKRead;
//import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
//import picard.sam.markduplicates.MarkDuplicates;
//import picard.sam.markduplicates.util.OpticalDuplicateFinder;
//import org.broadinstitute.hellbender.GATKBaseTest;
//import org.testng.Assert;
//import org.testng.annotations.DataProvider;
//import org.testng.annotations.Test;
//
//import java.io.File;
//
//public class MarkDuplicatesSparkUnitTest extends GATKBaseTest {
//    @DataProvider(name = "md")
//    public Object[][] loadReads() {
//        String dir = new File(toolsTestDir, "walkers/MarkDuplicatesGATK/").getAbsolutePath();
//        return new Object[][]{
//                {dir + "/example.chr1.1-1K.unmarkedDups.noDups.bam", 20, 0},
//                {dir + "/example.chr1.1-1K.unmarkedDups.bam", 90, 6},
//                {dir + "/example.chr1.1-1K.markedDups.bam", 90, 6},
//        };
//    }
//
//    @Test(dataProvider = "md", groups = "spark")
//    public void markDupesTest(final String input, final long totalExpected, final long dupsExpected) {
//        final GATKPath inputPathSpec = new GATKPath(input);
//        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
//
//        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
//        JavaRDD<GATKRead> reads = readSource.getParallelReads(inputPathSpec, null);
//        Assert.assertEquals(reads.count(), totalExpected);
//
//        SAMFileHeader header = readSource.getHeader(inputPathSpec, null);
//        OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();
//        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
//                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
//        JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.mark(reads, header, MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES, finder, 1, false, MarkDuplicates.DuplicateTaggingPolicy.DontTag);
//
//        Assert.assertEquals(markedReads.count(), totalExpected);
//        JavaRDD<GATKRead> dupes = markedReads.filter(GATKRead::isDuplicate);
//
//        Assert.assertEquals(dupes.count(), dupsExpected);
//    }
//
//}
