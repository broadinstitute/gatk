package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public final class PrintReadsIntegrationTest extends AbstractPrintReadsIntegrationTest {

    @Test(dataProvider="testingData")
    public void testFileToFileWithMD5(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, true);
    }

    @Test
    public void testNoConflictPG() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "print_reads_withPG.sam");
        final File outFile = GATKBaseTest.createTempFile("testNoConflictRG", ".sam");
        final String[] args = new String[] {
                "--input" , inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //Make sure contents are the same
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has GATK PrintReads not NOT GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));

        //output has both GATK PrintReads and GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));
    }

    @DataProvider
    public Object[][] getHttpPaths(){
        final String bam = "gs://hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam";
        final String bai = "gs://hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam.bai";
        final String cram = "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram";
        final String crai = "gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram.bai";
        final List<SimpleInterval> bamIntervals = Arrays.asList(new SimpleInterval("3",1_000_000, 1_000_001),
               new SimpleInterval("3", 1_000_003, 1_000_100),
                new SimpleInterval("20", 1_099_000, 1_100_000));

        final List<SimpleInterval> cramIntervals = Arrays.asList(new SimpleInterval("20", 9_999_902, 10_000_000));
        return new Object[][]{
                {BucketUtils.bucketPathToPublicHttpUrl(bam), BucketUtils.bucketPathToPublicHttpUrl(bai), bam, bai, bamIntervals, 528L},
                {BucketUtils.createSignedUrlToGcsObject(bam, 1L), BucketUtils.createSignedUrlToGcsObject(bai, 1L), bam, bai, bamIntervals, 528L},
                {BucketUtils.bucketPathToPublicHttpUrl(cram), BucketUtils.bucketPathToPublicHttpUrl(crai), cram, crai, cramIntervals, 112L},
                {BucketUtils.createSignedUrlToGcsObject(cram, 1L), BucketUtils.createSignedUrlToGcsObject(crai, 1L), cram, crai, cramIntervals, 112L}

        };
    }

    @Test(groups = {"cloud", "bucket"}, dataProvider = "getHttpPaths")
    public void testHttpPaths(String reads, String index, String nonHttpReads, String nonHttpIndex, List<SimpleInterval> intervals, long expectedNumberOfReads) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("out", ".bam");
        // this test reads tiny amounts of data from multiple places, if you don't set the prefetcher to a lower number
        // it loads large amounts of data that slows the test down significantly for no good reason
        args.addInput(reads)
                .add(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, 1)
                .add(StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, 1)
                .add("read-index", index)
                .addReference(GATKBaseTest.b37Reference)
                .addOutput(out);
        intervals.forEach(args::addInterval);
        runCommandLine(args);

        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        final File out2 = createTempFile("out", ".bam");
        args2.addInput(nonHttpReads)
                .add("read-index", nonHttpIndex)
                .add(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, 1)
                .add(StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, 1)
                .addReference(GATKBaseTest.b37Reference)
                .addOutput(out2);
        intervals.forEach(args2::addInterval);
        runCommandLine(args2);

        try(final ReadsDataSource reader = new ReadsPathDataSource(IOUtils.toGATKPath(out))){
            final long count = Utils.stream(reader).count();
            Assert.assertEquals( count, expectedNumberOfReads);
        }

        SamAssertionUtils.assertEqualBamFiles(out, out2, false, ValidationStringency.DEFAULT_STRINGENCY);
    }

}