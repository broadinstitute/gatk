package org.broadinstitute.hellbender.tools;

import com.google.cloud.storage.BlobInfo;
import com.google.cloud.storage.HttpMethod;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

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
        final List<SimpleInterval> largeFileIntervals = Arrays.asList(new SimpleInterval("3",1_000_000, 1_000_001),
               new SimpleInterval("3", 1_000_003, 1_000_100),
                new SimpleInterval("20", 1_099_000, 1_100_000));

        final List<SimpleInterval> cramIntervals = Arrays.asList(new SimpleInterval("20", 1_099_000, 1_100_000));
        return new Object[][]{
                {bucketPathToPublicHttpUrl(bam), bucketPathToPublicHttpUrl(bai), bam, bai, largeFileIntervals},
                {toSignedUrl(bam), toSignedUrl(bai), bam, bai, largeFileIntervals},
                {bucketPathToPublicHttpUrl(cram), bucketPathToPublicHttpUrl(crai), cram, crai, cramIntervals},
                {toSignedUrl(cram), toSignedUrl(crai), cram, crai, cramIntervals}

        };
    }

    @Test(groups = {"cloud", "bucket"}, dataProvider = "getHttpPaths")
    public void testHttpPaths(String reads, String index, String nonHttpReads, String nonHttpIndex, List<SimpleInterval> intervals) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("out", ".bam");
        System.out.println("reading http");
        args.addInput(reads)
                .add(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, 1)
                .add(StandardArgumentDefinitions.CLOUD_INDEX_PREFETCH_BUFFER_LONG_NAME, 1)
                .add("read-index", index)
                .addReference(GATKBaseTest.b37Reference)
                .addOutput(out);
        intervals.forEach(args::addInterval);
        runCommandLine(args);

        System.out.println("reading gs");
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

        SamAssertionUtils.assertEqualBamFiles(out, out2, false, ValidationStringency.DEFAULT_STRINGENCY);
    }

    public static String toSignedUrl(String path) {
        final Storage storage = StorageOptions.getDefaultInstance().getService();
        final BlobInfo info = BlobInfo.newBuilder(BucketUtils.getBucket(path), BucketUtils.getPathWithoutBucket(path)).build();
        return storage.signUrl(info, 1L, TimeUnit.HOURS,
                Storage.SignUrlOption.httpMethod(HttpMethod.GET)
        ).toString();
    }

    public static String bucketPathToPublicHttpUrl(String path){
        return String.format("https://storage.googleapis.com/%s/%s", BucketUtils.getBucket(path), BucketUtils.getPathWithoutBucket(path));
    }

}