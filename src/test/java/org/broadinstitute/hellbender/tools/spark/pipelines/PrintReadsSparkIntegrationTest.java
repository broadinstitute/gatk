package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SBIIndex;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.tools.AbstractPrintReadsIntegrationTest;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

@Test(groups = "spark")
public final class PrintReadsSparkIntegrationTest extends AbstractPrintReadsIntegrationTest {

    @Test()
    public void testUnSortedRoundTrip() throws Exception {
        // This is a technically incorrectly sam with a header indicating that it is coordinate sorted when it is actually
        // queryname sorted. If the ordering is the same after PrintReadsSpark then it means we aren't automatically sorting the output.
        final File inBam = new File(getTestDataDir(), "print_reads.mismatchedHeader.sam");
        try (ReadsDataSource ds = new ReadsPathDataSource(IOUtils.toGATKPath(inBam))){
            Assert.assertEquals(ds.getHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        }
        final File outBam = GATKBaseTest.createTempFile("print_reads", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.addRaw(inBam.getCanonicalPath());
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.addRaw(outBam.getCanonicalPath());

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outBam, inBam);
    }

    /**
     * Test GCS access using Spark NIO.
     *
     * For this to work, the settings in src/main/resources/core-site.xml must be correct,
     * and the project name and credential file it points to must be present.
     */
    @Test(dataProvider = "gcsTestingData", groups = "bucket")
    public void testGCSInputsAndOutputsWithSparkNio(final String gcsInput, final String outputExtension,
                                        final boolean outputToGCS, final File expectedOutput) throws IOException {
        final String gcsInputPath = getGCPTestInputPath() + gcsInput;
        final String outputPrefix = outputToGCS ? getGCPTestStaging() : "testGCSInputsAndOutputs";
        final String outputPath = BucketUtils.getTempFilePath(outputPrefix, outputExtension);

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.add("input", gcsInputPath)
                .add("output", outputPath)
                .add(GATKSparkTool.USE_NIO, true);
        runCommandLine(argBuilder);

        SamAssertionUtils.assertSamsEqual(IOUtils.getPath(outputPath), IOUtils.getPath(gcsInputPath), null);
    }

    @Test(groups = "spark")
    public void testReadAndWriteCRAMAndReferenceOnHDFS() throws Exception {
        final GATKPath testCram = new GATKPath(getTestFile("count_reads.cram").getAbsolutePath());
        final GATKPath testRef = new GATKPath(getTestFile("count_reads.fasta").getAbsolutePath());
        final GATKPath testRefDict = new GATKPath(getTestFile("count_reads.dict").getAbsolutePath());
        final GATKPath testRefIndex = new GATKPath(getTestFile("count_reads.fasta.fai").getAbsolutePath());

        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {
            final org.apache.hadoop.fs.Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);

            // copy all the test files to HDFS
            final org.apache.hadoop.fs.Path cramHDFSPath = new org.apache.hadoop.fs.Path(workingDirectory,"count_reads.cram");
            final org.apache.hadoop.fs.Path refHDFSPath = new org.apache.hadoop.fs.Path(workingDirectory, "count_reads.fasta");
            final org.apache.hadoop.fs.Path refHDFSDictPath = new org.apache.hadoop.fs.Path(workingDirectory,"count_reads.dict");
            final org.apache.hadoop.fs.Path refHDFSIndexPath = new org.apache.hadoop.fs.Path(workingDirectory, "count_reads.fasta.fai");
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(testCram.getURI()), cramHDFSPath);
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(testRef.getURI()), refHDFSPath);
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(testRefDict.getURI()), refHDFSDictPath);
            cluster.getFileSystem().copyFromLocalFile(new org.apache.hadoop.fs.Path(testRefIndex.getURI()), refHDFSIndexPath);

            // run PrintReadsSpark and print the contents of the HDFS cram test file to an output HDFS cram
            final GATKPath outputHDFSPath = new GATKPath(workingDirectory + "testCramOnHDFSOut.cram");
            final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
            argBuilder.addInput(cramHDFSPath.toUri().toString())
                    .addOutput( outputHDFSPath.getURI().toString())
                    .addReference(refHDFSPath.toUri().toString());
            runCommandLine(argBuilder);

            // compare the original test file with the HDFS cram written by PrintReadsSpark
            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
            final ReadsSparkSource readsSparkSource = new ReadsSparkSource(ctx);
            final List<GATKRead> localReads = readsSparkSource.getParallelReads(testCram, testRef).collect();
            final List<GATKRead> hdfsReads = readsSparkSource.getParallelReads(outputHDFSPath, testRef).collect();

            Assert.assertFalse(localReads.isEmpty());
            Assert.assertEquals(localReads, hdfsReads);
        });
    }

    @DataProvider
    public Object[][] getGranularityLevels(){
        return new Object[][] {{1L}, {100L}};
    }
    @Test(dataProvider = "getGranularityLevels")
    public void testSBIIndexGranularityArgument(long granularity) throws IOException {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        File outBam = createTempFile("prefix","out.bam");
        Path outBamPath = outBam.toPath();
        Path sbiPath = outBamPath.resolveSibling(outBamPath.getFileName().toString() +  ".sbi");
        args.addInput(getTestFile( "abam.bam"))
                .addOutput(outBam)
                .add((GATKSparkTool.SPLITTING_INDEX_GRANULARITY), granularity);

        runCommandLine(args);
        SBIIndex load = SBIIndex.load(sbiPath);
        Assert.assertEquals(load.getGranularity(), granularity);

    }
}
