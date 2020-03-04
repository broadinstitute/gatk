package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SBIIndex;
import htsjdk.samtools.SBIIndexMerger;
import htsjdk.samtools.SBIIndexWriter;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.tools.AbstractPrintReadsIntegrationTest;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

@Test(groups = "spark")
public final class PrintReadsSparkIntegrationTest extends AbstractPrintReadsIntegrationTest {

    @Test()
    public void testUnSortedRoundTrip() throws Exception {
        // This is a technically incorrectly sam with a header indicating that it is coordinate sorted when it is actually
        // queryname sorted. If the ordering is the same after PrintReadsSpark then it means we aren't automatically sorting the output.
        final File inBam = new File(getTestDataDir(), "print_reads.mismatchedHeader.sam");
        try (ReadsDataSource ds = new ReadsDataSource(inBam.toPath())){
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
