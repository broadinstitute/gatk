package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class SortSamSparkIntegrationTest extends CommandLineProgramTest {
    @DataProvider(name="sortbams")
    public Object[][] sortBAMData() {
        return new Object[][] {
                {"count_reads.sam", "count_reads_sorted.sam", null, ".sam", "coordinate"},
                {"count_reads.bam", "count_reads_sorted.bam", null, ".bam", "coordinate"},
                {"count_reads.cram", "count_reads_sorted.cram", "count_reads.fasta", ".bam", "coordinate"},
                {"count_reads.cram", "count_reads_sorted.cram", "count_reads.fasta", ".cram", "coordinate"},
                {"count_reads.bam", "count_reads_sorted.bam", "count_reads.fasta", ".cram", "coordinate"},

                {"count_reads.bam", "count_reads.bam", null, ".bam", "queryname"},
                {"count_reads.cram", "count_reads.cram", "count_reads.fasta", ".cram", "queryname"},
        };
    }

    @Test(dataProvider="sortbams", groups="spark")
    public void testSortBAMs(
            final String inputFileName,
            final String expectedOutputFileName,
            final String referenceFileName,
            final String outputExtension,
            final String sortOrderName) throws Exception {
        final File inputFile = new File(getTestDataDir(), inputFileName);
        final File expectedOutputFile = new File(getTestDataDir(), expectedOutputFileName);
        final File actualOutputFile = createTempFile("sort_sam", outputExtension);
        File referenceFile = null == referenceFileName ? null : new File(getTestDataDir(), referenceFileName);
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(inputFile);
        args.addOutput(actualOutputFile);
        if (null != referenceFile) {
            args.addReference(referenceFile);
        }
        args.addArgument(GATKSparkTool.NUM_REDUCERS_LONG_NAME, "1");
        args.addArgument(SortSamSpark.SORT_ORDER_LONG_NAME, sortOrderName);

        this.runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(actualOutputFile, expectedOutputFile, ValidationStringency.DEFAULT_STRINGENCY, referenceFile);
    }

    @Test(groups = "spark")
    public void test() throws Exception {
        final File unsortedBam = new File(getTestDataDir(), "count_reads.bam");
        final File sortedBam = new File(getTestDataDir(), "count_reads_sorted.bam");
        final File outputBam = createTempFile("sort_bam_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(unsortedBam);
        args.addOutput(outputBam);
        args.addArgument(GATKSparkTool.NUM_REDUCERS_LONG_NAME, "1");

        this.runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outputBam, sortedBam);
    }


    @DataProvider
    public Object[][] getInvalidSortOrders(){
        return new Object[][]{
                {SAMFileHeader.SortOrder.unknown},
                {SAMFileHeader.SortOrder.unsorted},
                {SAMFileHeader.SortOrder.duplicate}
        };
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class, dataProvider = "getInvalidSortOrders")
    public void testBadSortOrders(SAMFileHeader.SortOrder badOrder){
        final File unsortedBam = new File(getTestDataDir(), "count_reads.bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(unsortedBam);
        args.addOutput(createTempFile("sort_bam_spark", ".bam"));
        args.addArgument(SortSamSpark.SORT_ORDER_LONG_NAME, badOrder.toString());

        this.runCommandLine(args);
    }
}
