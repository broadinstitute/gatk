package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class SortReadFileSparkIntegrationTest extends CommandLineProgramTest {
    @DataProvider(name="sortbams")
    public Object[][] sortBAMData() {
        return new Object[][] {
                {"count_reads.sam", "count_reads_sorted.sam", null, ".sam", "coordinate"},
                {"count_reads.bam", "count_reads_sorted.bam", null, ".bam", "coordinate"},
                {"count_reads.cram", "count_reads_sorted.cram", "count_reads.fasta", ".bam", "coordinate"},
                {"count_reads.cram", "count_reads_sorted.cram", "count_reads.fasta", ".cram", "coordinate"},
                {"count_reads.bam", "count_reads_sorted.bam", "count_reads.fasta", ".cram", "coordinate"},

                //SortBamSpark is missing SORT_ORDER parameter  https://github.com/broadinstitute/gatk/issues/1260
//                {"count_reads.bam", "count_reads.bam", null, ".bam", "queryname"},
//                {"count_reads.cram", "count_reads.cram", "count_reads.fasta", ".cram", "queryname"},
        };
    }

    @Test(dataProvider="sortbams")
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
        args.add("--input"); args.add(inputFile.getCanonicalPath());
        args.add("--output"); args.add(actualOutputFile.getCanonicalPath());
        if (null != referenceFile) {
            args.add("--R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.add("--numReducers"); args.add("1");

        //https://github.com/broadinstitute/gatk/issues/1260
//        args.add("--SORT_ORDER");
//        args.add(sortOrderName);

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.samsEqualStringent(actualOutputFile, expectedOutputFile, ValidationStringency.DEFAULT_STRINGENCY, referenceFile);
    }

    @Test
    public void test() throws Exception {
        final File unsortedBam = new File(getTestDataDir(), "count_reads.bam");
        final File sortedBam = new File(getTestDataDir(), "count_reads_sorted.bam");
        final File outputBam = createTempFile("sort_bam_spark", ".bam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(unsortedBam.getCanonicalPath());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME); args.add(outputBam.getCanonicalPath());
        args.add("--numReducers"); args.add("1");

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(outputBam, sortedBam);
    }

}