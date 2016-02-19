package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class SortSamIntegrationTest extends CommandLineProgramTest {

    @DataProvider(name="sortbams")
    public Object[][] sortBAMData() {
        return new Object[][] {
                {"count_reads.bam", "count_reads_sorted.bam", null, ".bam", "coordinate"},
                {"count_reads.bam", "count_reads_sorted.bam", "count_reads.fasta", ".cram", "coordinate"},
                {"count_reads.bam", "count_reads.bam", null, ".bam", "queryname"},
                {"count_reads.cram", "count_reads_sorted.cram", "count_reads.fasta", ".cram", "coordinate"},
                {"count_reads.cram", "count_reads_sorted.cram", "count_reads.fasta", ".bam", "coordinate"},
                {"count_reads.cram", "count_reads.cram", "count_reads.fasta", ".cram", "queryname"}
        };
    }

    @Test(dataProvider="sortbams")
    public void testSortBAMs(
            final String inputFileName,
            final String expectedOutputFileName,
            final String referenceFileName,
            final String outputExtension,
            final String sortOrderName) throws Exception
    {
        final File inputBam = new File(getTestDataDir(), inputFileName);
        final File expectedBam = new File(getTestDataDir(), expectedOutputFileName);
        final File outputBam = createTempFile("sort_sam", outputExtension);
        File referenceFile = null == referenceFileName ? null : new File(getTestDataDir(), referenceFileName);
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input"); args.add(inputBam.getCanonicalPath());
        args.add("--output"); args.add(outputBam.getCanonicalPath());
        if (null != referenceFile) {
            args.add("--R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.add("--SORT_ORDER");
        args.add(sortOrderName);

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.samsEqualStringent(expectedBam, outputBam, ValidationStringency.DEFAULT_STRINGENCY, referenceFile);
    }
}

