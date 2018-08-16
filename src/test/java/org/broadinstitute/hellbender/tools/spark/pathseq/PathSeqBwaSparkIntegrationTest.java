package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqBwaSparkIntegrationTest extends CommandLineProgramTest {

    private final String IMAGE_PATH = publicTestDir + "hg19mini.fasta.img";
    private final String REF_PATH = publicTestDir + "hg19mini.fasta";

    @DataProvider(name = "pathseqBwaTestData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"basic_input.bam", "basic_output.paired.bam"}
        };
    }

    @Override
    public String getTestedClassName() {
        return PathSeqBwaSpark.class.getSimpleName();
    }

    @Test(dataProvider = "pathseqBwaTestData")
    private void testBwaTool(final String inputBamFilename, final String expectedBamFilename) throws Exception {
        final File inputBamFile = getTestFile(inputBamFilename);
        final File pairedOutputBamFile = createTempFile("paired_output", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addFileArgument(PathSeqBwaSpark.PAIRED_INPUT_LONG_NAME, inputBamFile);
        args.addFileArgument(PathSeqBwaSpark.PAIRED_OUTPUT_LONG_NAME, pairedOutputBamFile);
        args.addArgument(PSBwaArgumentCollection.MICROBE_BWA_IMAGE_LONG_NAME, IMAGE_PATH);
        args.addArgument(PSBwaArgumentCollection.MICROBE_FASTA_LONG_NAME, REF_PATH);
        this.runCommandLine(args.getArgsArray());
        SamAssertionUtils.assertSamsEqual(pairedOutputBamFile, getTestFile(expectedBamFilename), ValidationStringency.LENIENT);
    }
}
