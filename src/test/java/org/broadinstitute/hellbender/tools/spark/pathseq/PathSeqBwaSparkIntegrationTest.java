package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqBwaSparkIntegrationTest extends CommandLineProgramTest {

    private final String IMAGE_PATH = "src/test/resources/" + PathSeqBwaSpark.class.getPackage().getName().replace(".", "/") + "/hg19mini.fasta.img";
    private final String REF_PATH = "src/test/resources/hg19mini.fasta";

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
        final File outputTempDir = createTempDir("tmp_pathseqBwaTest");
        final File outputBamFile = new File(outputTempDir.getAbsolutePath(), "output");
        final String outputBasePath = outputBamFile.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("pairedInput", inputBamFile.getAbsolutePath());
        args.addArgument("output", outputBasePath);
        args.addArgument("pathogenBwaImage", IMAGE_PATH);
        args.addArgument("pathogenFasta", REF_PATH);
        this.runCommandLine(args.getArgsArray());
        SamAssertionUtils.assertSamsEqual(new File(outputBasePath + ".paired.bam"), getTestFile(expectedBamFilename), ValidationStringency.LENIENT);
    }
}
