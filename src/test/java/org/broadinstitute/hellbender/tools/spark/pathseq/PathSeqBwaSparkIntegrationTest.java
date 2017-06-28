package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqBwaSparkIntegrationTest extends CommandLineProgramTest {

    private final String IMAGE_PATH = "src/test/resources/" + PathSeqBwaSpark.class.getPackage().getName().replace(".", "/") + "/hg19mini.fasta.bwa_image";
    private final String REF_PATH = "src/test/resources/hg19mini.fasta";

    @Override
    public String getTestedClassName() {
        return PathSeqBwaSpark.class.getSimpleName();
    }

    private void testFilterTool(final ArgumentsBuilder args, final String inputFilename, final String expectedFilename) throws Exception {
        final File inputFile = getTestFile(inputFilename);
        final File tmpDir = createTempDir("tmp");
        if (!tmpDir.delete()) {
            Assert.fail();
        }
        final File outputFile = new File(tmpDir.getAbsolutePath(), "output");
        final String outputBasePath = outputFile.getAbsolutePath();
        args.addArgument("pairedInput", inputFile.getAbsolutePath());
        args.addArgument("output", outputBasePath);
        args.addArgument("bwamemIndexImage", IMAGE_PATH);
        args.addArgument("referencePath", REF_PATH);

        this.runCommandLine(args.getArgsArray());

        SamAssertionUtils.assertSamsEqual(new File(outputBasePath + ".paired.bam"), getTestFile(expectedFilename), ValidationStringency.LENIENT);
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testAlignBasic() throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        testFilterTool(args, "basic_input.bam", "basic_output.paired.bam");
    }

}
