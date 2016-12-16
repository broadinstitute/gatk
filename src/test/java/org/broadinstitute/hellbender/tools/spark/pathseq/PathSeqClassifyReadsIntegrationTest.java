package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PathSeqClassifyReadsIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return PathSeqClassifyReadsSpark.class.getSimpleName();
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void test() throws Exception {
        final File expectedFile = getTestFile("expected.txt");
        final File inputFile = getTestFile("alignment.bam");
        final File taxFile = getTestFile("tax.db");
        final File output = createTempFile("test", ".txt");
        final File bamOut = createTempFile("output", ".bam");
        if (!output.delete() || !bamOut.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("pairedInput", inputFile.getAbsolutePath());
        args.addOutput(bamOut);
        args.addFileArgument("taxonomicDatabasePath", taxFile);
        args.addFileArgument("scoresOutputPath", output);
        this.runCommandLine(args.getArgsArray());

        final byte[] input_expected = FileUtils.readFileToByteArray(expectedFile);
        final byte[] input_test = FileUtils.readFileToByteArray(output);

        Assert.assertEquals(input_test, input_expected);
    }

    @SuppressWarnings("unchecked")
    @Test(groups = "spark")
    public void testUnpaired() throws Exception {
        final File expectedFile = getTestFile("expected_unpaired.txt");
        final File inputFile = getTestFile("alignment.bam");
        final File taxFile = getTestFile("tax.db");
        final File output = createTempFile("test", ".txt");
        final File warnings = createTempFile("warnings", ".txt");
        if (!output.delete() || !warnings.delete()) {
            Assert.fail();
        }
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument("unpairedInput", inputFile.getAbsolutePath());
        args.addFileArgument("taxonomicDatabasePath", taxFile);
        args.addFileArgument("scoresOutputPath", output);
        args.addFileArgument("warningsFile",warnings);
        this.runCommandLine(args.getArgsArray());

        final byte[] input_expected = FileUtils.readFileToByteArray(expectedFile);
        final byte[] input_test = FileUtils.readFileToByteArray(output);

        Assert.assertEquals(input_test, input_expected);
    }

}