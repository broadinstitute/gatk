package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

public class PrintReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final String myTestDir = toolsTestDir + "sv/";

    @Test
    public void testRD() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.INTERVALS_SHORT_NAME, "chr1:10414291-10416290");
        argsBuilder.add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, myTestDir + "chr1.dict");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, myTestDir + "test.HG00096.rd.txt.gz");
        argsBuilder.add(PrintReadCounts.OUTPUT_PREFIX_ARGNAME, myTestDir + "testCounts.");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(), Collections.emptyList());
        testSpec.executeTest("testRD", this);
        checkOutput(myTestDir + "testCounts.HG00096.counts.tsv",
                myTestDir + "expect.HG00096.counts.tsv");
    }

    @Test
    public void testCounts() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.INTERVALS_SHORT_NAME, "chr1:10414291-10416290");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, myTestDir + "test.HG00096.counts.tsv.gz");
        argsBuilder.add(PrintReadCounts.OUTPUT_PREFIX_ARGNAME, myTestDir + "testCounts.");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(), Collections.emptyList());
        testSpec.executeTest("testCounts", this);
        checkOutput(myTestDir + "testCounts.HG00096.counts.tsv",
                myTestDir + "expect.HG00096.counts.tsv");
    }

    @Test
    public void testBCI() throws IOException {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.INTERVALS_SHORT_NAME, "chr1:10414291-10416290");
        argsBuilder.add(StandardArgumentDefinitions.FEATURE_SHORT_NAME, myTestDir + "merged.rd.bci");
        argsBuilder.add(PrintReadCounts.OUTPUT_PREFIX_ARGNAME, myTestDir + "testCounts.");
        final IntegrationTestSpec testSpec =
                new IntegrationTestSpec(argsBuilder.getString(), Collections.emptyList());
        testSpec.executeTest("testBCI", this);
        checkOutput(myTestDir + "testCounts.HG00096.counts.tsv",
                myTestDir + "expect.HG00096.counts.tsv");
        checkOutput(myTestDir + "testCounts.HG00129.counts.tsv",
                myTestDir + "expect.HG00129.counts.tsv");
        checkOutput(myTestDir + "testCounts.HG00140.counts.tsv",
                myTestDir + "expect.HG00140.counts.tsv");
    }

    private void checkOutput( final String actualPath, final String expectedPath ) throws IOException {
        final File actual = new File(actualPath);
        try {
            IntegrationTestSpec.assertEqualTextFiles(actual, new File(expectedPath));
        } finally {
            actual.delete();
        }
    }
}
