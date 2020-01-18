package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithReference;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ReadWalkerIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return ExampleReadWalkerWithReference.class.getSimpleName();
    }

    @Test
    public void testManuallySpecifiedIndices() throws IOException {
        final String BAM_PATH = publicTestDir + "org/broadinstitute/hellbender/engine/readIndexTest/";
        final String INDEX_PATH = BAM_PATH + "indices/";
        final File outFile = createTempFile("testManuallySpecifiedIndices", ".txt");
        final File expectedFile = new File(publicTestDir + "org/broadinstitute/hellbender/engine/expected_ReadWalkerIntegrationTest_testManuallySpecifiedIndices.txt");

        final String[] args = new String[] {
            "-I", BAM_PATH + "reads_data_source_test1.bam",
            "-I", BAM_PATH + "reads_data_source_test2.bam",
            "--read-index", INDEX_PATH + "reads_data_source_test1.bam.bai",
            "--read-index", INDEX_PATH + "reads_data_source_test2.bam.bai",
            "-O", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(outFile, expectedFile);
    }

    @Test(expectedExceptions = UserException.class)
    public void testManuallySpecifiedIndicesWrongNumberOfIndices() throws IOException {
        final String BAM_PATH = publicTestDir + "org/broadinstitute/hellbender/engine/readIndexTest/";
        final String INDEX_PATH = BAM_PATH + "indices/";

        final String[] args = new String[] {
                "-I", BAM_PATH + "reads_data_source_test1.bam",
                "-I", BAM_PATH + "reads_data_source_test2.bam",
                "--read-index", INDEX_PATH + "reads_data_source_test1.bam.bai"
        };
        runCommandLine(args);
    }

    @DataProvider
    private Object[][] provideForTestReadStartFilter() {

        final int readsToLinesMultiplier = 3;

        return new Object[][] {
                {
                        Collections.emptyList(),
                        8 * readsToLinesMultiplier
                },
                {
                        Collections.singletonList(new SimpleInterval("chr7", 1,1)),
                        3 * readsToLinesMultiplier
                },
                {
                        Collections.singletonList(new SimpleInterval("chr7", 300,404)),
                        2 * readsToLinesMultiplier
                },
                {
                        Collections.singletonList(new SimpleInterval("chr7", 303,404)),
                        0
                },
                {
                        Collections.singletonList(new SimpleInterval("chr7", 120,404)),
                        2 * readsToLinesMultiplier
                },
                {
                        Collections.singletonList(new SimpleInterval("chr7", 15,17)),
                        1 * readsToLinesMultiplier
                },
                {
                        Collections.singletonList(new SimpleInterval("chr7", 19,30)),
                        1 * readsToLinesMultiplier
                },
                {
                        Arrays.asList(
                                new SimpleInterval("chr7", 1,1),
                                new SimpleInterval("chr7", 21, 21)
                        ),
                        4 * readsToLinesMultiplier
                },
        };
    }

    @Test(dataProvider = "provideForTestReadStartFilter")
    public void testReadStartFilter(final List<SimpleInterval> intervals, final long expectedNumLines) throws IOException {
        final File inFile = new File(getTestDataDir(), "print_reads_withPG.bam");
        final File outFile = GATKBaseTest.createTempFile("testNoConflictRG", ".txt");

        final List<String> argList = new ArrayList<>();

        argList.add("--input");
        argList.add(inFile.getAbsolutePath());
        argList.add("--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD);
        argList.add("--output");
        argList.add(outFile.getAbsolutePath());
        argList.add("--reads-must-start-within-intervals");

        for ( final SimpleInterval interval : intervals ) {
            argList.add("-L");
            argList.add(interval.getContig() + ":" + interval.getStart() + "-" + interval.getEnd());
        }

        runCommandLine(argList);

        // ================

        // The number of lines in the output file will be 3*the number of expected reads:
        final long numLines = Files.lines(IOUtils.getPath(outFile.getAbsolutePath())).count();
        Assert.assertEquals( numLines, expectedNumLines );

    }
}
