package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Created by David Benjamin on 2/29/16.
 */
public class CalculateCoverageZScoresIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File INPUT_COUNTS = new File(TEST_DIR, "exome-read-counts.output");

    @Test
    public void testZScores() throws IOException {
        final File outputFile = createTempFile("test", ".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_COUNTS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        final ReadCountCollection input = ReadCountCollectionUtils.parse(INPUT_COUNTS);
        final ReadCountCollection output = ReadCountCollectionUtils.parse(outputFile);
        final ReadCountCollection expectedOutput = input.zScoreCounts();
        Assert.assertEquals(expectedOutput.targets(), output.targets());
        Assert.assertEquals(expectedOutput.columnNames(), output.columnNames());
        Assert.assertEquals(expectedOutput.counts().subtract(output.counts()).getNorm(), 0, 1e-8);
    }
}

