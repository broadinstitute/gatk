package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Integration tests for {@link CombineReadCounts}.  The main functionality is tested in
 * {@link ReadCountCollectionUnitTest}.  This class merely tests that the behavior of {@link ReadCountCollection}
 * occurs.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class NormalizeBySampleDepthIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File INPUT_WITH_INTERVALS = new File(TEST_DIR, "exome-read-counts.output");
    private static final File INPUT_WITHOUT_INTERVALS = new File(TEST_DIR, "exome-read-counts-no-intervals.output");

    @Test(expectedExceptions = UserException.class)
    public void testWeightingNoIntervals() {
        final File outputFile = createTempFile("test", ".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_WITHOUT_INTERVALS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + NormalizeBySampleDepth.WEIGHTED_AVERAGE_SHORT_NAME
        };
        runCommandLine(arguments);
    }

    @Test
    public void testNoWeighting() {
        final File outputFile = createTempFile("test", ".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_WITHOUT_INTERVALS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        try {
            final ReadCountCollection input = ReadCountCollectionUtils.parse(INPUT_WITHOUT_INTERVALS);
            final ReadCountCollection output = ReadCountCollectionUtils.parse(outputFile);
            final ReadCountCollection expectedOutput = input.normalizeByColumnAverages(false);
            Assert.assertEquals(expectedOutput.targets(), output.targets());
            Assert.assertEquals(expectedOutput.columnNames(), output.columnNames());
            Assert.assertEquals(expectedOutput.counts().subtract(output.counts()).getNorm(), 0, 1e-8);
        } catch (final Exception e) {
        }
    }

    @Test
    public void testWeighting() {
        final File outputFile = createTempFile("test", ".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_WITH_INTERVALS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + NormalizeBySampleDepth.WEIGHTED_AVERAGE_SHORT_NAME
        };
        runCommandLine(arguments);
        try {
            final ReadCountCollection input = ReadCountCollectionUtils.parse(INPUT_WITH_INTERVALS);
            final ReadCountCollection output = ReadCountCollectionUtils.parse(outputFile);
            final ReadCountCollection expectedOutput = input.normalizeByColumnAverages(true);
            Assert.assertEquals(expectedOutput.targets(), output.targets());
            Assert.assertEquals(expectedOutput.columnNames(), output.columnNames());
            Assert.assertEquals(expectedOutput.counts().subtract(output.counts()).getNorm(), 0, 1e-8);
        } catch (final Exception e) {
        }
    }
}