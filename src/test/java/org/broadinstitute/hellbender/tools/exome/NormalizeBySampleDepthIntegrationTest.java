package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration tests for {@link CombineReadCounts}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class NormalizeBySampleDepthIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File INPUT_WITH_INTERVALS = new File(TEST_DIR,"exome-read-counts.output");
    private static final File INPUT_WITHOUT_INTERVALS = new File(TEST_DIR,"exome-read-counts-no-intervals.output");

    @Test(expectedExceptions = UserException.class)
    public void testWeightingNoIntervals() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_WITHOUT_INTERVALS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + NormalizeBySampleDepth.WEIGHTED_AVERAGE_SHORT_NAME
        };
        runCommandLine(arguments);
    }

    @Test
    public void testNoWeighting() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_WITHOUT_INTERVALS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        try {
            final ReadCountCollection input = ReadCountCollectionUtils.parse(INPUT_WITHOUT_INTERVALS);
            final ReadCountCollection output = ReadCountCollectionUtils.parse(outputFile);
            Assert.assertEquals(input.targets(), output.targets());
            Assert.assertEquals(input.columnNames(), output.columnNames());

            final RealMatrix inputCounts = input.counts();
            final RealMatrix outputCounts = output.counts();

            Assert.assertEquals(inputCounts.getRowDimension(), outputCounts.getRowDimension());
            Assert.assertEquals(inputCounts.getColumnDimension(), outputCounts.getColumnDimension());

            for (int col = 0; col < inputCounts.getColumnDimension(); col++) {
                final RealVector inputColumn = inputCounts.getColumnVector(col);
                final RealVector outputColumn = outputCounts.getColumnVector(col);
                Assert.assertEquals(inputColumn.cosine(outputColumn), 1, 1e-8);
                final double inputAverage = inputColumn.getL1Norm()/inputColumn.getDimension();
                Assert.assertEquals(outputColumn.mapMultiply(inputAverage).getL1Distance(inputColumn), 0, 1e-10);
            }
        } catch (final IOException e) {
            //throw new Exception(e);
        }
    }

    @Test
    public void testWeighting() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_WITH_INTERVALS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + NormalizeBySampleDepth.WEIGHTED_AVERAGE_SHORT_NAME
        };
        runCommandLine(arguments);

        try {
            final ReadCountCollection input = ReadCountCollectionUtils.parse(INPUT_WITH_INTERVALS);
            final ReadCountCollection output = ReadCountCollectionUtils.parse(outputFile);
            Assert.assertEquals(input.targets(), output.targets());
            Assert.assertEquals(input.columnNames(), output.columnNames());

            final RealVector targetLengths = new ArrayRealVector(input.targets().stream()
                    .mapToDouble(t -> t.getEnd() - t.getStart() + 1).toArray());

            final RealMatrix inputCounts = input.counts();
            final RealMatrix outputCounts = output.counts();

            Assert.assertEquals(inputCounts.getRowDimension(), outputCounts.getRowDimension());
            Assert.assertEquals(inputCounts.getColumnDimension(), outputCounts.getColumnDimension());

            for (int col = 0; col < inputCounts.getColumnDimension(); col++) {
                final RealVector inputColumn = inputCounts.getColumnVector(col);
                final RealVector outputColumn = outputCounts.getColumnVector(col);
                final double inputWeightedAverage = inputColumn.dotProduct(targetLengths) / targetLengths.getL1Norm();
                Assert.assertEquals(inputColumn.cosine(outputColumn), 1, 1e-8);
                Assert.assertEquals(outputColumn.mapMultiply(inputWeightedAverage).getL1Distance(inputColumn), 0, 1e-10);
            }
        } catch (final IOException e) {
            //throw new Exception(e);
        }
    }
}