package org.broadinstitute.hellbender.tools.exome;


import au.com.bytecode.opencsv.CSVWriter;
import com.opencsv.CSVReader;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.hdf5.PoNTestUtils;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;
import org.broadinstitute.hellbender.utils.svd.SVDTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DecomposeSingularValuesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File CONTROL_PCOV_FULL_FILE = new File(TEST_FILE_DIR, "create-pon-control-full.pcov");

    @Test
    public void testBasicSparkSVD() {
        final File outputFileV = createTempFile("svd-",".tsv");
        final File outputFileS = createTempFile("svd-",".tsv");
        final File outputFileU = createTempFile("svd-",".tsv");

        final List<String> arguments = createDefaultArguments(outputFileV, outputFileS, outputFileU);

        runCommandLine(arguments);

        assertOutputFile(outputFileU, 10000, 100);
        assertOutputFile(outputFileV, 100, 100);
        assertOutputFile(outputFileS, 100, 100);
        assertSVDValues(outputFileV, outputFileS, outputFileU);

        assertSVDValues(outputFileV, outputFileS, outputFileU);
    }

    private List<String> createDefaultArguments(final File outputFileV, final File outputFileS, final File outputFileU) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + DecomposeSingularValues.INPUT_FILE_SHORT_NAME);
        arguments.add(CONTROL_PCOV_FULL_FILE.getAbsolutePath());

        arguments.add("-" + DecomposeSingularValues.OUTPUT_V_FILE_SHORT_NAME);
        arguments.add(outputFileV.toString());

        arguments.add("-" + DecomposeSingularValues.OUTPUT_S_FILE_SHORT_NAME);
        arguments.add(outputFileS.toString());

        arguments.add("-" + DecomposeSingularValues.OUTPUT_U_FILE_SHORT_NAME);
        arguments.add(outputFileU.toString());

        arguments.add("--verbosity");
        arguments.add("WARNING");
        return arguments;
    }

    @Test
    public void testBasicSparkDisabledSVD() {
        final File outputFileV = createTempFile("svd-",".tsv");
        final File outputFileS = createTempFile("svd-",".tsv");
        final File outputFileU = createTempFile("svd-",".tsv");

        final List<String> arguments = createDefaultArguments(outputFileV, outputFileS, outputFileU);

        arguments.add("-" + SparkToggleCommandLineProgram.DISABLE_SPARK_SHORT_NAME);

        runCommandLine(arguments);

        assertOutputFile(outputFileU, 10000, 100);
        assertOutputFile(outputFileV, 100, 100);
        assertOutputFile(outputFileS, 100, 100);

        assertSVDValues(outputFileV, outputFileS, outputFileU);
    }

    private void assertSVDValues(final File outputFileV, final File outputFileS, final File outputFileU) {
        try {
            final ReadCountCollection rcc = ReadCountCollectionUtils.parse(CONTROL_PCOV_FULL_FILE);
            final SVD svd = SVDFactory.createSVD(rcc.counts());
            final RealMatrix sDiag = new DiagonalMatrix(svd.getSingularValues());
            assertOutputFileValues(outputFileU, svd.getU());
            assertOutputFileValues(outputFileS, sDiag);
            assertOutputFileValues(outputFileV, svd.getV());
            SVDTestUtils.assertUnitaryMatrix(svd.getV());
            SVDTestUtils.assertUnitaryMatrix(svd.getU());
            Assert.assertTrue(MatrixUtils.isSymmetric(sDiag, 1e-32));

        } catch (final IOException ioe) {
            Assert.fail("Could not open test file: " + CONTROL_PCOV_FULL_FILE, ioe);
        }
    }

    private void assertOutputFileValues(final File outputFile, final RealMatrix expected){

        final RealMatrix actualResults = PoNTestUtils.readTsvIntoMatrix(outputFile);

        Assert.assertEquals(actualResults.getColumnDimension(), expected.getColumnDimension());
        Assert.assertEquals(actualResults.getRowDimension(), expected.getRowDimension());

        for (int i = 0; i < actualResults.getRowDimension(); i ++) {
            for (int j = 0; j < actualResults.getColumnDimension(); j ++) {
                Assert.assertEquals(Math.abs(actualResults.getEntry(i,j)), Math.abs(expected.getEntry(i,j)), 1e-7);
            }
        }
    }

    private void assertOutputFile(final File outputFile, final int gtNumRows, final int gtNumCols){
        Assert.assertTrue(outputFile.exists(), "Output file does not exist.");
        Assert.assertTrue(outputFile.length() > 0, "Output file is empty: " + outputFile.getAbsolutePath());
        final RealMatrix actualResults = PoNTestUtils.readTsvIntoMatrix(outputFile);
        Assert.assertEquals(actualResults.getRowDimension(), gtNumRows);
        Assert.assertEquals(actualResults.getColumnDimension(), gtNumCols);
    }
}
