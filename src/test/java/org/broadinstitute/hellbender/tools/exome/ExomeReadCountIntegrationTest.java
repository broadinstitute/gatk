package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Integration tests for {@link ExomeReadCounts}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ExomeReadCountIntegrationTest extends CommandLineProgramTest {

    private static File testFile(final String fileName) {
        return new File(ExomeToolsTestUtils.getTestDataDir(),fileName);
    }


    //////////////////////
    // Test input files //
    //////////////////////

    private final static File INTERVALS_LIST = testFile("exome-read-counts-intervals.list");

    private final static File NA12878_BAM = testFile("exome-read-counts-NA12878.bam");

    private final static File NA12778_BAM = testFile("exome-read-counts-NA12778.bam");

    private final static File NA12872_BAM = testFile("exome-read-counts-NA12872.bam");

    private final static File[] ALL_BAMS = new File[] {NA12878_BAM,NA12778_BAM,NA12872_BAM};

    ////////////////////////
    //  Test output files //
    ////////////////////////

    private final static File NA12878_COHORT_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-NA12878.output");

    private final static File NA12878_SAMPLE_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-sample-NA12878.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-cohort.output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT =
            testFile("exome-read-counts-cohort.row-output");

    private final static File COHORT_COUNT_EXPECTED_COLUMN_OUTPUT =
            testFile("exome-read-counts-cohort.column-output");

    private final static File SAMPLE_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-sample.output");

    private final static File SAMPLE_COUNT_EXPECTED_ROW_OUTPUT =
            testFile("exome-read-counts-sample.row-output");

    private final static File SAMPLE_COUNT_EXPECTED_COLUMN_OUTPUT =
            testFile("exome-read-counts-sample.column-output");

    private final static File READ_GROUP_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-read-group.output");

    private final static File READ_GROUP_COUNT_EXPECTED_ROW_OUTPUT =
            testFile("exome-read-counts-read-group.row-output");

    private final static File READ_GROUP_COUNT_EXPECTED_COLUMN_OUTPUT =
            testFile("exome-read-counts-read-group.column-output");

    @Override
    public String getTestedClassName() {
        return ExomeReadCounts.class.getSimpleName();
    }

    /**
     * Creates a temp file making sure is going to be remove after testing VM finishes.
     * @param id temporal file name identifiable part.
     * @return never {@code null}.
     */
    private File createTempFile(final String id) {
        try {
            final File result = File.createTempFile("GATK4-ERCTest-" + id, ".tmp", null);
            result.deleteOnExit();
            return result;
        } catch (IOException ex) {
            Assert.fail("problems creating temporal file for output.",ex);
            throw new IllegalStateException("this code should never be reached");
        }
    }

    @DataProvider(name="correctRunData")
    public Object[][] correctRunData() {
        return new Object[][] {
                {       ALL_BAMS,
                        COHORT_COUNT_EXPECTED_OUTPUT,
                        COHORT_COUNT_EXPECTED_ROW_OUTPUT,
                        COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                        new String[0]
                },
                {       ALL_BAMS,
                        SAMPLE_COUNT_EXPECTED_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_ROW_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_COLUMN_OUTPUT,
                        new String[] { "-" + ExomeReadCounts.GROUP_BY_SHORT_NAME,ExomeReadCounts.GroupBy.SAMPLE.name()}
                },
                {       ALL_BAMS,
                        READ_GROUP_COUNT_EXPECTED_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_ROW_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_COLUMN_OUTPUT,
                        new String[] { "-" + ExomeReadCounts.GROUP_BY_SHORT_NAME,ExomeReadCounts.GroupBy.READ_GROUP.name()}
                }
        };
    }

    @Test(dataProvider = "correctRunData")
    public void testCorrectRun(final File[] bamFiles, final File expectedOutputFile, final File expectedRowOutputFile,
                               final File expectedColumnOutputFile,
                                final String[] additionalArguments) {
        final File outputFile = createTempFile("cohort-output");
        final File rowOutputFile = createTempFile("cohort-row-output");
        final File columnOutputFile = createTempFile("cohort-column-output");

        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-L", INTERVALS_LIST.getAbsolutePath(),
                "-" + ExomeReadCounts.ROW_SUMMARY_OUTPUT_SHORT_NAME, rowOutputFile.getAbsolutePath(),
                "-" + ExomeReadCounts.COLUMN_SUMMARY_OUTPUT_SHORT_NAME, columnOutputFile.getAbsolutePath()
        ));
        Arrays.asList(bamFiles).forEach(bam -> {
            arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
            arguments.add(bam.getAbsolutePath());
        });
        arguments.addAll(Arrays.asList(additionalArguments));
        System.err.println("COMMAND LINE: " + Arrays.toString(arguments.toArray()));
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists(), "output file does not exist: " + outputFile);
        Assert.assertTrue(outputFile.isFile(), "output file is not a regular file: " + outputFile);
        Assert.assertTrue(outputFile.canRead(), "output file cannot be read: " + outputFile);
        compareTableFiles(outputFile, expectedOutputFile, " matrix output ");
        if (expectedRowOutputFile != null) compareTableFiles(rowOutputFile, expectedRowOutputFile, " row output ");
        if (expectedColumnOutputFile != null) compareTableFiles(columnOutputFile, expectedColumnOutputFile, " column output ");

        checkTableConsistency(outputFile, rowOutputFile, columnOutputFile);
    }

    private void checkTableConsistency(final File outputFile, final File rowOutputFile, final File columnOutputFile) {
        final Table matrix;
        final Table row;
        final Table column;

        try {
            matrix = Table.fromFile(outputFile);
            row = rowOutputFile == null ? null : Table.fromFile(rowOutputFile);
            column = columnOutputFile == null ? null : Table.fromFile(columnOutputFile);

        } catch (final IOException ex) {
            Assert.fail("problem opening output files",ex);
            throw new IllegalStateException("unreachable test code");
        }

        if (row != null) {
            Assert.assertEquals(matrix.rowCount, row.rowCount,
                    "main output row count must match row output rows");

            for (int r = 0; r < matrix.rowCount; r++) {
                // Check coords agree:
                Assert.assertEquals(matrix.getString(r, 0), row.getString(r, 0));
                Assert.assertEquals(matrix.getLong(r, 1), row.getLong(r, 1));
                Assert.assertEquals(matrix.getLong(r, 2), row.getLong(r, 2));
                // Check row sum agree:
                long rowCount = 0;
                for (int c = 3; c < matrix.columnCount; c++)
                    rowCount += matrix.getLong(r, c);
                Assert.assertEquals(rowCount, row.getLong(r, 3));
                long size = matrix.getLong(r, 2) - matrix.getLong(r, 1) + 1;
                // Check consistent average per bp and column epsilon of 0.005 assume two decimal precision (~ %.2f).
                Assert.assertEquals(rowCount / ((double) size * (matrix.columnCount - 3)), row.getDouble(r, 4), 0.01);

            }
        }

        if (column != null) {
            Assert.assertEquals(matrix.columnCount - 3, column.rowCount,
                    "main output count column count must match column output rows");

            for (int c = 3; c < matrix.columnCount; c++) {
                Assert.assertEquals(matrix.columnNames[c], column.getString(c - 3, 0));
                long totalSum = 0;
                long totalBp = 0;
                for (int r = 0; r < matrix.rowCount; r++) {
                    totalBp += matrix.getLong(r, 2) - matrix.getLong(r, 1) + 1;
                    totalSum += matrix.getLong(r, c);
                }
                Assert.assertEquals(column.getLong(c - 3, 1), totalSum);
                Assert.assertEquals(column.getDouble(c - 3, 2), totalSum / (double) totalBp, 0.005);

            }
        }

    }

    private void compareTableFiles(File outputFile, File expectedOutputFile, final String role) {
        try {
            final Table actualTable = Table.fromFile(outputFile);
            final Table expectedTable = Table.fromFile(expectedOutputFile);
            Table.assertEquals(actualTable, expectedTable);
        } catch (final IOException ex) {
            Assert.fail("Failed to read actual or expected " +  role + " table files ", ex);
        }
    }

    public static class Table {

        private String[] values;
        private final int columnCount;
        private final int rowCount;
        private String[] columnNames;

        public static Table fromFile(final File file) throws IOException {
            final FileReader fr = new FileReader(file);
            final Table result = new Table(fr);
            fr.close();
            return result;
        }

        public Table(final Reader reader) throws IOException {
            if (reader == null) {
                throw new IllegalArgumentException("the reader cannot be null");
            }
            final BufferedReader lineReader = new BufferedReader(reader);
            final List<String> headerLines = new ArrayList<>(10);
            String lastLine;
            while ((lastLine = lineReader.readLine()) != null) {
                if (lastLine.matches("^#.*$"))
                    headerLines.add(lastLine);
                else
                    break;
            }
            final String[] header = lastLine == null ? null : lastLine.split("\\t");
            if (header == null) {
                throw new IllegalArgumentException("the table does not have a header");
            }
            columnCount = header.length;

            final List<String> values = new ArrayList<>(100);
            int lineNumber = 1 + 1 + headerLines.size();
            while ((lastLine = lineReader.readLine()) != null) {
                final String[] lineValues = lastLine.split("\\t");
                if (lineValues.length != columnCount)
                    throw new IllegalArgumentException("line " + lineNumber + " does not has the expected number of columns " + columnCount);
                values.addAll(Arrays.asList(lineValues));
            }

            this.values = values.toArray(new String[values.size()]);
            rowCount = values.size() / columnCount;
            columnNames = header;
        }

        public static void assertEquals(final Table tb1, final Table tb2) {

            // For now we ignore the header as is not a vital.
            // Assert.assertEquals(tb1.headerLines,tb2.headerLines),"Different header lines");
            Assert.assertEquals(tb1.columnCount,tb2.columnCount,"Different number of columns");
            Assert.assertEquals(tb1.rowCount,tb2.rowCount,"Different number of rows");
            Assert.assertEquals(tb1.columnNames,tb2.columnNames,"Differences in header");
            assertValuesAreEqual(tb1.values, tb2.values, 0.01, "Difference in values");

        }

        private static void assertValuesAreEqual(final String[] v1, final String[] v2, double epsilon, String message) {
            Assert.assertEquals(v1.length,v2.length,message);
            for (int i = 0; i < v1.length; i++) {
                final boolean v1IsDouble = v1[i].matches("^\\d*(\\.\\d*)?$");
                final boolean v2IsDouble = v2[i].matches("^\\d*(\\.\\d*)?$");
                Assert.assertEquals(v1IsDouble,v2IsDouble,message);
                if (v1IsDouble) {
                    Assert.assertEquals(Double.valueOf(v1[i]), Double.valueOf(v2[i]), epsilon, message);
                } else {
                    Assert.assertEquals(v1[i],v2[i],message);
                }
            }
        }

        public String getString(int r, int c) {
            return values[r * columnCount + c];
        }

        public long getLong(int r, int c) {
            return Long.valueOf(values[r * columnCount + c]);
        }

        public double getDouble(int r, int c) {
            return Double.valueOf(values[r * columnCount + c]);
        }
    }
}
