package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Integration tests for {@link CalculateTargetCoverage}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CalculateTargetCoverageIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/calculatetargetcoverage");
    private static File testFile(final String fileName) {
        return new File(TEST_DIR, fileName);
    }

    //////////////////////
    // Test input files //
    //////////////////////

    private final static File INTERVALS_LIST = testFile("exome-read-counts-intervals.list");

    private final static File INTERVALS_BED = testFile("exome-read-counts-intervals.tsv");

    private final static File VCF_FEATURE_FILE = testFile("random-variant-file.vcf");

    private final static File INTERVALS_BED_MISSING_NAMES = testFile("exome-read-counts-intervals-missing-names.tsv");

    private final static File NA12878_BAM = testFile("exome-read-counts-NA12878.bam");

    private final static File NA12778_BAM = testFile("exome-read-counts-NA12778.bam");

    private final static File NA12872_BAM = testFile("exome-read-counts-NA12872.bam");

    private final static File[] ALL_BAMS = new File[] {NA12878_BAM,NA12778_BAM,NA12872_BAM};

    private final static File[] DUPREADS_BAM = new File[] {testFile("dupReadsMini.bam")};

    private final static File INTERVALS_LIST_DUPS = testFile("exome-read-counts-intervals_dups.list");

    private final static File INTERVALS_BED_DUPS = testFile("exome-read-counts-intervals_dups.tsv");

    ////////////////////////
    //  Test output files //
    ////////////////////////

    private final static File COHORT_COUNT_EXPECTED_OUTPUT_WITH_NAMES =
            testFile("exome-read-counts-cohort-with-names.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT_WITH_NAMES_ONLY =
            testFile("exome-read-counts-cohort-with-names-only.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES =
            testFile("exome-read-counts-cohort-with-BED-names.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES_DUPS =
            testFile("exome-read-counts-cohort-with-BED-names_dups.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_MISSING_NAMES =
            testFile("exome-read-counts-cohort-with-BED-missing-names.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES_ONLY =
            testFile("exome-read-counts-cohort-with-BED-names-only.output");

    private final static File COHORT_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-cohort.output");

    private final static File COHORT_COUNT_EXPECTED_PCOV_OUTPUT =
            testFile("exome-read-counts-cohort.pcov-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT =
            testFile("exome-read-counts-cohort.row-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_NAMES =
            testFile("exome-read-counts-cohort-with-names.row-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_NAMES_ONLY =
            testFile("exome-read-counts-cohort-with-names-only.row-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES =
            testFile("exome-read-counts-cohort-with-BED-names.row-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES_DUPS =
            testFile("exome-read-counts-cohort-with-BED-names_dups.row-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_MISSING_NAMES =
            testFile("exome-read-counts-cohort-with-BED-missing-names.row-output");

    private final static File COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES_ONLY =
            testFile("exome-read-counts-cohort-with-BED-names-only.row-output");

    private final static File COHORT_COUNT_EXPECTED_COLUMN_OUTPUT =
            testFile("exome-read-counts-cohort.column-output");

    private final static File COHORT_COUNT_EXPECTED_COLUMN_OUTPUT_DUPS =
            testFile("exome-read-counts-cohort_dups.column-output");

    private final static File SAMPLE_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-sample.output");

    private final static File SAMPLE_COUNT_EXPECTED_PCOV_OUTPUT =
            testFile("exome-read-counts-sample.pcov-output");

    private final static File SAMPLE_COUNT_EXPECTED_ROW_OUTPUT =
            testFile("exome-read-counts-sample.row-output");

    private final static File SAMPLE_COUNT_EXPECTED_COLUMN_OUTPUT =
            testFile("exome-read-counts-sample.column-output");

    private final static File READ_GROUP_COUNT_EXPECTED_OUTPUT =
            testFile("exome-read-counts-read-group.output");

    private final static File READ_GROUP_COUNT_EXPECTED_PCOV_OUTPUT =
            testFile("exome-read-counts-read-group.pcov-output");

    private final static File READ_GROUP_COUNT_EXPECTED_ROW_OUTPUT =
            testFile("exome-read-counts-read-group.row-output");

    private final static File READ_GROUP_COUNT_EXPECTED_COLUMN_OUTPUT =
            testFile("exome-read-counts-read-group.column-output");

    @Override
    public String getTestedClassName() {
        return CalculateTargetCoverage.class.getSimpleName();
    }

    /**
     * Creates a temp file making sure is going to be remove after testing VM finishes.
     * @param id temporal file name identifiable part.
     * @return never {@code null}.
     */
    private File createTempFile(final String id) {
        return BaseTest.createTempFile("GATK4-ERCTest-" + id,".tmp");
    }

    @DataProvider(name="correctRunData")
    public Object[][] correctRunData() {
        return new Object[][] {
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        COHORT_COUNT_EXPECTED_OUTPUT,
                        COHORT_COUNT_EXPECTED_ROW_OUTPUT,
                        COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.RAW,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[0]
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        SAMPLE_COUNT_EXPECTED_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_ROW_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.RAW,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.GROUP_BY_SHORT_NAME, CalculateTargetCoverage.GroupBy.SAMPLE.name()}
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        READ_GROUP_COUNT_EXPECTED_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_ROW_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.RAW,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.GROUP_BY_SHORT_NAME, CalculateTargetCoverage.GroupBy.READ_GROUP.name()}
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        COHORT_COUNT_EXPECTED_PCOV_OUTPUT,
                        COHORT_COUNT_EXPECTED_ROW_OUTPUT,
                        COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.PCOV,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[0]
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        SAMPLE_COUNT_EXPECTED_PCOV_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_ROW_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.PCOV,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.GROUP_BY_SHORT_NAME, CalculateTargetCoverage.GroupBy.SAMPLE.name()}
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        READ_GROUP_COUNT_EXPECTED_PCOV_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_ROW_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.PCOV,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.GROUP_BY_SHORT_NAME, CalculateTargetCoverage.GroupBy.READ_GROUP.name()}
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        COHORT_COUNT_EXPECTED_OUTPUT,
                        COHORT_COUNT_EXPECTED_ROW_OUTPUT,
                        COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.RAW,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME,INTERVALS_BED.getAbsolutePath()},
                },

                {       ALL_BAMS,
                        INTERVALS_LIST,
                        SAMPLE_COUNT_EXPECTED_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_ROW_OUTPUT,
                        SAMPLE_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.RAW,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.GROUP_BY_SHORT_NAME, CalculateTargetCoverage.GroupBy.SAMPLE.name(),
                                "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME,INTERVALS_BED.getAbsolutePath()}
                },
                {       ALL_BAMS,
                        INTERVALS_LIST,
                        READ_GROUP_COUNT_EXPECTED_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_ROW_OUTPUT,
                        READ_GROUP_COUNT_EXPECTED_COLUMN_OUTPUT,
                        CalculateTargetCoverage.Transform.RAW,
                        CalculateTargetCoverage.TargetOutInfo.COORDS,
                        new String[] { "-" + CalculateTargetCoverage.GROUP_BY_SHORT_NAME, CalculateTargetCoverage.GroupBy.READ_GROUP.name(),
                                "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME,INTERVALS_BED.getAbsolutePath()}
                },{     ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.FULL,
                new String[0]
        },{     ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_NAMES_ONLY,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_NAMES_ONLY,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.NAME,
                new String[0]
        },{     ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.FULL,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_BED.getAbsolutePath() }
        },{     ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES_ONLY,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES_ONLY,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.NAME,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_BED.getAbsolutePath() }
        },{     ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.FULL,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_BED_MISSING_NAMES.getAbsolutePath() }
        },{     DUPREADS_BAM,
                INTERVALS_LIST_DUPS,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES_DUPS,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES_DUPS,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT_DUPS,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.FULL,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_BED_DUPS.getAbsolutePath(),
                        "-" + StandardArgumentDefinitions.DISABLE_READ_FILTER_SHORT_NAME, ReadFilterLibrary.NOT_DUPLICATE.getClass().getSimpleName() },
        }
        };
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingTargetNameRun() {
        testCorrectRun(ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.NAME,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_BED_MISSING_NAMES.getAbsolutePath() }
        );
    }

    @Test(expectedExceptions = UserException.class)
    public void testNoIntervalFile() {
        testCorrectRun(ALL_BAMS,
                null,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.NAME,
                new String[0]
        );
    }

    @Test
    public void testTargetFileOnly() {
        testCorrectRun(ALL_BAMS,
                null,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.FULL,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_BED.getAbsolutePath() }
        );
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testVcfTargetFile() {
        testCorrectRun(ALL_BAMS,
                null,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.FULL,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, VCF_FEATURE_FILE.getAbsolutePath() }
        );
    }

    @Test(expectedExceptions = UserException.class)
    public void testBadTargetFile() {
        testCorrectRun(ALL_BAMS,
                INTERVALS_LIST,
                COHORT_COUNT_EXPECTED_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_ROW_OUTPUT_WITH_BED_MISSING_NAMES,
                COHORT_COUNT_EXPECTED_COLUMN_OUTPUT,
                CalculateTargetCoverage.Transform.RAW,
                CalculateTargetCoverage.TargetOutInfo.NAME,
                new String[] { "-" + CalculateTargetCoverage.TARGET_FILE_SHORT_NAME, INTERVALS_LIST.getAbsolutePath() }
        );
    }


    @Test(dataProvider = "correctRunData")
    public void testCorrectRun(final File[] bamFiles, final File intervalFile, final File expectedOutputFile, final File expectedRowOutputFile,
                               final File expectedColumnOutputFile, final CalculateTargetCoverage.Transform transform,
                               final CalculateTargetCoverage.TargetOutInfo targetOutInfo,
                               final String[] additionalArguments) {
        final File outputFile = createTempFile("cohort-output");
        final File rowOutputFile = createTempFile("cohort-row-output");
        final File columnOutputFile = createTempFile("cohort-column-output");

        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + CalculateTargetCoverage.ROW_SUMMARY_OUTPUT_SHORT_NAME, rowOutputFile.getAbsolutePath(),
                "-" + CalculateTargetCoverage.COLUMN_SUMMARY_OUTPUT_SHORT_NAME, columnOutputFile.getAbsolutePath(),
                "-" + CalculateTargetCoverage.TRANSFORM_SHORT_NAME, transform.name(),
                "-" + CalculateTargetCoverage.TARGET_OUT_INFO_SHORT_NAME, targetOutInfo.name()
        ));
        if (intervalFile != null) {
            arguments.add("-L");
            arguments.add(intervalFile.getAbsolutePath());
        }
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
        compareTableFiles(outputFile, expectedOutputFile, "matrix output: " + outputFile + " " + expectedOutputFile);
        if (expectedRowOutputFile != null) compareTableFiles(rowOutputFile, expectedRowOutputFile, " row output ");
        if (expectedColumnOutputFile != null) compareTableFiles(columnOutputFile, expectedColumnOutputFile, " column output ");

        switch (transform) {
            case RAW:
                checkTableRawConsistency(targetOutInfo, outputFile, rowOutputFile, columnOutputFile);
                break;
            case PCOV:
                checkTablePcovConsistency(targetOutInfo, outputFile, rowOutputFile, columnOutputFile);
                break;
            default:
                Assert.fail("bug in testing code: unexpected transform " + transform.name());
        }
    }

    private void checkTableRawConsistency(final CalculateTargetCoverage.TargetOutInfo targetOutInfo, final File outputFile, final File rowOutputFile, final File columnOutputFile) {
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
                // Check coordinates agree:
                Assert.assertEquals(matrix.getString(r, 0), row.getString(r, 0));
                if (targetOutInfo != CalculateTargetCoverage.TargetOutInfo.NAME) {
                    Assert.assertEquals(matrix.getLong(r, 1), row.getLong(r, 1));
                    Assert.assertEquals(matrix.getLong(r, 2), row.getLong(r, 2));
                }
                // Check row sum agree:
                long rowCount = 0;
                for (int c = targetOutInfo.columnCount(); c < matrix.columnCount; c++)
                    rowCount += matrix.getLong(r, c);
                Assert.assertEquals(rowCount, row.getLong(r, targetOutInfo.columnCount()));
                // Check consistent average per bp and column epsilon of 0.005 assume two decimal precision (~ %.2f).
                if (targetOutInfo != CalculateTargetCoverage.TargetOutInfo.NAME) {
                    final long size = matrix.getLong(r, 2) - matrix.getLong(r, 1) + 1;
                    Assert.assertEquals(rowCount / ((double) size * (matrix.columnCount - targetOutInfo.columnCount())), row.getDouble(r, targetOutInfo.columnCount() + 1), 0.01);
                }
            }
        }

        if (column != null) {
            Assert.assertEquals(matrix.columnCount - targetOutInfo.columnCount(), column.rowCount,
                    "main output count column count must match column output rows");

            for (int c = targetOutInfo.columnCount(); c < matrix.columnCount; c++) {
                Assert.assertEquals(matrix.columnNames[c], column.getString(c - targetOutInfo.columnCount(), 0));
                long totalSum = 0;
                for (int r = 0; r < matrix.rowCount; r++) {
                    totalSum += matrix.getLong(r, c);
                }
                Assert.assertEquals(column.getLong(c - targetOutInfo.columnCount(), 1), totalSum);
                if (targetOutInfo != CalculateTargetCoverage.TargetOutInfo.NAME) {
                    long totalBp = 0;
                    for (int r = 0; r < matrix.rowCount; r++) {
                        totalBp += matrix.getLong(r, 2) - matrix.getLong(r, 1) + 1;
                    }
                    Assert.assertEquals(column.getDouble(c - targetOutInfo.columnCount(), 2), totalSum / (double) totalBp, 0.005);
                }
            }
        }

    }

    private void checkTablePcovConsistency(final CalculateTargetCoverage.TargetOutInfo targetOutInfo, final File outputFile, final File rowOutputFile, final File columnOutputFile) {
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

        // without the column table we cannot do much in terms of comparisons
        // as the main matrix output is a proportion of those column totals.
        if (column == null) {
            return;
        }

        for (int i = targetOutInfo.columnCount(); i < matrix.columnCount; i++) {
            if (column.getLong(i- targetOutInfo.columnCount(),1) > 0) {
                double sumPcov = 0;
                for (int j = 0; j < matrix.rowCount; j++) {
                    sumPcov += matrix.getDouble(j,i);
                }
                Assert.assertEquals(sumPcov,1.0,0.0005);
            } else {
                for (int j = 0; j < matrix.rowCount; j++) {
                    Assert.assertTrue(Double.isNaN(matrix.getDouble(j, i)));
                }
            }
        }

        Assert.assertEquals(matrix.columnCount - targetOutInfo.columnCount(), column.rowCount,
                "main output count column count must match column output rows");

        for (int c = targetOutInfo.columnCount(); c < matrix.columnCount; c++) {
            Assert.assertEquals(matrix.columnNames[c], column.getString(c - targetOutInfo.columnCount(), 0));
            long totalBp = 0;
            for (int r = 0; r < matrix.rowCount; r++) {
                totalBp += matrix.getLong(r, 2) - matrix.getLong(r, 1) + 1;
            }
            Assert.assertEquals(column.getDouble(c - targetOutInfo.columnCount(), 2), column.getLong(c - targetOutInfo.columnCount(), 1) / (double) totalBp, 0.005);
        }

        if (row != null) {
            Assert.assertEquals(matrix.rowCount, row.rowCount,
                    "main output row count must match row output rows");

            for (int r = 0; r < matrix.rowCount; r++) {
                // Check coordinates agree:
                Assert.assertEquals(matrix.getString(r, 0), row.getString(r, 0));
                if (targetOutInfo != CalculateTargetCoverage.TargetOutInfo.NAME) {
                    Assert.assertEquals(matrix.getLong(r, 1), row.getLong(r, 1));
                    Assert.assertEquals(matrix.getLong(r, 2), row.getLong(r, 2));
                }
                // Check row sum agree:
                long rowCount = 0;
                for (int c = targetOutInfo.columnCount(); c < matrix.columnCount; c++) {
                    rowCount += Math.round(matrix.getDouble(r, c) * column.getLong(c - targetOutInfo.columnCount(), 1));
                }
                Assert.assertEquals(rowCount, row.getLong(r, targetOutInfo.columnCount()));
                // Check consistent average per bp and column epsilon of 0.005 assume two decimal precision (~ %.2f).
                if (targetOutInfo != CalculateTargetCoverage.TargetOutInfo.NAME) {
                    final long size = matrix.getLong(r, 2) - matrix.getLong(r, 1) + 1;
                    Assert.assertEquals(rowCount / ((double) size * (matrix.columnCount - targetOutInfo.columnCount())), row.getDouble(r, targetOutInfo.columnCount() + 1), 0.01);
                }
            }
        }
    }

    private void compareTableFiles(File outputFile, File expectedOutputFile, final String role) {
        try {
            final Table actualTable = Table.fromFile(outputFile);
            final Table expectedTable = Table.fromFile(expectedOutputFile);
            Table.assertEquals(actualTable, expectedTable, outputFile, expectedOutputFile);
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
            try (final FileReader fr = new FileReader(file)) {
                return new Table(fr);
            }
        }

        public Table(final Reader reader) throws IOException {
            Utils.nonNull(reader, "the reader cannot be null");
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
            Utils.nonNull(header, "the table does not have a header");
            columnCount = header.length;

            final List<String> values = new ArrayList<>(100);
            int lineNumber = 1 + 1 + headerLines.size();
            while ((lastLine = lineReader.readLine()) != null) {
                final String[] lineValues = lastLine.split("\\t");
                if (lineValues.length != columnCount) {
                    throw new IllegalArgumentException("line " + lineNumber + " does not has the expected number of columns " + columnCount + ": " + lastLine);
                }
                values.addAll(Arrays.asList(lineValues));
            }

            this.values = values.toArray(new String[values.size()]);
            rowCount = values.size() / columnCount;
            columnNames = header;
        }

        public static void assertEquals(final Table tb1, final Table tb2, final File f1, final File f2) {

            // For now we ignore the header as is not a vital.
            // Assert.assertEquals(tb1.headerLines,tb2.headerLines),"Different header lines");
            Assert.assertEquals(tb1.columnCount,tb2.columnCount,"Different number of columns");
            Assert.assertEquals(tb1.rowCount,tb2.rowCount,"Different number of rows");
            Assert.assertEquals(tb1.columnNames,tb2.columnNames,"Differences in header");
            assertValuesAreEqual(tb1.values, tb2.values, 0.01, "Difference in values between " + f1 + " and " + f2);

        }

        private static void assertValuesAreEqual(final String[] v1, final String[] v2, double epsilon, String message) {
            Assert.assertEquals(v1.length, v2.length, message);
            for (int i = 0; i < v1.length; i++) {
                final boolean v1IsDouble = isADouble(v1[i]);
                final boolean v2IsDouble = isADouble(v2[i]);
                Assert.assertEquals(v1IsDouble,v2IsDouble,message);
                if (v1IsDouble) {
                    if (Double.isNaN(Double.valueOf(v1[i]))) {
                        Assert.assertTrue(Double.isNaN(Double.valueOf(v2[i])));
                    } else {
                        Assert.assertEquals(Double.valueOf(v1[i]), Double.valueOf(v2[i]), epsilon, message);
                    }
                } else {
                    Assert.assertEquals(v1[i],v2[i],message);
                }
            }
        }

        private static boolean isADouble(final String s) {
            if (s == null) {
                return false;
            } else {
                try {
                    Double.parseDouble(s);
                } catch (final NumberFormatException n) {
                    return false;
                }
                return true;
            }
        }

        public String getString(int r, int c) {
            return values[r * columnCount + c];
        }

        public long getLong(int r, int c) {
            return Long.parseLong(values[r * columnCount + c]);
        }

        public double getDouble(int r, int c) {
            return Double.parseDouble(values[r * columnCount + c]);
        }
    }
}
