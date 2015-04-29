package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.recalibration.BQSRGatherer;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class BQSRGathererTest extends BaseTest {

    private static File recal1 = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.sg1.table.gz");
    private static File recal2 = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.sg2.table.gz");
    private static File recal3 = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.sg3.table.gz");
    private static File recal4 = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.sg4.table.gz");
    private static File recal5 = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.sg5.table.gz");
    private static File recalEmpty = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.empty.table.gz");

    private static File recal_original = new File(CommandLineProgramTest.getTestDataDir(), "HiSeq.1mb.1RG.noSG.table.gz");

    @Test
    public void testManyObservations() throws IOException {
        File recal = new File(new File(CommandLineProgramTest.getTestDataDir(), "BQSR"), "bqsr.manyObservations.piece.table.gz");

        final File output = BaseTest.createTempFile("BQSRgathererTest", ".table.gz");

        List<File> recalFiles = new LinkedList<>();
        for ( int i=0; i < 5; i++ )
            recalFiles.add(recal);

        BQSRGatherer gatherer = new BQSRGatherer();
        gatherer.gather(recalFiles, output);

        GATKReport originalReport = new GATKReport(new File(new File(CommandLineProgramTest.getTestDataDir(), "BQSR"), "bqsr.manyObservations.full.table.gz"));
        GATKReport calculatedReport = new GATKReport(output);

        testReports(originalReport, calculatedReport);
    }

    @Test
    public void testGatherBQSR() throws IOException {
        BQSRGatherer gatherer = new BQSRGatherer();
        List<File> recalFiles = new LinkedList<>();
        final File output = BaseTest.createTempFile("BQSRgathererTest", ".table.gz");

        recalFiles.add(recal1);
        recalFiles.add(recal2);
        recalFiles.add(recal3);
        recalFiles.add(recal4);
        recalFiles.add(recal5);
        gatherer.gather(recalFiles, output);

        GATKReport originalReport = new GATKReport(recal_original);
        GATKReport calculatedReport = new GATKReport(output);

        testReports(originalReport, calculatedReport);
    }

    @Test
    public void testGatherBQSRWithEmptyFile() throws IOException {
        BQSRGatherer gatherer = new BQSRGatherer();
        List<File> recalFiles = new LinkedList<>();
        final File output = BaseTest.createTempFile("BQSRgathererTest", ".table.gz");

        recalFiles.add(recal1);
        recalFiles.add(recal2);
        recalFiles.add(recal3);
        recalFiles.add(recal4);
        recalFiles.add(recal5);
        recalFiles.add(recalEmpty);
        gatherer.gather(recalFiles, output);

        GATKReport originalReport = new GATKReport(recal_original);
        GATKReport calculatedReport = new GATKReport(output);

        testReports(originalReport, calculatedReport);
    }

    private void testReports(final GATKReport originalReport, final GATKReport calculatedReport) {

        // test the Arguments table
        List<String> columnsToTest = Arrays.asList(RecalUtils.ARGUMENT_COLUMN_NAME, RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
        GATKReportTable originalTable = originalReport.getTable(RecalUtils.ARGUMENT_REPORT_TABLE_TITLE);
        GATKReportTable calculatedTable = calculatedReport.getTable(RecalUtils.ARGUMENT_REPORT_TABLE_TITLE);
        testTablesWithColumns(originalTable, calculatedTable, columnsToTest);

        // test the Quantized table
        columnsToTest = Arrays.asList(RecalUtils.QUALITY_SCORE_COLUMN_NAME, RecalUtils.QUANTIZED_COUNT_COLUMN_NAME, RecalUtils.QUANTIZED_VALUE_COLUMN_NAME);
        originalTable = originalReport.getTable(RecalUtils.QUANTIZED_REPORT_TABLE_TITLE);
        calculatedTable = calculatedReport.getTable(RecalUtils.QUANTIZED_REPORT_TABLE_TITLE);
        testTablesWithColumns(originalTable, calculatedTable, columnsToTest);

        // test the RecalTable0 table
        columnsToTest = Arrays.asList(RecalUtils.READGROUP_COLUMN_NAME, RecalUtils.EVENT_TYPE_COLUMN_NAME, RecalUtils.ESTIMATED_Q_REPORTED_COLUMN_NAME, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, RecalUtils.NUMBER_ERRORS_COLUMN_NAME);
        originalTable = originalReport.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE);
        calculatedTable = calculatedReport.getTable(RecalUtils.READGROUP_REPORT_TABLE_TITLE);
        testTablesWithColumns(originalTable, calculatedTable, columnsToTest);

        // test the RecalTable1 table
        columnsToTest = Arrays.asList(RecalUtils.READGROUP_COLUMN_NAME, RecalUtils.QUALITY_SCORE_COLUMN_NAME, RecalUtils.EVENT_TYPE_COLUMN_NAME, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, RecalUtils.NUMBER_ERRORS_COLUMN_NAME);
        originalTable = originalReport.getTable(RecalUtils.QUALITY_SCORE_REPORT_TABLE_TITLE);
        calculatedTable = calculatedReport.getTable(RecalUtils.QUALITY_SCORE_REPORT_TABLE_TITLE);
        testTablesWithColumns(originalTable, calculatedTable, columnsToTest);

        // test the RecalTable2 table
        columnsToTest = Arrays.asList(RecalUtils.READGROUP_COLUMN_NAME, RecalUtils.QUALITY_SCORE_COLUMN_NAME, RecalUtils.COVARIATE_VALUE_COLUMN_NAME, RecalUtils.COVARIATE_NAME_COLUMN_NAME, RecalUtils.EVENT_TYPE_COLUMN_NAME, RecalUtils.NUMBER_OBSERVATIONS_COLUMN_NAME, RecalUtils.NUMBER_ERRORS_COLUMN_NAME);
        originalTable = originalReport.getTable(RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE);
        calculatedTable = calculatedReport.getTable(RecalUtils.ALL_COVARIATES_REPORT_TABLE_TITLE);
        testTablesWithColumns(originalTable, calculatedTable, columnsToTest);
    }

    /**
     * Common testing functionality given the columns to test and the multiplication factor to the expected result
     *
     * @param original the original table
     * @param calculated the calculated table
     * @param columnsToTest list of columns to test. All columns will be tested with the same criteria (equality given factor)
     */
    private void testTablesWithColumns(GATKReportTable original, GATKReportTable calculated, List<String> columnsToTest) {
        for (int row = 0; row < original.getNumRows(); row++ ) {
            for (String column : columnsToTest) {
                Object actual = calculated.get(new Integer(row), column);
                Object expected = original.get(row, column);
                //if ( !actual.equals(expected) )
                //    System.out.println("Row=" + row + " Table=" + original.getTableName() + " Column=" + column + " Expected=" + expected + " Actual=" + actual);
                Assert.assertEquals(actual, expected, "Row: " + row + " Original Table: " + original.getTableName() + " Column=" + column);
            }
        }
        
    }

    @Test
    public void testGatherMissingReadGroup() {
        // TODO: This test data is NOT private but privateTestDir offers:
        // TODO:   - Doesn't end up in protected / private github
        // TODO:   - IS otherwise available for local debugging unlike /humgen NFS mounts
        // Hand modified subset of problematic gather inputs submitted by George Grant
        File input1 = new File(CommandLineProgramTest.getTestDataDir(), "NA12878.rg_subset.chr1.recal_data.table.gz");
        File input2 = new File(CommandLineProgramTest.getTestDataDir(), "NA12878.rg_subset.chrY_Plus.recal_data.table.gz");

        GATKReport report12 = BQSRGatherer.gatherReport(Arrays.asList(input1, input2));
        GATKReport report21 = BQSRGatherer.gatherReport(Arrays.asList(input2, input1));

        Assert.assertTrue(report12.equals(report21), "GATK reports are different when gathered in a different order.");
    }
}
