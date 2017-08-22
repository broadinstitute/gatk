package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;


public final class GATKReportUnitTest extends GATKBaseTest {

    private static final String TABLE_NAME = "TableName";

    @Test
    public void testParse() throws Exception {
        String reportPath = publicTestDir + "exampleGATKReportv2.tbl";
        GATKReport report = new GATKReport(reportPath);
        Assert.assertEquals(report.getVersion(), GATKReportVersion.V1_1);
        Assert.assertEquals(report.getTables().size(), 5);

        GATKReportTable countVariants = report.getTable("CountVariants");
        Assert.assertEquals(countVariants.get(0, "nProcessedLoci"), "63025520");
        Assert.assertEquals(countVariants.get(0, "nNoCalls"), "0");
        Assert.assertEquals(countVariants.get(0, "heterozygosity"), 4.73e-06);

        GATKReportTable validationReport = report.getTable("ValidationReport");
        Assert.assertEquals(validationReport.get(2, "PPV"), Double.NaN);
    }

    @DataProvider(name = "rightAlignValues")
    public Object[][] getRightAlignValues() {
        return new Object[][]{
                new Object[]{null, true},
                new Object[]{"null", true},
                new Object[]{"NA", true},
                new Object[]{"0", true},
                new Object[]{"0.0", true},
                new Object[]{"-0", true},
                new Object[]{"-0.0", true},
                new Object[]{String.valueOf(Long.MAX_VALUE), true},
                new Object[]{String.valueOf(Long.MIN_VALUE), true},
                new Object[]{String.valueOf(Float.MIN_NORMAL), true},
                new Object[]{String.valueOf(Double.MAX_VALUE), true},
                new Object[]{String.valueOf(Double.MIN_VALUE), true},
                new Object[]{String.valueOf(Double.POSITIVE_INFINITY), true},
                new Object[]{String.valueOf(Double.NEGATIVE_INFINITY), true},
                new Object[]{String.valueOf(Double.NaN), true},
                new Object[]{"hello", false}
        };
    }

    @Test(dataProvider = "rightAlignValues")
    public void testIsRightAlign(String value, boolean expected) {
        Assert.assertEquals(GATKReportColumn.isRightAlign(value), expected, "right align of '" + value + "'");
    }

    private GATKReport getTableWithRandomValues() {
        Random number = new Random(123L);
        final int MAX_VALUE = 20;

        GATKReport report = GATKReport.newSimpleReportWithDescription(TABLE_NAME, "table with random values sorted by columns", GATKReportTable.Sorting.SORT_BY_COLUMN, "col1", "col2", "col3");

        final int NUM_ROWS = 100;
        for (int x = 0; x < NUM_ROWS; x++) {
            report.addRow(number.nextInt(MAX_VALUE), number.nextInt(MAX_VALUE), number.nextInt(MAX_VALUE));
        }
        return report;
    }

    @Test
    public void testSortingByColumn() throws IOException {
        assertTableSortsCorrectly(getTableWithRandomValues());
    }

    private void assertTableSortsCorrectly(GATKReport inputReport) throws IOException {
        File testingSortingTableFile = createTempFile("testSortingFile",".txt");

        try (final PrintStream ps = new PrintStream(testingSortingTableFile)) {
            inputReport.print(ps);
        }

        final GATKReport reloaded = new GATKReport(testingSortingTableFile);
        final GATKReportTable reloadedTable = reloaded.getTable(TABLE_NAME);
        final int numRows = reloadedTable.getNumRows();
        for( int row = 0; row < numRows-1; row++){
            Assert.assertTrue(compareRows(reloadedTable, row, row+1) <= 0);
        }
    }

    private static int compareRows(GATKReportTable table, int firstRow, int secondRow){
        for( int i = 0; i < table.getNumColumns(); i++){
            final Integer first = Integer.parseInt((String) table.get(firstRow, i));
            final Integer second = Integer.parseInt((String) table.get(secondRow, i));
            final int comparison = first.compareTo(second);
            if (comparison != 0)  {
                return comparison;
            }
        }
        return 0;
    }


    private GATKReportTable makeBasicTable() {
        GATKReport report = GATKReport.newSimpleReport(TABLE_NAME, GATKReportTable.Sorting.SORT_BY_COLUMN, "sample", "value");
        GATKReportTable table = report.getTable(TABLE_NAME);
        report.addRow("foo.1", "hello");
        report.addRow("foo.2", "world");
        return table;
    }

    @Test
    public void testDottedSampleName() {
        GATKReportTable table = makeBasicTable();
        Assert.assertEquals(table.get(0, "value"), "hello");
        Assert.assertEquals(table.get(1, "value"), "world");
    }

    @Test
    public void testSimpleGATKReport() throws FileNotFoundException {
        // Create a new simple GATK report named "TableName" with columns: Roger, is, and Awesome
        GATKReport report = GATKReport.newSimpleReport(TABLE_NAME, GATKReportTable.Sorting.SORT_BY_COLUMN, "Roger", "is", "Awesome");

        // Add data to simple GATK report
        report.addRow(12, 23.45, true);
        report.addRow("ans", '3', 24.5);
        report.addRow("hi", "", 2.3);

        File file = createTempFile("GATKReportGatherer-UnitTest", ".tbl");
        try (PrintStream ps = new PrintStream(file)) {
            report.print(ps);
            GATKReport inputRead = new GATKReport(file);
            Assert.assertTrue(report.isSameFormat(inputRead));
        }

    }

    @Test
    public void testGATKReportGatherer() throws FileNotFoundException {

        GATKReport report1, report2, report3;
        report1 = new GATKReport();
        report1.addTable(TABLE_NAME, "Description", 2);
        report1.getTable(TABLE_NAME).addColumn("colA", "%s");
        report1.getTable(TABLE_NAME).addColumn("colB", "%c");
        report1.getTable(TABLE_NAME).set(0, "colA", "NotNum");
        report1.getTable(TABLE_NAME).set(0, "colB", (char) 64);

        report2 = new GATKReport();
        report2.addTable(TABLE_NAME, "Description", 2);
        report2.getTable(TABLE_NAME).addColumn("colA", "%s");
        report2.getTable(TABLE_NAME).addColumn("colB", "%c");
        report2.getTable(TABLE_NAME).set(0, "colA", "df3");
        report2.getTable(TABLE_NAME).set(0, "colB", 'A');

        report3 = new GATKReport();
        report3.addTable(TABLE_NAME, "Description", 2);
        report3.getTable(TABLE_NAME).addColumn("colA", "%s");
        report3.getTable(TABLE_NAME).addColumn("colB", "%c");
        report3.getTable(TABLE_NAME).set(0, "colA", "df5f");
        report3.getTable(TABLE_NAME).set(0, "colB", 'c');

        report1.concat(report2);
        report1.concat(report3);

        report1.addTable("Table2", "To contain some more data types", 3);
        GATKReportTable table = report1.getTable("Table2");
        table.addColumn("SomeInt", "%d");
        table.addColumn("SomeFloat", "%.16E");
        table.addColumn("TrueFalse", "%B");
        table.addRowIDMapping("12df", 0);
        table.addRowIDMapping("5f", 1);
        table.addRowIDMapping("RZ", 2);
        table.set("12df", "SomeInt", Byte.MAX_VALUE);
        table.set("12df", "SomeFloat", 34.0);
        table.set("12df", "TrueFalse", true);
        table.set("5f", "SomeInt", Short.MAX_VALUE);
        table.set("5f", "SomeFloat", Double.MAX_VALUE);
        table.set("5f", "TrueFalse", false);
        table.set("RZ", "SomeInt", Long.MAX_VALUE);
        table.set("RZ", "SomeFloat", 535646345.657453464576);
        table.set("RZ", "TrueFalse", true);

        report1.addTable("Table3", "blah", 1, GATKReportTable.Sorting.SORT_BY_ROW);
        report1.getTable("Table3").addColumn("a", "%s");
        report1.getTable("Table3").addRowIDMapping("q", 2);
        report1.getTable("Table3").addRowIDMapping("5", 3);
        report1.getTable("Table3").addRowIDMapping("573s", 0);
        report1.getTable("Table3").addRowIDMapping("ZZZ", 1);
        report1.getTable("Table3").set("q", "a", "34");
        report1.getTable("Table3").set("5", "a", "c4g34");
        report1.getTable("Table3").set("573s", "a", "fDlwueg");
        report1.getTable("Table3").set("ZZZ", "a", "Dfs");


        File file = createTempFile("GATKReportGatherer-UnitTest", ".tbl");
        try (final PrintStream ps = new PrintStream(file)) {
            report1.print(ps);
            GATKReport inputRead = new GATKReport(file);
            Assert.assertTrue(report1.isSameFormat(inputRead));
            Assert.assertTrue(report1.equals(inputRead));
        }

    }
}
