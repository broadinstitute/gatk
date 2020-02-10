package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

public class GATKReportTableTest extends GATKBaseTest {

    @Test
    public void testWritingNonFiniteValues() {
        final String TABLE_NAME = "Non-finites";
        final GATKReportTable nonFinites = new GATKReportTable(TABLE_NAME, "Infs, NaNs, and stuff", 3, GATKReportTable.Sorting.DO_NOT_SORT);
        nonFinites.addColumn("name", "");
        nonFinites.addColumn("doubleRepresentation", "%.8f");
        final List<String> names = Arrays.asList("PositiveInfinity","NegativeInfinity","NotANumber");
        final double[] values = {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NaN};
        for (int i = 0; i < values.length; i++) {
            nonFinites.addRowIDMapping(names.get(i), i, true);
            nonFinites.set(i, 1, values[i]);
        }
        final GATKReport report = new GATKReport();
        report.addTable(nonFinites);

        File tableOutput = createTempFile("output.withNonFiniteDoubles", "table");

        try(PrintStream modelReporter = new PrintStream(tableOutput)) {
            report.print(modelReporter);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        final GATKReport reportIn = new GATKReport(tableOutput);
        final GATKReportTable readNonFinites = reportIn.getTable(TABLE_NAME);
        for (final GATKReportColumn reportColumn : readNonFinites.getColumnInfo() ) {
            if (!reportColumn.getColumnName().equals("name")) {
                for (int row = 0; row < readNonFinites.getNumRows(); row++) {
                    if (Double.isNaN(values[row])) {  //NaN isn't equal to NaN, so special case that match
                        Assert.assertTrue(Double.isNaN((Double) readNonFinites.get(row, reportColumn.getColumnName())));
                    } else {
                        Assert.assertTrue(values[row] == (Double) readNonFinites.get(row, reportColumn.getColumnName()));
                    }
                }
            }
        }
    }
}