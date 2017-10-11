package org.broadinstitute.hellbender.utils.test;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountFileHeaderKey;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * This is a testing utility class to check the correctness of coverage collection results for CNV analysis.
 * It assumes that results come in a TSV format, and checks for exact string matches between two tables. Therefore,
 * this class is NOT intended to be used with output formats containing floats, as it doesn't allow for rounding error.
 */
public class IntegerReadCountFileComparator {
    private final List<String[]> values;
    private final int columnCount;
    private final int rowCount;
    private final String[] columnNames;
    private final Map<String, String> headerValuesMap;

    private final static String TAB_SEPARATOR_REGEX = "\\t";
    private final static String COMMENT_PREFIX_START_REGEX = "^#.*$";

    /**
     * Private constructor
     *
     * @param file input file
     * @throws IOException if there is any problems while reading the file
     * @throws UserException.BadInput if non-integer value is encountered below column names line or file
     * is not formatted correctly
     */
    private IntegerReadCountFileComparator(final File file) throws IOException {
        Utils.nonNull(file, "the file cannot be null");
        final BufferedReader lineReader = new BufferedReader(new FileReader(file));
        headerValuesMap = new HashMap<>(10);
        String lastLine;

        //read in the header values
        String headerKeyValuePattern = ReadCountFileHeaderKey.headerKeyValuePattern;
        while ((lastLine = lineReader.readLine()) != null) {
            if (lastLine.matches(COMMENT_PREFIX_START_REGEX)) {
                headerValuesMap.put(lastLine.replaceAll(headerKeyValuePattern, "$1"), lastLine.replaceAll(headerKeyValuePattern, "$2"));
            } else {
                break;
            }
        }

        //read in table's column values
        columnNames = lastLine == null ? null : lineReader.readLine().split(TAB_SEPARATOR_REGEX);
        Utils.nonNull(columnNames, "the table does not have a header");
        columnCount = columnNames.length;

        //read in table's values
        values = new ArrayList<>();
        int currentRowCount = 0;
        while ((lastLine = lineReader.readLine()) != null) {
            final String[] lineValues = lastLine.split(TAB_SEPARATOR_REGEX);
            if (lineValues.length != columnCount) {
                throw new UserException.BadInput(String.format("line %d of the table has different number of values than the column name line", currentRowCount + 1));
            }
            values.add(new String[columnCount]);
            for (int index = 0; index < columnCount; index ++) {
                values.get(currentRowCount)[index] = lineValues[index];

            }
            currentRowCount++;
        }

        rowCount = currentRowCount;
        lineReader.close();
    }

    /**
     * Check if two integer read counts tables are equivalent. In particular we check that all the table values are identical,
     * and that the headers have identical lists of key value pairs.
     *
     * @param f1 first TSV file to compare
     * @param f2 second TSV file to compare
     * @param headerKeysToCompare list of keys to check in both table headers
     */
    public static void AssertEquals(final File f1, final File f2, final List<String> headerKeysToCompare) {
        try {
            IntegerReadCountFileComparator firstTable = new IntegerReadCountFileComparator(f1);
            IntegerReadCountFileComparator secondTable = new IntegerReadCountFileComparator(f2);
            //check that mandatory header key value pairs are the same.
            if (headerKeysToCompare != null) {
                headerKeysToCompare.stream().forEach(key -> {
                    Assert.assertTrue(firstTable.headerValuesMap.containsKey(key));
                    Assert.assertTrue(secondTable.headerValuesMap.containsKey(key));
                    Assert.assertEquals(firstTable.headerValuesMap.get(key), secondTable.headerValuesMap.get(key));
                });
            }
            //check that tables have the same dimensions and columns
            Assert.assertEquals(firstTable.columnCount, secondTable.columnCount, "Different number of columns");
            Assert.assertEquals(firstTable.rowCount, secondTable.rowCount, "Different number of rows");
            Assert.assertEquals(firstTable.columnNames, secondTable.columnNames, "Different column names");
            //check that tables have same values row by row
            IntStream.range(0, firstTable.rowCount).forEach(row -> Assert.assertEquals(firstTable.values.get(row), secondTable.values.get(row)));
        } catch (final IOException ex) {
            Assert.fail("Failed to read actual or expected table files ", ex);
        }
    }
}
