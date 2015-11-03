package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link TableReader}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TableReaderUnitTest extends BaseTest {

    private static TestTuple[] ORDINARY_VALUE_TEST_TUPLES = new TestTuple[]{
            new TestTuple("str1", 1, 1.1),
            new TestTuple("str2", 2, 2.2E-2),
            new TestTuple("str3", -3, Double.NaN),
            new TestTuple("", Integer.MAX_VALUE, Double.NEGATIVE_INFINITY)
    };

    protected static class TestTuple {
        public final String strValue;
        public final int intValue;
        public final double dblValue;

        public TestTuple(final String s, final int i, final double d) {
            strValue = s;
            intValue = i;
            dblValue = d;
        }

        public String toTabFileLine() {
            return String.join("" + TableUtils.COLUMN_SEPARATOR, strValue, Integer.toString(intValue), Double.toString(dblValue));
        }

        public String toTabFileLineWithAlterInt(final String replace) {
            return String.join("" + TableUtils.COLUMN_SEPARATOR, strValue, replace, Double.toString(dblValue));
        }

        public String toTabFileLineWithAlterDouble(final String replace) {
            return String.join("" + TableUtils.COLUMN_SEPARATOR, strValue, Integer.toString(intValue), replace);
        }

        public int hashCode() {
            return (strValue.hashCode() * 11) + (13 * Integer.hashCode(intValue)) + (31 * Double.hashCode(dblValue));
        }

        public boolean equals(final Object other) {
            if (other instanceof TestTuple) {
                final TestTuple otherTuple = (TestTuple) other;
                return strValue.equals(otherTuple.strValue) &&
                        intValue == otherTuple.intValue &&
                        (dblValue == otherTuple.dblValue ||
                                Double.isNaN(dblValue) == Double.isNaN(otherTuple.dblValue));
            } else {
                return false;
            }
        }
    }

    protected static class TestTupleReader extends TableReader<TestTuple> {

        public TestTupleReader(final File file) throws IOException {
            super(file);
        }

        public TestTupleReader(final FileReader reader) throws IOException {
            super(reader);
        }

        public TestTupleReader(final String sourceName, final FileReader reader) throws IOException {
            super(sourceName, reader);
        }

        @Override
        protected void processColumns(final TableColumnCollection columns) {
            if (columns.columnCount() != 3)
                throw formatException("bad header, must have 3 columns but has " + columns.columnCount() + " instead");
            if (!columns.nameAt(0).equals("col1.str"))
                throw formatException("bad header, first column bad name: " + columns.nameAt(0));
            if (!columns.nameAt(1).equals("col2.int"))
                throw formatException("bad header, second column bad name: " + columns.nameAt(0));
            if (!columns.nameAt(2).equals("col3.dbl"))
                throw formatException("bad header, second column bad name: " + columns.nameAt(0));
        }

        @Override
        protected TestTuple createRecord(final DataLine dataLine) {
            return new TestTuple(
                    dataLine.get("col1.str"),
                    dataLine.getInt("col2.int"),
                    dataLine.getDouble("col3.dbl")
            );
        }
    }

    @Test(dataProvider = "ordinaryValuesData")
    public void testIterator(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TableReader<TestTuple> reader = new TestTupleReader(testFile);
        final List<TestTuple> actual = new ArrayList<>();
        final Iterator<TestTuple> it = reader.iterator();
        while (it.hasNext()) {
            Assert.assertTrue(it.hasNext());
            final TestTuple tuple = it.next();
            Assert.assertNotNull(tuple);
            actual.add(tuple);
        }
        Assert.assertEquals(actual, Arrays.asList(ORDINARY_VALUE_TEST_TUPLES));
        reader.close();
    }

    @Test(dataProvider = "ordinaryValuesData")
    public void testForEach(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TableReader<TestTuple> reader = new TestTupleReader(testFile);
        final List<TestTuple> actual = new ArrayList<>();
        for (final TestTuple tuple : reader) {
            Assert.assertNotNull(tuple);
            actual.add(tuple);
        }
        Assert.assertEquals(actual, Arrays.asList(ORDINARY_VALUE_TEST_TUPLES));
        reader.close();
    }

    @Test(dataProvider = "ordinaryValuesData")
    public void testStream(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TableReader<TestTuple> reader = new TestTupleReader(testFile);
        final List<TestTuple> actual = reader.stream()
                .collect(Collectors.toList());
        Assert.assertEquals(actual, Arrays.asList(ORDINARY_VALUE_TEST_TUPLES));
        reader.close();
    }

    @Test(dataProvider = "ordinaryValuesData")
    public void testStandardValuesUsingReader(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TableReader<TestTuple> reader = new TestTupleReader(new FileReader(testFile));
        // read the expected three records:
        final TestTuple firstRecord = reader.readRecord();
        final TestTuple secondRecord = reader.readRecord();
        final TestTuple thirdRecord = reader.readRecord();
        final TestTuple forthRecord = reader.readRecord();
        Assert.assertNotNull(firstRecord);
        Assert.assertNotNull(secondRecord);
        Assert.assertNotNull(thirdRecord);
        Assert.assertNotNull(forthRecord);
        // beyond end reading:
        final TestTuple fifthRecord = reader.readRecord();
        final TestTuple sixthRecord = reader.readRecord();
        Assert.assertNull(fifthRecord);
        Assert.assertNull(sixthRecord);

        // Check record contents:

        //  "str1", "1", "1.1"
        Assert.assertEquals(firstRecord, ORDINARY_VALUE_TEST_TUPLES[0]);
        Assert.assertEquals(secondRecord, ORDINARY_VALUE_TEST_TUPLES[1]);
        Assert.assertEquals(thirdRecord, ORDINARY_VALUE_TEST_TUPLES[2]);
        Assert.assertEquals(forthRecord, ORDINARY_VALUE_TEST_TUPLES[3]);
    }


    @Test(dataProvider = "ordinaryValuesData")
    public void testStandardValues(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TableReader<TestTuple> reader = new TestTupleReader(testFile);
        // read the expected three records:
        final TestTuple firstRecord = reader.readRecord();
        final TestTuple secondRecord = reader.readRecord();
        final TestTuple thirdRecord = reader.readRecord();
        final TestTuple forthRecord = reader.readRecord();
        Assert.assertNotNull(firstRecord);
        Assert.assertNotNull(secondRecord);
        Assert.assertNotNull(thirdRecord);
        Assert.assertNotNull(forthRecord);
        // beyond end reading:
        final TestTuple fifthRecord = reader.readRecord();
        final TestTuple sixthRecord = reader.readRecord();
        Assert.assertNull(fifthRecord);
        Assert.assertNull(sixthRecord);

        // Check record contents:

        //  "str1", "1", "1.1"
        Assert.assertEquals(firstRecord, ORDINARY_VALUE_TEST_TUPLES[0]);
        Assert.assertEquals(secondRecord, ORDINARY_VALUE_TEST_TUPLES[1]);
        Assert.assertEquals(thirdRecord, ORDINARY_VALUE_TEST_TUPLES[2]);
        Assert.assertEquals(forthRecord, ORDINARY_VALUE_TEST_TUPLES[3]);
    }

    @Test(dataProvider = "commentsOnlyData", expectedExceptions = UserException.BadInput.class)
    public void testEmptyInput(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        new TestTupleReader(testFile);
    }

    @Test(dataProvider = "headerOnlyData")
    public void testHeaderOnly(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TestTupleReader reader = new TestTupleReader(testFile);
        final TestTuple testTuple = reader.readRecord();
        Assert.assertNull(testTuple);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadHeader() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1")
        );
        new TestTupleReader(testFile);
    }

    @Test(dependsOnMethods = "testBadHeader", expectedExceptions = UserException.BadInput.class)
    public void testSourceNameFromFile() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1")
        );
        try {
            new TestTupleReader(testFile);
        } catch (final Throwable tb) {
            Assert.assertTrue(tb.getMessage().matches(".*" + Pattern.quote(testFile.getPath()) + ".*$"));
            throw tb;
        }
    }

    @Test(dependsOnMethods = "testBadHeader", expectedExceptions = UserException.BadInput.class)
    public void testArbitrarySourceName() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1")
        );
        try {
            new TestTupleReader("funny-name", new FileReader(testFile));
        } catch (final Throwable tb) {
            Assert.assertTrue(tb.getMessage().matches(".*" + Pattern.quote("funny-name") + ".*$"), tb.getMessage());
            throw tb;
        }
    }

    @Test
    public void testColumnNames() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );

        final boolean[] columnNamesAvailableWhenCreatingRecords = new boolean[1];


        final TableReader<String[]> reader = new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertNotNull(columns());
                Assert.assertTrue(columns().matchesExactly("col1", "col2", "col3"));
                columnNamesAvailableWhenCreatingRecords[0] = true;
                return dataLine.toArray();
            }
        };

        final TableColumnCollection columns = reader.columns();
        Assert.assertTrue(columns.matchesExactly("col1", "col2", "col3"));

        Assert.assertTrue(columnNamesAvailableWhenCreatingRecords[0], "the readDataLine code did not get executed");
    }

    @Test
    public void testReadFromString() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3")
        );

        final TableReader<String[]> reader = new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertNotNull(columns());
                return dataLine.toArray();
            }
        };

        final String line = String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3");
        final String[] record = reader.readRecord(line);
        Assert.assertEquals(record, new String[]{"1", "2", "3"});
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadFromStringTooManyValues() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3")
        );

        final TableReader<String[]> reader = new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertNotNull(columns());
                return dataLine.toArray();
            }
        };

        final String line = String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3", "4");
        reader.readRecord(line);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadFromStringTooFewValues() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3")
        );

        final TableReader<String[]> reader = new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertNotNull(columns());
                return dataLine.toArray();
            }
        };

        final String line = String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2");
        reader.readRecord(line);
    }

    @Test
    public void testColumnIndex() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );

        final boolean[] tested = new boolean[1];

        new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertNotNull(columns());
                Assert.assertTrue(columns().matchesAll(0, "col1", "col2", "col3"));
                Assert.assertEquals(columns().indexOf("no-col"), -1);
                tested[0] = true;
                return dataLine.toArray();
            }
        };

        Assert.assertTrue(tested[0], "the readDataLine code did not get executed");
    }

    @Test
    public void testColumnValueAsInt() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );

        final boolean[] tested = new boolean[1];

        new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertEquals(dataLine.getInt("col1"), 1);
                Assert.assertEquals(dataLine.getInt("col2"), 2);
                Assert.assertEquals(dataLine.getInt("col3"), 3);
                try {
                    dataLine.getInt("no-col");
                    Assert.fail();
                } catch (final IllegalArgumentException ex) {
                    // expected.
                } catch (final RuntimeException ex) {
                    Assert.fail();
                }
                tested[0] = true;
                return dataLine.toArray();
            }
        };

        Assert.assertTrue(tested[0], "the readDataLine code did not get exected");
    }

    @Test
    public void testColumnValueAsDouble() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );

        final boolean[] tested = new boolean[1];

        new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertEquals(dataLine.getDouble("col1"), 1.0);
                Assert.assertEquals(dataLine.getDouble("col2"), 2.0);
                Assert.assertEquals(dataLine.getDouble("col3"), 3.0);
                try {
                    dataLine.getDouble("no-col");
                    Assert.fail();
                } catch (final IllegalArgumentException ex) {
                    // expected.
                } catch (final RuntimeException ex) {
                    Assert.fail();
                }
                tested[0] = true;
                return dataLine.toArray();
            }
        };

        Assert.assertTrue(tested[0], "the readDataLine code did not get exected");
    }

    @Test
    public void testColumnValueAsString() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );

        final boolean[] tested = new boolean[1];

        new TableReader<String[]>(testFile) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                Assert.assertEquals(dataLine.get("col1"), "1");
                Assert.assertEquals(dataLine.get("col2"), "2");
                Assert.assertEquals(dataLine.get("col3"), "3");
                try {
                    dataLine.get("no-col");
                    Assert.fail();
                } catch (final IllegalArgumentException ex) {
                    // expected.
                } catch (final RuntimeException ex) {
                    Assert.fail();
                }
                tested[0] = true;
                return dataLine.toArray();
            }
        };

        Assert.assertTrue(tested[0], "the readDataLine code did not get exected");
    }


    @Test(dependsOnMethods = "testBadHeader", expectedExceptions = UserException.BadInput.class)
    public void testNoSourceName() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1")
        );
        try {
            new TestTupleReader(null, new FileReader(testFile));
        } catch (final Throwable tb) {
            Assert.assertTrue(tb.getMessage().matches("^.*format error at line.*$"), tb.getMessage());
            throw tb;
        }
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testTooManyValues() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3", "4"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );
        final TableReader<String> reader = firstOfThreeColumnReader(testFile);
        final String firstTestTuple = reader.readRecord();
        Assert.assertNotNull(firstTestTuple);
        reader.readRecord();  // will cause the exception.
    }

    /**
     * Creates a silly reader use for testing somewhere else that
     * returns the first of value out of three columns.
     *
     * @param testFile
     * @return never {@code null}.
     * @throws IOException
     */
    private TableReader<String> firstOfThreeColumnReader(final File testFile) throws IOException {
        return new TableReader<String>(testFile) {

            @Override
            protected String createRecord(final DataLine dataLine) {
                return dataLine.get(0);
            }
        };
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testTooFewValues() throws IOException {
        final File testFile = createTestInput(
                String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2"),
                String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3")
        );
        final TableReader<String> reader = firstOfThreeColumnReader(testFile);
        final String firstTestTuple = reader.readRecord();
        Assert.assertNotNull(firstTestTuple);
        reader.readRecord(); // will cause the exception.
    }

    @Test(dataProvider = "dataTypeConversionErrorData", expectedExceptions = UserException.BadInput.class)
    public void testDataTypeConversionError(final String[] lines) throws IOException {
        final File testFile = createTestInput(lines);
        final TestTupleReader reader = new TestTupleReader(testFile);
        while (reader.readRecord() != null) {
            // eventually will cause the exception.
        }
    }

    @Test
    public void testIgnoreHeaderRepetitions() throws IOException {
        final File testFile = createTestInput(
                TableUtils.COMMENT_PREFIX + "comment1",
                "Column1",
                "",
                TableUtils.COMMENT_PREFIX + "comment2",
                "abc",
                "Column1",
                "cbb"
        );
        final TableReader<String> tableReader = new TableReader<String>(testFile) {

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                Assert.assertEquals(columns.columnCount(), 1);
                Assert.assertEquals(columns.nameAt(0), "Column1");
            }

            @Override
            protected String createRecord(final DataLine dataLine) {
                return dataLine.get(0);
            }
        };

        final String record = tableReader.readRecord();
        final String secondRecord = tableReader.readRecord();
        final String thirdRecord = tableReader.readRecord();
        final String forthRecord = tableReader.readRecord();
        Assert.assertNotNull(record);
        Assert.assertEquals(record, "");
        Assert.assertEquals(thirdRecord, "cbb");
        Assert.assertNotNull(secondRecord,"abc");
        Assert.assertNull(forthRecord);
    }

    @Test
    public void testSingleColumnWithEmptyDataLines() throws IOException {
        final File testFile = createTestInput(
                TableUtils.COMMENT_PREFIX + "comment1",
                "Column1",
                "",
                TableUtils.COMMENT_PREFIX + "comment2"
        );
        final TableReader<String> tableReader = new TableReader<String>(testFile) {

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                Assert.assertEquals(columns.columnCount(), 1);
                Assert.assertEquals(columns.nameAt(0), "Column1");
            }

            @Override
            protected String createRecord(final DataLine dataLine) {
                return dataLine.get(0);
            }
        };

        final String record = tableReader.readRecord();
        final String secondRecord = tableReader.readRecord();
        final String thirdRecord = tableReader.readRecord();
        Assert.assertNotNull(record);
        Assert.assertEquals(record, "");
        Assert.assertNull(thirdRecord);
        Assert.assertNull(secondRecord);
    }

    private File createTestInput(final String... lines) throws IOException {
        final File testFile = createTempFile("test", ".tab");
        final PrintWriter testWriter = new PrintWriter(new FileWriter(testFile));
        for (final String line : lines) {
            testWriter.println(line);
        }
        testWriter.close();
        return testFile;
    }

    @DataProvider(name = "ordinaryValuesData")
    public Object[][] ordinaryValuesData() {
        final String comment1 = TableUtils.COMMENT_PREFIX + "comment1";
        final String comment2 = TableUtils.COMMENT_PREFIX + "comment2";
        final String header = String.join("" + TableUtils.COLUMN_SEPARATOR, "col1.str", "col2.int", "col3.dbl");
        final String comment3 = TableUtils.COMMENT_PREFIX + "comment3";
        final String dataLine1 = ORDINARY_VALUE_TEST_TUPLES[0].toTabFileLine();
        final String comment4 = TableUtils.COMMENT_PREFIX + "comment4";
        final String comment5 = TableUtils.COMMENT_PREFIX + "comment5";
        final String dataLine2 = ORDINARY_VALUE_TEST_TUPLES[1].toTabFileLine();
        final String dataLine3 = ORDINARY_VALUE_TEST_TUPLES[2].toTabFileLine();
        final String comment6 = TableUtils.COMMENT_PREFIX + "comment6";
        final String dataLine4 = ORDINARY_VALUE_TEST_TUPLES[3].toTabFileLine();
        return new Object[][]{
                {new String[]{comment1, comment2, header, comment3, dataLine1, comment4, comment5, dataLine2, dataLine3, comment6, dataLine4}},
                {new String[]{header, dataLine1, dataLine2, dataLine3, dataLine4}},
                {new String[]{header, dataLine1, dataLine2, dataLine3, dataLine4, comment1, comment2}}
        };
    }

    @DataProvider(name = "commentsOnlyData")
    public Object[][] commentsOnlyData() {
        final String comment1 = TableUtils.COMMENT_PREFIX + "comment1";
        final String comment2 = TableUtils.COMMENT_PREFIX + "comment2";
        final String comment3 = TableUtils.COMMENT_PREFIX + "comment3";
        final String comment4 = TableUtils.COMMENT_PREFIX + "comment4";
        final String comment5 = TableUtils.COMMENT_PREFIX + "comment5";
        final String comment6 = TableUtils.COMMENT_PREFIX + "comment6";
        return new Object[][]{
                {new String[]{comment1, comment2, comment3, comment4, comment5, comment6}},
                {new String[]{comment1, comment2}},
                {new String[]{}},
                {new String[]{comment1}}
        };
    }

    @DataProvider(name = "dataTypeConversionErrorData")
    public Object[][] dataTypeConversionErrorData() {
        final String header = String.join("" + TableUtils.COLUMN_SEPARATOR, "col1.str", "col2.int", "col3.dbl");
        final String dataLine1 = ORDINARY_VALUE_TEST_TUPLES[0].toTabFileLine();
        final String dataLine2 = ORDINARY_VALUE_TEST_TUPLES[1].toTabFileLine();
        final String dataLine2withBadInt = ORDINARY_VALUE_TEST_TUPLES[1].toTabFileLineWithAlterInt("no-int");
        final String dataLine2withBadInt2 = ORDINARY_VALUE_TEST_TUPLES[1].toTabFileLineWithAlterInt("2.2");
        final String dataLine2withBadDouble = ORDINARY_VALUE_TEST_TUPLES[1].toTabFileLineWithAlterDouble("a121");
        return new Object[][]{
                {new String[]{header, dataLine1, dataLine2withBadDouble}},
                {new String[]{header, dataLine1, dataLine2withBadInt}},
                {new String[]{header, dataLine1, dataLine2withBadInt2}},
                {new String[]{header, dataLine2withBadDouble}}
        };
    }

    @DataProvider(name = "headerOnlyData")
    public Object[][] headerOnlyData() {
        final String header = String.join("" + TableUtils.COLUMN_SEPARATOR, "col1.str", "col2.int", "col3.dbl");
        final String comment1 = TableUtils.COMMENT_PREFIX + "comment1";
        final String comment2 = TableUtils.COMMENT_PREFIX + "comment2";
        final String comment3 = TableUtils.COMMENT_PREFIX + "comment3";
        final String comment4 = TableUtils.COMMENT_PREFIX + "comment4";
        final String comment5 = TableUtils.COMMENT_PREFIX + "comment5";
        final String comment6 = TableUtils.COMMENT_PREFIX + "comment6";
        return new Object[][]{
                {new String[]{comment1, comment2, comment3, header, comment4, comment5, comment6}},
                {new String[]{header, comment1, comment2}},
                {new String[]{header}},
                {new String[]{comment1, header}}
        };
    }


}
