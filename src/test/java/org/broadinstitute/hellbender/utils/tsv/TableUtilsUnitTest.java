package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Objects;

/**
 * Unit tests for {@link TableUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TableUtilsUnitTest extends GATKBaseTest {

    private static File CORRECT_TEST_FILE;

    private static File INVALID_HEADER_FILE;

    private static File INVALID_RECORD_FILE;

    @BeforeClass
    public void setUp() throws IOException {
        CORRECT_TEST_FILE = createTestFile();
        INVALID_HEADER_FILE = createInvalidHeaderTestFile();
        INVALID_RECORD_FILE = createInvalidRecordHeaderTestFile();
    }

    @AfterClass
    public void tearDown() throws IOException {
        CORRECT_TEST_FILE.delete();
        INVALID_HEADER_FILE.delete();
        INVALID_RECORD_FILE.delete();
    }

    @Test(dataProvider = "correctFileReaders")
    public void testReader(final TableReader<TestTuple> reader) {
        final TestTuple[] actualTuples = reader.stream().toArray(TestTuple[]::new);
        Assert.assertEquals(actualTuples, TEST_RECORD);
    }

    @Test(dependsOnMethods = "testReader")
    public void testWriter0() throws IOException {
        final File testFile = createTempFile("test",".tsv");
        final TableWriter<TestTuple> writer = TableUtils.writer(testFile,TEST_COLUMNS,
                (tuple,dataLine) -> {
                    dataLine.append(tuple.strValue).append(tuple.intValue).append(tuple.dblValue);
                });
        for (final TestTuple tuple : TEST_RECORD) {
            writer.writeRecord(tuple);
        }
        writer.close();
        final TableReader<TestTuple> reader = TableUtils.reader(testFile,(columns, exceptionFactory) -> {
            if (!columns.matchesExactly("col1.str","col2.int","col3.dbl"))
                Assert.fail();
            return (dataLine) -> new TestTuple(dataLine.get(0),dataLine.getInt(1),dataLine.getDouble(2));
        });
        final TestTuple[] actualTuples = reader.stream().toArray(TestTuple[]::new);
        Assert.assertEquals(actualTuples, TEST_RECORD);
        reader.close();
    }

    @Test(dependsOnMethods = "testReader")
    public void testWriter1() throws IOException {
        final File testFile = createTempFile("test",".tsv");
        final TableWriter<TestTuple> writer = TableUtils.writer(new FileWriter(testFile),TEST_COLUMNS,
                (tuple,dataLine) -> {
                    dataLine.append(tuple.strValue).append(tuple.intValue).append(tuple.dblValue);
                });
        for (final TestTuple tuple : TEST_RECORD) {
            writer.writeRecord(tuple);
        }
        writer.close();
        final TableReader<TestTuple> reader = TableUtils.reader(testFile,(columns, exceptionFactory) -> {
            if (!columns.matchesExactly("col1.str","col2.int","col3.dbl"))
                Assert.fail();
            return (dataLine) -> new TestTuple(dataLine.get(0),dataLine.getInt(1),dataLine.getDouble(2));
        });
        final TestTuple[] actualTuples = reader.stream().toArray(TestTuple[]::new);
        Assert.assertEquals(actualTuples, TEST_RECORD);
        reader.close();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReaderInvalidFactory0() throws IOException {
        TableUtils.reader(CORRECT_TEST_FILE, null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReaderInvalidFactory1() throws IOException {
        TableUtils.reader(new FileReader(CORRECT_TEST_FILE), null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReaderInvalidFactory2() throws IOException {
        TableUtils.reader("some-name", new FileReader(CORRECT_TEST_FILE), null);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testReaderInvalidFactoryStateNull0() throws IOException {
        TableUtils.reader(CORRECT_TEST_FILE, (c, e) -> null );
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testReaderInvalidFactoryReturnsNull1() throws IOException {
        TableUtils.reader(new FileReader(CORRECT_TEST_FILE),  (c, e) -> null);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testReaderInvalidFactoryReturnsNull2() throws IOException {
        TableUtils.reader("some-name",new FileReader(CORRECT_TEST_FILE), (c, e) -> null);
    }

    @Test(dataProvider = "invalidRecordReaders", expectedExceptions = UserException.BadInput.class)
    public void testInvalidRecordReader(final TableReader<TestTuple> reader) {
        final TestTuple[] actualTuples = reader.stream().toArray(TestTuple[]::new);
        Assert.assertEquals(actualTuples, TEST_RECORD);
    }



    @Test(expectedExceptions = UserException.BadInput.class)
    public void testInvalidHeaderReader0() throws IOException {
        TableUtils.<TestTuple>reader(INVALID_HEADER_FILE,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("col1.str", "col2.int", "col3.dbl"))
                        throw formatExceptionFactory.apply("Bad header");
                    return (dataLine) -> new TestTuple(dataLine.get(0), dataLine.getInt(1), dataLine.getDouble(2));
                });
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testInvalidHeaderReader1() throws IOException {
        TableUtils.<TestTuple>reader(new FileReader(INVALID_HEADER_FILE),
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("col1.str", "col2.int", "col3.dbl"))
                        throw formatExceptionFactory.apply("Bad header");
                    return (dataLine) -> new TestTuple(dataLine.get(0), dataLine.getInt(1), dataLine.getDouble(2));
                });
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testInvalidHeaderReader2() throws IOException {
        TableUtils.<TestTuple>reader("source-name",new FileReader(INVALID_HEADER_FILE),
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("col1.str", "col2.int", "col3.dbl"))
                        throw formatExceptionFactory.apply("Bad header");
                    return (dataLine) -> new TestTuple(dataLine.get(0), dataLine.getInt(1), dataLine.getDouble(2));
                });
    }



    @DataProvider(name = "correctFileReaders")
    public Object[][] correctFileReaders() throws IOException {
        return readers(CORRECT_TEST_FILE);
    }

    @DataProvider(name = "invalidRecordReaders")
    public Object[][] invalidRecordReaders() throws IOException {
        return readers(INVALID_RECORD_FILE);
    }


    private Object[][] readers(final File inputFile) throws IOException {
        return new Object[][]{
                {
                     TableUtils.<TestTuple>reader(inputFile,
                        (columns, formatExceptionFactory) -> {
                            if (!columns.matchesExactly("col1.str", "col2.int", "col3.dbl"))
                                throw formatExceptionFactory.apply("Bad header");
                            return (dataLine) -> new TestTuple(dataLine.get(0), dataLine.getInt(1), dataLine.getDouble(2));
                        })
                }, {
                    TableUtils.<TestTuple>reader(new FileReader(inputFile),
                        (columns, formatExceptionFactory) -> {
                            if (!columns.matchesExactly("col1.str", "col2.int", "col3.dbl"))
                                throw formatExceptionFactory.apply("Bad header");
                            return (dataLine) -> new TestTuple(dataLine.get(0), dataLine.getInt(1), dataLine.getDouble(2));
                        })
                }, {
                    TableUtils.<TestTuple>reader("source-name",new FileReader(inputFile),
                        (columns, formatExceptionFactory) -> {
                            if (!columns.matchesExactly("col1.str", "col2.int", "col3.dbl"))
                                throw formatExceptionFactory.apply("Bad header");
                            return (dataLine) -> new TestTuple(dataLine.get(0), dataLine.getInt(1), dataLine.getDouble(2));
                        })

                }
        };
    }

    private final static class TestTuple {
        public final String strValue;
        public final int intValue;
        public final double dblValue;

        private TestTuple(final String strValue, final int intValue, final double dblValue) {
            this.strValue = strValue;
            this.intValue = intValue;
            this.dblValue = dblValue;
        }

        private String toFileString() {
            return String.join(TableUtils.COLUMN_SEPARATOR_STRING, strValue, Integer.toString(intValue), Double.toString(dblValue));
        }

        public boolean equals(final Object other) {
            if (!(other instanceof TestTuple))
                return false;
            else {
                final TestTuple otherCasted = (TestTuple) other;
                return Objects.equals(this.strValue,otherCasted.strValue) && this.intValue == otherCasted.intValue
                        && Math.abs(this.dblValue - otherCasted.dblValue) < 0.000001;
            }
        }

        public int hashCode() {
            return strValue.hashCode();
        }
    }

    private final static TableColumnCollection TEST_COLUMNS =
            new TableColumnCollection("col1.str","col2.int","col3.dbl");

    private final static TestTuple[] TEST_RECORD = new TestTuple[] {
        new TestTuple("1",1,1.1),
        new TestTuple("2",2,2.2),
        new TestTuple("3",3,2.2),
    };

    public File createTestFile() throws IOException {
        final File result = createTempFile("test",".tsv");
        final PrintWriter writer = new PrintWriter(new FileWriter(result));
        writer.println("#comment1");
        writer.println("#comment2");
        writer.println(String.join(TableUtils.COLUMN_SEPARATOR_STRING,TEST_COLUMNS.names().stream().toArray(i -> new String[i])));
        writer.println("#comment3");
        writer.println(TEST_RECORD[0].toFileString());
        writer.println(TEST_RECORD[1].toFileString());
        writer.println("#comment4");
        writer.println(TEST_RECORD[2].toFileString());
        writer.println("#comment5");
        writer.close();
        return result;
    }

    private File createInvalidHeaderTestFile() throws IOException {
        final File result = createTempFile("test",".tsv");
        final PrintWriter writer = new PrintWriter(new FileWriter(result));
        writer.println("#comment1");
        writer.println("#comment2");
        writer.println(String.join(TableUtils.COLUMN_SEPARATOR_STRING,"col1.str","col2.str"));
        writer.println("#comment3");
        writer.println(TEST_RECORD[0].toFileString());
        writer.println(TEST_RECORD[1].toFileString());
        writer.println("#comment4");
        writer.println(TEST_RECORD[2].toFileString());
        writer.println("#comment5");
        writer.close();
        return result;
    }

    private File createInvalidRecordHeaderTestFile() throws IOException {
        final File result = createTempFile("test",".tsv");
        final PrintWriter writer = new PrintWriter(new FileWriter(result));
        writer.println("#comment1");
        writer.println("#comment2");
        writer.println(String.join(TableUtils.COLUMN_SEPARATOR_STRING,"col1.str","col2.int","col3.dbl"));
        writer.println("#comment3");
        writer.println(TEST_RECORD[0].toFileString());
        writer.println(TEST_RECORD[1].toFileString());
        writer.println(String.join(TableUtils.COLUMN_SEPARATOR_STRING,"10","aaas","10.1"));
        writer.println("#comment4");
        writer.println(TEST_RECORD[2].toFileString());
        writer.println("#comment5");
        writer.close();
        return result;
    }


}
