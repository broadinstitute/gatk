package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link DataLine}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class DataLineUnitTest extends BaseTest {

    @Test(dataProvider = "tableColumnsData")
    public void testCreation(final TableColumnCollection columns) {
        new DataLine(columns, IllegalArgumentException::new);
    }

    @Test(dataProvider = "tableColumnsData")
    public void testToArray(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        Assert.assertEquals(subject.toArray(), new String[columns.columnCount()]);
        Assert.assertNotSame(subject.toArray(), subject.toArray());
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testSetStringByName(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), "" + i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testAppendString(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.append("" + i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testSetInt(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testSetIntUsingIndex(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(i, i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testAppendInt(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.append(i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testSetDoubleByName(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), (double) i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + (double) i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testAppendDouble(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.append((double) i), subject);
        }
        final String[] array = subject.toArray();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(Double.parseDouble(array[i]), (double) i);
        }
    }

    @Test(dataProvider = "tableColumnsData", dependsOnMethods = "testToArray")
    public void testUnpack(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), "" + i), subject);
        }
        final String[] array = subject.unpack();
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(array[i], "" + i);
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testUnpackUnprepared() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.append("a", "b", "c").append(new String[0]).append(1, 2, 3).append(new int[0]).append(1.1);
        subject.unpack();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSetColumnStartingFirstColumn() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.append("a", "b", "c").append(1, 2, 3).append(1.1, 2.2);
        try {
            subject.unpack();
        } catch (final RuntimeException ex) {
            Assert.fail();
        }
        subject.set(0, TableUtils.COMMENT_PREFIX + subject.get(0));

    }

    @Test(dependsOnMethods = "testToArray")
    public void testAppendDiverse() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.append("a", "b", "c").append(new String[0]).append(1, 2, 3).append(new int[0]).append(1.1, 2.2).append(new double[0]);
        Assert.assertEquals(subject.toArray(), new String[]{"a", "b", "c", "1", "2", "3", "" + 1.1, "" + 2.2});
        subject.seek(3).append(4, 5, 6);
        Assert.assertEquals(subject.toArray(), new String[]{"a", "b", "c", "4", "5", "6", "" + 1.1, "" + 2.2});
    }

    @Test(dependsOnMethods = "testToArray", expectedExceptions = IllegalStateException.class)
    public void testAppendBeyondEnd() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, GATKException::new);
        subject.append("a", "b", "c").append(new String[0]).append(1, 2, 3).append(new int[0]).append(1.1, 2.2).append(new double[0]).append("extra");
    }

    @Test(dependsOnMethods = "testToArray")
    public void testSeekByName() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.seek("col1.dbl").append(1.1, 2.2).seek("col1.int").append(1, 2, 3).seek("col1.str").append("a", "b", "c");
        Assert.assertEquals(subject.toArray(), new String[]{"a", "b", "c", "1", "2", "3", "" + 1.1, "" + 2.2});
    }

    private enum TestEnum {
        COL1_STR,
        COL2_STR,
        COL3_STR,
        COL1_INT,
        COL2_INT,
        COL3_INT,
        COL1_DBL,
        COL2_DBL
    }

    @Test(dependsOnMethods = "testToArray")
    public void testSeekByEnum() {
        final TableColumnCollection columns = new TableColumnCollection(TestEnum.class);
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.seek(TestEnum.COL1_DBL).append(1.1, 2.2).seek(TestEnum.COL1_INT).append(1, 2, 3).seek(TestEnum.COL1_STR).append("a", "b", "c");
        Assert.assertEquals(subject.toArray(), new String[]{"a", "b", "c", "1", "2", "3", "" + 1.1, "" + 2.2});
    }

    @Test()
    public void testGetDoubleByEnum() {
        final TableColumnCollection columns = new TableColumnCollection(TestEnum.class);
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(i, (double) i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.getDouble(TestEnum.values()[i]), (double) i, 0.00000001);
        }
    }

    @Test(dependsOnMethods = "testToArray")
    public void testSeekByIndex() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.seek(columns.indexOf("col1.dbl")).append(1.1, 2.2).seek(columns.indexOf("col1.int")).append(1, 2, 3).seek(columns.indexOf("col1.str")).append("a", "b", "c");
        Assert.assertEquals(subject.toArray(), new String[]{"a", "b", "c", "1", "2", "3", "" + 1.1, "" + 2.2});
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSeekToNegative() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.seek(-1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSeekBeyondEnd() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.seek(columns.columnCount() + 1);
    }

    @Test
    public void testSeekToJustBeyondEnd() {
        final TableColumnCollection columns = new TableColumnCollection("col1.str", "col2.str", "col3.str", "col1.int", "col2.int", "col3.int", "col1.dbl", "col2.dbl");
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.seek(columns.columnCount());
        try {
            subject.append("test_value");
            Assert.fail();
        } catch (final IllegalStateException ex) {
            // fine.
        }
    }

    @Test(dataProvider = "tableColumnsData")
    public void testColumns(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        Assert.assertSame(subject.columns(), columns);
    }

    @Test(dataProvider = "tableColumnsData")
    public void testGetInt(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.getInt(i), i);
        }
    }

    @Test(dataProvider = "tableColumnsData")
    public void testGetDouble(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), (double) i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.getDouble(i), (double) i, 0.00000001);
        }
    }

    @Test(dataProvider = "tableColumnsData")
    public void testGetString(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), "" + i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.get(i), "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData")
    public void testGetIntByName(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.getInt(columns.nameAt(i)), i);
        }
    }

    @Test(dataProvider = "tableColumnsData")
    public void testGetDoubleByName(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), (double) i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.getDouble(columns.nameAt(i)), (double) i, 0.00000001);
        }
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = GATKException.class)
    public void testGetNoValidInt(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, GATKException::new);
        subject.set(0, "no-int");
        subject.getInt(0);
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = GATKException.class)
    public void testGetNoValidDouble(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, GATKException::new);
        subject.set(0, "no-dbl");
        subject.getDouble(0);
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalArgumentException.class)
    public void testGetNegativeIndex(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.get(-1);
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalArgumentException.class)
    public void testGetBeyondLastColumn(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.get(columns.columnCount());
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalStateException.class)
    public void testGetUndefinedValue(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.get(columns.nameAt(0));
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalArgumentException.class)
    public void testGetUnknownColumn(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        subject.get("no-column");
    }

    @Test(dataProvider = "tableColumnsData")
    public void testGetStringByName(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertSame(subject.set(columns.nameAt(i), "" + i), subject);
        }
        for (int i = 0; i < columns.columnCount(); i++) {
            Assert.assertEquals(subject.get(columns.nameAt(i)), "" + i);
        }
    }

    @Test(dataProvider = "tableColumnsData")
    public void testSetAll(final TableColumnCollection columns) {
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        final String[] values = new String[columns.columnCount()];
        for (int i = 0; i < values.length; i++) {
            values[i] = "" + i;
        }
        Assert.assertSame(subject.setAll(values), subject);
        Assert.assertEquals(subject.toArray(), values.clone());
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalArgumentException.class)
    public void testSetAllTooShort(final TableColumnCollection columns) {
        // just in case.
        if (columns.columnCount() == 0)
            return;
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        final String[] values = new String[columns.columnCount() - 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = "" + i;
        }
        subject.setAll(values);
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalArgumentException.class)
    public void testSetAllTooLong(final TableColumnCollection columns) {
        // just in case.
        if (columns.columnCount() == 0)
            return;
        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        final String[] values = new String[columns.columnCount() + 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = (i == 0) ?  "#" + i : "" + i;
        }
        subject.setAll(values);
    }

    @Test(dataProvider = "tableColumnsData", expectedExceptions = IllegalArgumentException.class)
    public void testSetAllStartingWithComment(final TableColumnCollection columns) {

        final DataLine subject = new DataLine(columns, IllegalArgumentException::new);
        final String[] values = new String[columns.columnCount() + 1];
        for (int i = 0; i < values.length; i++) {
            values[i] = "" + i;
        }
        subject.setAll(values);
    }

    @DataProvider(name = "tableColumnsData")
    public Object[][] tableColumnsData() {
        return new Object[][]{
                new Object[]{new TableColumnCollection("col1")},
                new Object[]{new TableColumnCollection("col1", "col2", "col3")},
                new Object[]{new TableColumnCollection("col1", "col2", "col3", "col4", "col5", "col6", "col7")}
        };
    }
}
