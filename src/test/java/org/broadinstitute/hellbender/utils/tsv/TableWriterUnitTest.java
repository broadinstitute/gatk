package org.broadinstitute.hellbender.utils.tsv;

import java.nio.file.Files;
import java.nio.file.Path;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link TableWriter}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TableWriterUnitTest extends GATKBaseTest {

    @Test
    public void testNoRecordsOutput() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<Object> writer = new TableWriter<Object>(testPath, columnNames) {
            @Override
            protected void composeLine(Object record, final DataLine dataLine) {
                dataLine.append(record.toString());
            }
        };
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 1);
        Assert.assertEquals(outLines.get(0), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
    }

    @Test
    public void testNoRecordsOutputWithCommentsAndDuplicatedAttemptsToWriterHeader() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<Object> writer = new TableWriter<Object>(testPath, columnNames) {
            @Override
            protected void composeLine(Object record, final DataLine dataLine) {
                dataLine.append(record.toString());
            }
        };
        writer.writeComment("commentLine1");
        writer.writeComment("commentLine2");
        writer.writeHeaderIfApplies();
        writer.writeComment("commentLine3");
        writer.writeHeaderIfApplies();
        writer.writeComment("commentLine4");
        writer.writeHeaderIfApplies();
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 5, outLines.toString());
        Assert.assertEquals(outLines.get(0), TableUtils.COMMENT_PREFIX + "commentLine1");
        Assert.assertEquals(outLines.get(1), TableUtils.COMMENT_PREFIX + "commentLine2");
        Assert.assertEquals(outLines.get(2), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(3), TableUtils.COMMENT_PREFIX + "commentLine3");
        Assert.assertEquals(outLines.get(4), TableUtils.COMMENT_PREFIX + "commentLine4");
    }

    @Test
    public void testNoRecordsOutputWithComments() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<Object> writer = new TableWriter<Object>(testPath, columnNames) {
            @Override
            protected void composeLine(Object record, final DataLine dataLine) {
                dataLine.append(record.toString());
            }
        };
        writer.writeComment("commentLine1");
        writer.writeComment("commentLine2");
        writer.writeHeaderIfApplies();
        writer.writeComment("commentLine3");
        writer.writeComment("commentLine4");
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 5, outLines.toString());
        Assert.assertEquals(outLines.get(0), TableUtils.COMMENT_PREFIX + "commentLine1");
        Assert.assertEquals(outLines.get(1), TableUtils.COMMENT_PREFIX + "commentLine2");
        Assert.assertEquals(outLines.get(2), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(3), TableUtils.COMMENT_PREFIX + "commentLine3");
        Assert.assertEquals(outLines.get(4), TableUtils.COMMENT_PREFIX + "commentLine4");
    }

    @Test
    public void testSomeRecordsOutputWithCommentsAndLateAttemptToWriteHeader() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<String[]> writer = new TableWriter<String[]>(testPath, columnNames) {
            @Override
            protected void composeLine(final String[] record, final DataLine dataLine) {
                dataLine.setAll(record);
            }
        };
        writer.writeComment("commentLine1");
        writer.writeComment("commentLine2");
        writer.writeRecord(new String[]{"1", "2", "3"});
        writer.writeHeaderIfApplies();
        writer.writeRecord(new String[]{"4", "5", "6"});
        writer.writeHeaderIfApplies();
        writer.writeComment("commentLine3");
        writer.writeHeaderIfApplies();
        writer.writeRecord(new String[]{"-1", "-2", "-3"});
        writer.writeHeaderIfApplies();
        writer.writeComment("commentLine4");
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 8, outLines.toString());
        Assert.assertEquals(outLines.get(0), TableUtils.COMMENT_PREFIX + "commentLine1");
        Assert.assertEquals(outLines.get(1), TableUtils.COMMENT_PREFIX + "commentLine2");
        Assert.assertEquals(outLines.get(2), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(3), String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"));
        Assert.assertEquals(outLines.get(4), String.join("" + TableUtils.COLUMN_SEPARATOR, "4", "5", "6"));
        Assert.assertEquals(outLines.get(5), TableUtils.COMMENT_PREFIX + "commentLine3");
        Assert.assertEquals(outLines.get(6), String.join("" + TableUtils.COLUMN_SEPARATOR, "-1", "-2", "-3"));
        Assert.assertEquals(outLines.get(7), TableUtils.COMMENT_PREFIX + "commentLine4");
    }

    @Test
    public void testSomeRecordsOutputWithComments() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<String[]> writer = new TableWriter<String[]>(testPath, columnNames) {
            @Override
            protected void composeLine(final String[] record, final DataLine dataLine) {
                dataLine.setAll(record);
            }
        };
        writer.writeComment("commentLine1");
        writer.writeComment("commentLine2");
        writer.writeRecord(new String[]{"1", "2", "3"});
        writer.writeRecord(new String[]{"4", "5", "6"});
        writer.writeComment("commentLine3");
        writer.writeRecord(new String[]{"-1", "-2", "-3"});
        writer.writeComment("commentLine4");
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 8, outLines.toString());
        Assert.assertEquals(outLines.get(0), TableUtils.COMMENT_PREFIX + "commentLine1");
        Assert.assertEquals(outLines.get(1), TableUtils.COMMENT_PREFIX + "commentLine2");
        Assert.assertEquals(outLines.get(2), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(3), String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"));
        Assert.assertEquals(outLines.get(4), String.join("" + TableUtils.COLUMN_SEPARATOR, "4", "5", "6"));
        Assert.assertEquals(outLines.get(5), TableUtils.COMMENT_PREFIX + "commentLine3");
        Assert.assertEquals(outLines.get(6), String.join("" + TableUtils.COLUMN_SEPARATOR, "-1", "-2", "-3"));
        Assert.assertEquals(outLines.get(7), TableUtils.COMMENT_PREFIX + "commentLine4");
    }

    @Test
    public void testMetadata() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<String[]> writer = new TableWriter<String[]>(testPath, columnNames) {
            @Override
            protected void composeLine(final String[] record, final DataLine dataLine) {
                dataLine.setAll(record);
            }
        };

        final String key1 = "key1";
        final String value1 = "val1";
        final String key2 = "key2";
        final String value2 = "val2";

        writer.writeMetadata(key1, value1);
        writer.writeComment("commentLine1");
        writer.writeMetadata(key2, value2);
        writer.writeRecord(new String[]{"1", "2", "3"});
        writer.writeComment("commentLine2");
        writer.writeRecord(new String[]{"4", "5", "6"});
        writer.writeComment("commentLine3");
        writer.writeRecord(new String[]{"-1", "-2", "-3"});
        writer.writeComment("commentLine4");
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.get(0), TableUtils.COMMENT_PREFIX + TableWriter.METADATA_TAG + key1 + "=" + value1);
        Assert.assertEquals(outLines.get(1), TableUtils.COMMENT_PREFIX + "commentLine1");
        Assert.assertEquals(outLines.get(2), TableUtils.COMMENT_PREFIX + TableWriter.METADATA_TAG + key2 + "=" + value2);
        Assert.assertEquals(outLines.get(3), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(4), String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"));

    }

    @Test
    public void testWriteRecordFromIterable() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<String[]> writer = new TableWriter<String[]>(testPath, columnNames) {
            @Override
            protected void composeLine(final String[] record, final DataLine dataLine) {
                dataLine.setAll(record);
            }
        };
        final List<String[]> records = new ArrayList<>(3);
        records.add(new String[]{"1", "2", "3"});
        records.add(new String[]{"4", "5", "6"});
        records.add(new String[]{"-1", "-2", "-3"});
        writer.writeAllRecords(records);
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 4, outLines.toString());
        Assert.assertEquals(outLines.get(0), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(1), String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"));
        Assert.assertEquals(outLines.get(2), String.join("" + TableUtils.COLUMN_SEPARATOR, "4", "5", "6"));
        Assert.assertEquals(outLines.get(3), String.join("" + TableUtils.COLUMN_SEPARATOR, "-1", "-2", "-3"));
    }

    @Test
    public void testWriteRecordFromArray() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableColumnCollection columnNames = new TableColumnCollection("col1", "col2", "col3");
        final TableWriter<String[]> writer = new TableWriter<String[]>(testPath, columnNames) {
            @Override
            protected void composeLine(final String[] record, final DataLine dataLine) {
                dataLine.setAll(record);
            }
        };
        final List<String[]> records = new ArrayList<>(3);
        records.add(new String[]{"1", "2", "3"});
        records.add(new String[]{"4", "5", "6"});
        records.add(new String[]{"-1", "-2", "-3"});
        writer.writeAllRecords(records);
        writer.close();

        final List<String> outLines = outputLines(testPath);
        Assert.assertEquals(outLines.size(), 4, outLines.toString());
        Assert.assertEquals(outLines.get(0), String.join("" + TableUtils.COLUMN_SEPARATOR, "col1", "col2", "col3"));
        Assert.assertEquals(outLines.get(1), String.join("" + TableUtils.COLUMN_SEPARATOR, "1", "2", "3"));
        Assert.assertEquals(outLines.get(2), String.join("" + TableUtils.COLUMN_SEPARATOR, "4", "5", "6"));
        Assert.assertEquals(outLines.get(3), String.join("" + TableUtils.COLUMN_SEPARATOR, "-1", "-2", "-3"));
    }

    @Test
    public void testTestTupleWriting() throws IOException {
        final Path testPath = createTempPath("test", ".tab");
        final TableWriter<TableReaderUnitTest.TestTuple>
                writer = new TableWriter<TableReaderUnitTest.TestTuple>(testPath, new TableColumnCollection("col1.str", "col2.int", "col3.dbl")) {
            @Override
            protected void composeLine(final TableReaderUnitTest.TestTuple record, final DataLine dataLine) {
                dataLine.set("col1.str", record.strValue)
                        .set("col2.int", record.intValue)
                        .set("col3.dbl", record.dblValue);
            }
        };
    }

    private List<String> outputLines(Path testPath) throws IOException {
        try (BufferedReader bufferedReader = Files.newBufferedReader(testPath)) {
            return bufferedReader.lines().collect(Collectors.toList());
        }
    }
}
