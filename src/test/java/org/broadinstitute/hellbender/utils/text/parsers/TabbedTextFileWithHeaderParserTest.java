package org.broadinstitute.hellbender.utils.text.parsers;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public final class TabbedTextFileWithHeaderParserTest {
    @Test
    public void basicParsingTest() throws Exception {
        final String[][] data = new String[][] {
                new String[] {"FOO", "BAR", "SPLAT"},
                new String[] {"1", "2", "3"},
                new String[] {"a", "b", "c"},
                new String[] {"foo", "bar", "splat"},
        };

        final File tmp = BaseTest.createTempFile("tabbedTextTest.", ".txt");
        try (final BufferedWriter out = new BufferedWriter(new FileWriter(tmp))) {

            for (final String[] fields : data) {
                out.write(StringUtil.join("\t", fields));
                out.newLine();
            }
        }

        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(tmp);
        for (final String col : data[0]) Assert.assertTrue(parser.hasColumn(col));

        int i=1;
        for (final TabbedTextFileWithHeaderParser.Row row : parser) {
            final String[] expected = data[i++];
            Assert.assertEquals(row.getFields(), expected);
            Assert.assertEquals(row.getCurrentLine(), StringUtil.join("\t", expected));
        }

        Assert.assertEquals(i, data.length);
    }

    @Test
    public void parsingWithColumnHeadersTest() throws Exception {
        final String[][] data = new String[][] {
                new String[] {"1", "2", "3"},
                new String[] {"a", "b", "2"},
                new String[] {"foo", "bar", ""},
        };

        final String[] headers = {"STRING", "STRING2", "NUMBER"};

        final File tmp = BaseTest.createTempFile("tabbedTextTest.", ".txt");

        try (final BufferedWriter out = new BufferedWriter(new FileWriter(tmp))) {

            for (final String[] fields : data) {
                out.write(StringUtil.join("\t", fields));
                out.newLine();
            }
        }

        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(tmp, headers);
        for (final String col : headers) Assert.assertTrue(parser.hasColumn(col));

        int i=0;
        for (final TabbedTextFileWithHeaderParser.Row row : parser) {
            final String[] expected = data[i++];
            final String[] actual = row.getFields();
            for (int j = 0; j < expected.length; j++) {
                Assert.assertTrue((expected[j].equals("") && actual[j] == null) || expected[j].equals(actual[j]));
            }
            Assert.assertEquals(row.getCurrentLine(), StringUtil.join("\t", expected));
            try {
                row.getField(headers[0]);
                row.getField(headers[1]);
                row.getIntegerField(headers[2]);
            }
            catch(Exception e) {
                Assert.fail("Failed to parse one of the fields in " + row.getCurrentLine() + ": " + e.getMessage());
            }
        }

        Assert.assertEquals(i, data.length);
    }

}
