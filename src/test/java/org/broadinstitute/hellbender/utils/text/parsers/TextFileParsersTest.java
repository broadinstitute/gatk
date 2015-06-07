package org.broadinstitute.hellbender.utils.text.parsers;

import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

public final class TextFileParsersTest {
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/utils/text/parsers");
    private static final File testFile1 = new File(TEST_DIR, "whitespace_text_file.txt");
    private static final File testFile2 = new File(TEST_DIR, "all_ones_text_file.txt");
    private static final File testFile3 = new File(TEST_DIR, "no_grouping_file.txt");
    private static final File testFile4 = new File(TEST_DIR, "tabbed_text_file.txt");
    // There is a comment in the file data that should be skipped by the parser, so it is not included below
    private static final Object[][] testFile1Data = {
        { "Now", "is", "the", "time" },
        { "for", "all", "good", "men" },
        { "to", "come", "to", "the" },
        { "aid", "of", "their", "country." },
        { 15.0d, 23, 55, 67.88888}
    };

    @Test(dataProvider = "basicInputParserData")
    public void testTextFileParser(Object fileOrStream) throws IOException {
        FormatUtil format = new FormatUtil();

        List<String> expected = new ArrayList<>();
        if (fileOrStream instanceof File) {
            try (BufferedReader reader = new BufferedReader(new FileReader((File)fileOrStream))) {
                String line = null;
                while ((line = reader.readLine()) != null) {
                    if (!line.startsWith("#")) expected.add(line);
                }
            }
        }

        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, (File)fileOrStream )
            : new BasicInputParser(true, (InputStream)fileOrStream);
        int index = 0;
        while (parser.hasNext())
        {
            String parts[] = parser.next();
            if (fileOrStream instanceof File) {
                // Can't have the parser and the reader workking with an InputStream at the same time
                // so we only do this test with the file
                Assert.assertEquals(parser.getCurrentLine(), expected.get(index));
            }
            // Line 4 is a comment, so there's a gap in the line numbers
            Assert.assertEquals(parser.getCurrentLineNumber(), index <= 2 ? index + 1 : index + 2);
            Assert.assertEquals(parts.length, 4);
            if (index < 4) {
                for (int i = 0; i < parts.length; i++) {
                    Assert.assertEquals(parts, testFile1Data[index]);
                }
            }
            else {
                Assert.assertEquals(testFile1Data[index][0], format.parseDouble(parts[0]));
                Assert.assertEquals(testFile1Data[index][1], format.parseInt(parts[1]));
                Assert.assertEquals(testFile1Data[index][2], format.parseInt(parts[2]));
                Assert.assertEquals(testFile1Data[index][3], format.parseDouble(parts[3]));
            }
            index++;
        }
    }

    @DataProvider(name = "basicInputParserData")
    private Object[][] getBasicInputParserData()
    {
        return new Object[][] {
                {testFile1},
                {IOUtil.openFileForReading(testFile1)}
        };
    }

    @Test(dataProvider = "multiFileParsingData")
    public void testMultiFileParsing(Object fileOrStream1, Object fileOrStream2) throws IOException {
        FormatUtil format = new FormatUtil();

        List<String> expected = new ArrayList<>();
        if (fileOrStream1 instanceof File) {
            try (final BufferedReader reader = new BufferedReader(new FileReader((File)fileOrStream1))) {
                String line = null;
                while ((line = reader.readLine()) != null) {
                    if (!line.startsWith("#")) expected.add(line);
                }
            }
            try (final BufferedReader reader = new BufferedReader(new FileReader((File)fileOrStream2))){
                String line = null;
                while ((line = reader.readLine()) != null) {
                    if (!line.startsWith("#")) expected.add(line);
                }
            }
        }

        BasicInputParser parser = fileOrStream1 instanceof File
            ? new BasicInputParser(true, (File)fileOrStream1, (File)fileOrStream2 )
            : new BasicInputParser(true, (InputStream)fileOrStream1, (InputStream)fileOrStream2);
        int index = 0;
        // Line 4 is a comment, so there's a gap in the line numbers
        int expectedLineNumbers[] = {1,2,3,5,6,1,2,3,5,6};
        while (parser.hasNext())
        {
            String parts[] = parser.next();
            if (fileOrStream1 instanceof File) {
                // Can't have the parser and the reader working with an InputStream at the same time
                // so we only test the files
                Assert.assertEquals(parser.getCurrentLine(), expected.get(index));
            }
            Assert.assertEquals(parser.getCurrentLineNumber(), expectedLineNumbers[index]);
            Assert.assertEquals(parts.length, 4);
            int indexIntoTestData = (index<5) ? index : index-5;
            if (index != 4 && index != 9) {
                for (int i = 0; i < parts.length; i++) {
                    Assert.assertEquals(parts, testFile1Data[indexIntoTestData]);
                }
            }
            else {
                Assert.assertEquals(testFile1Data[indexIntoTestData][0], format.parseDouble(parts[0]));
                Assert.assertEquals(testFile1Data[indexIntoTestData][1], format.parseInt(parts[1]));
                Assert.assertEquals(testFile1Data[indexIntoTestData][2], format.parseInt(parts[2]));
                Assert.assertEquals(testFile1Data[indexIntoTestData][3], format.parseDouble(parts[3]));
            }
            index++;
        }
    }

    @DataProvider(name = "multiFileParsingData")
    private Object[][] getMultiFileParsingData()
    {
        return new Object[][] {
                {testFile1, testFile1},
                {IOUtil.openFileForReading(testFile1), IOUtil.openFileForReading(testFile1)}
        };
    }

    @Test(dataProvider = "noGroupingData")
    public void testTextFileParserNoGrouping(Object fileOrStream) {
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, (File)fileOrStream)
            : new BasicInputParser(true, (InputStream)fileOrStream);
        parser.setTreatGroupedDelimitersAsOne(false);
        while (parser.hasNext()) {
            String parts[] = parser.next();
            for (int i = 0; i < parts.length; i++) {
                if (parts[i] != null) {
                    Assert.assertEquals(Integer.parseInt(parts[i]), i + 1);
                }
            }
        }
    }

    @DataProvider(name = "noGroupingData")
    private Object[][] getNoGroupingData()
    {
        return new Object[][] {
                {testFile3},
                {IOUtil.openFileForReading(testFile3)}
        };
    }


    @Test(dataProvider = "leadingWhiteSpaceData")
    public void testTextFileParserLeadingWhitespace(Object fileOrStream) {
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, (File)fileOrStream)
            : new BasicInputParser(true, (InputStream)fileOrStream);
        while (parser.hasNext())
        {
            String parts[] = parser.next();
            Assert.assertEquals(parts.length, 1);
            Assert.assertEquals("1", parts[0]);
        }
    }

    @DataProvider(name = "leadingWhiteSpaceData")
    private Object[][] getLeadingWhiteSpaceData()
    {
        return new Object[][] {
                {testFile2},
                {IOUtil.openFileForReading(testFile2)}
        };
    }


    @Test(expectedExceptions= RuntimeIOException.class, dataProvider = "tooManyWordsData")
    public void testTooManyWords(Object fileOrStream) {
        BasicInputParser parser = fileOrStream instanceof File
            ? new BasicInputParser(true, 3, (File)fileOrStream)
            : new BasicInputParser(true, 3, (InputStream)fileOrStream);
        if (parser.hasNext()) {
            String parts[] = parser.next();
        }
        Assert.fail("Attempt to parse extra-long file should have failed but didn't.");
    }

    @DataProvider(name = "tooManyWordsData")
    private Object[][] getTooManyWordsData()
    {
        return new Object[][] {
                {testFile1},
                {IOUtil.openFileForReading(testFile1)}
        };
    }

    @Test(dataProvider = "tabbedData")
    public void testTabbedFileParser(Object fileOrStream) {
        TabbedInputParser parser = fileOrStream instanceof File
            ? new TabbedInputParser(false, (File)fileOrStream)
            : new TabbedInputParser(false, (InputStream)fileOrStream);
        while (parser.hasNext()) {
            String parts[] = parser.next();
            for (int i = 0; i < parts.length; i++) {
                if (parts[i] != null && !parts[i].equals("")) {
                    Assert.assertEquals(parts[i].trim(), String.valueOf(i + 1));
                }
            }
        }
    }

    @DataProvider(name = "tabbedData")
    private Object[][] getTabbedData()
    {
        return new Object[][] {
                {testFile4},
                {IOUtil.openFileForReading(testFile4)}
        };
    }

    @Test(dataProvider="data")
    public void testWordCountCalculation(String line, boolean groupDelimiters, String name) {

        WordCountTestParser parser = new WordCountTestParser();
        parser.setDelimiter("\t ");
        parser.setTreatGroupedDelimitersAsOne(groupDelimiters);
        parser.calculateWordCount(line.getBytes());
        Assert.assertEquals(parser.getWordCount(), 3, name);
    }

    @DataProvider(name = "data")
    private Object[][] getWordCountCalculationData()
    {
        return new Object[][]{
                {"1\t2\t3", false, "Tabs with all fields filled."},
                {"1\t2\t", false, "Tabs with no final field."},
                {"\t2\t3", false, "Tabs with no first field."},
                {"\t2\t", false, "Tabs with no first or final field."},
                {"1  2  3", true, "Spaces with all fields filled  (grouping on)."},
                {"1  2  3  ", true, "Spaces with no final field (grouping on)."},
                {"   2   3   4", true, "Spaces with no first field (grouping on)."},
                {" 2 ", false, "Spaces with no first or final field."}
        };
    }

    /**
     * Toy class for testing the word count functionality
     */
    private static class WordCountTestParser extends AbstractInputParser {

        private char delimiters[] = null;
        
        public WordCountTestParser() {
        }

        public void setDelimiter(String delim) {
            delimiters = delim.toCharArray();
        }

        protected boolean isDelimiter(final byte b) {
            for (int i = 0; i < delimiters.length; i++) {
                if (b == delimiters[i]) {
                    return true;
                }
            }
            return false;
        }

        protected byte[] readNextLine() {  return new byte[0]; }
        public String getFileName() { return null; }
        public void close() {}
    }
}
