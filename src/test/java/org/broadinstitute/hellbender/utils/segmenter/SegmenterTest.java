package org.broadinstitute.hellbender.utils.segmenter;

import org.testng.annotations.Test;
import org.testng.Assert;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import org.apache.commons.csv.*;

public class SegmenterTest {
    public static File createTempFile(String prefix, String suffix) {
        try {
            final File result = File.createTempFile(prefix, suffix, null);
            result.deleteOnExit();
            return result;
        } catch (IOException e) {
            Assert.fail("Failure on creating temp file.");
            return null;
        }
    }

    @Test
    public void testHCC1143Reduced() throws IOException {
        float min_log_value = -10;
        final File INPUT_FILE = new File("src/test/resources/segmenter/input/HCC1143_reduced.tsv");
        final File EXPECTED = new File("src/test/resources/segmenter/output/HCC1143_reduced_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        String sampleName = "HCC1143";
        File segmentFile = RCBSSegmenter.segment(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), min_log_value);
        Assert.assertTrue(segmentFile.exists(), "R library was not written to temp file: " + output);

        final CSVParser parser = new CSVParser(new FileReader(output), CSVFormat.TDF.withHeader());
        final CSVParser expectedParser = new CSVParser(new FileReader(EXPECTED), CSVFormat.TDF.withHeader());
        final List<CSVRecord> resultRecord = parser.getRecords();
        final List<CSVRecord> expectedRecord = expectedParser.getRecords();
        Assert.assertEquals(resultRecord.size(), expectedRecord.size());
        for (int i=0; i<resultRecord.size(); i++) {
            Assert.assertEquals(Double.parseDouble(resultRecord.get(i).get("Segment_Mean")),
                    Double.parseDouble(expectedRecord.get(i).get("Segment_Mean")), 0.0000000000001,
                    "Different on line: "+(i+2));
        }
        parser.close();
        expectedParser.close();
    }

    @Test
    public void testHCC1143Short() throws IOException {
        float min_log_value = -10;
        final File INPUT_FILE = new File("src/test/resources/segmenter/input/HCC1143_short.tsv");
        final File EXPECTED = new File("src/test/resources/segmenter/output/HCC1143_short_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        String sampleName = "HCC1143";
        File segmentFile = RCBSSegmenter.segment(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), min_log_value);
        Assert.assertTrue(segmentFile.exists(), "R library was not written to temp file: " + output);

        final CSVParser parser = new CSVParser(new FileReader(output), CSVFormat.TDF.withHeader());
        final CSVParser expectedParser = new CSVParser(new FileReader(EXPECTED), CSVFormat.TDF.withHeader());
        final List<CSVRecord> resultRecord = parser.getRecords();
        final List<CSVRecord> expectedRecord = expectedParser.getRecords();
        Assert.assertEquals(resultRecord.size(), expectedRecord.size());
        for (int i=0; i<resultRecord.size(); i++) {
            Assert.assertEquals(Double.parseDouble(resultRecord.get(i).get("Segment_Mean")),
                    Double.parseDouble(expectedRecord.get(i).get("Segment_Mean")), 0.0000000000001,
                    "Different on line: "+(i+2));
        }
        parser.close();
        expectedParser.close();
    }

    @Test
    public void testSimple() throws IOException {
        float min_log_value = -10;
        final File INPUT_FILE = new File("src/test/resources/segmenter/input/Simple.tsv");
        final File EXPECTED = new File("src/test/resources/segmenter/output/Simple_result.seg");
        final File output = createTempFile("recapseg.HCC1143", ".seg");
        String sampleName = "Simple";
        File segmentFile = RCBSSegmenter.segment(sampleName, INPUT_FILE.getAbsolutePath(), output.getAbsolutePath(), min_log_value);
        Assert.assertTrue(segmentFile.exists(), "R library was not written to temp file: " + output);

        final CSVParser parser = new CSVParser(new FileReader(output), CSVFormat.TDF.withHeader());
        final CSVParser expectedParser = new CSVParser(new FileReader(EXPECTED), CSVFormat.TDF.withHeader());
        final List<CSVRecord> resultRecord = parser.getRecords();
        final List<CSVRecord> expectedRecord = expectedParser.getRecords();
        Assert.assertEquals(resultRecord.size(), expectedRecord.size());
        for (int i=0; i<resultRecord.size(); i++) {
            Assert.assertEquals(Double.parseDouble(resultRecord.get(i).get("Segment_Mean")),
                    Double.parseDouble(expectedRecord.get(i).get("Segment_Mean")), 0.0000000000001,
                    "Different on line: "+(i+2));
        }
        parser.close();
        expectedParser.close();
    }
}