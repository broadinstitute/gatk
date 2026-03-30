package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.MessageTypeParser;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.IngestUtils;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class VrsAlleleCreatorUnitTest {

    private static final long   SAMPLE_ID  = 100L;
    private static final String SAMPLE_IDENTIFIER = "test";
    private static final String PROJECT_ID  = "proj";
    private static final String DATASET_NAME = "ds";

    private static final String ALLELE_ID   = "ga4gh:VA.testAllele";
    private static final String LOCATION_ID = "ga4gh:SL.testLoc";
    private static final String REFGET      = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";

    private static final int sampleTableNumber =
            IngestUtils.getTableNumber(SAMPLE_ID, IngestConstants.partitionPerTable);
    private static final String TABLE_NUMBER = String.format("%03d", sampleTableNumber);

    public static final MessageType VRS_PARQUET_SCHEMA = MessageTypeParser.parseMessageType("""
            message VrsAlleleRow {
                required binary vrs_allele_id (UTF8);
                required binary vrs_location_id (UTF8);
                required binary refget_accession (UTF8);
                required int64 ref_genome_coord_start;
                required int64 ref_genome_coord_end;
                required binary state_type (UTF8);
                optional binary state_sequence (UTF8);
                optional int64 state_length;
                optional int64 state_repeat_subunit_length;
            }
            """);

    private File outputDir;

    @BeforeMethod
    public void setUp() throws IOException {
        outputDir = Files.createTempDirectory("vrs_allele_creator_test").toFile();
    }

    @AfterMethod
    public void tearDown() {
        if (outputDir != null && outputDir.exists()) {
                        File[] files = outputDir.listFiles();
                        if (files != null) {
                                for (File f : files) {
                                        f.delete();
                                }
            }
            outputDir.delete();
        }
    }

        private VrsAlleleCreator newCreator(CommonCode.OutputType outputType) {
                return new VrsAlleleCreator(
                                SAMPLE_IDENTIFIER, SAMPLE_ID, TABLE_NUMBER, outputDir,
                                outputType, PROJECT_ID, DATASET_NAME, VRS_PARQUET_SCHEMA);
        }

        private File expectedTsvFile() {
                return new File(outputDir,
                                "vrs_allele_" + TABLE_NUMBER + "_" + SAMPLE_IDENTIFIER + IngestConstants.FILETYPE);
        }

        private List<String> readLines(File file) throws IOException {
                return Files.readAllLines(file.toPath());
        }

    // ------------------------------------------------------------------
    // Output file naming
    // ------------------------------------------------------------------

    @Test
    public void testTsvOutputFileNaming() {
        String filename = VrsAlleleCreator.getOutputFileName(
                TABLE_NUMBER, SAMPLE_ID, SAMPLE_IDENTIFIER, CommonCode.OutputType.TSV);
        Assert.assertEquals(filename, "vrs_allele_001_100_test.tsv");
    }

    @Test
    public void testParquetOutputFileNaming() {
        String filename = VrsAlleleCreator.getOutputFileName(
                TABLE_NUMBER, SAMPLE_ID, SAMPLE_IDENTIFIER, CommonCode.OutputType.PARQUET);
        Assert.assertEquals(filename, "vrs_allele_001_100_test.parquet");
    }

    @Test
    public void testOutputFileNamingDifferentSampleId() {
        String filename = VrsAlleleCreator.getOutputFileName(
                "002", 200L, "mysample", CommonCode.OutputType.PARQUET);
        Assert.assertEquals(filename, "vrs_allele_002_200_mysample.parquet");
    }

    // ------------------------------------------------------------------
    // apply() no-ops
    // ------------------------------------------------------------------

    @Test
    public void testApplyWithNullListIsNoop() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);
        creator.apply(null);  // must not throw
        creator.closeTool();
        // TSV file is created on construction; it should be empty (header only)
        File tsvFile = expectedTsvFile();
        Assert.assertTrue(tsvFile.exists());
        List<String> lines = readLines(tsvFile);
        Assert.assertEquals(lines.size(), 1, "Expected only header line when no records written");
    }

    @Test
    public void testApplyWithEmptyListIsNoop() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);
        creator.apply(Collections.emptyList());  // must not throw
        creator.closeTool();
        File tsvFile = expectedTsvFile();
        List<String> lines = readLines(tsvFile);
        Assert.assertEquals(lines.size(), 1, "Expected only header line when no records written");
    }

    // ------------------------------------------------------------------
    // TSV file creation and content
    // ------------------------------------------------------------------

    @Test
    public void testTsvCreatesFileWithHeader() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);
        creator.closeTool();

        File tsvFile = expectedTsvFile();
        Assert.assertTrue(tsvFile.exists(), "TSV file should be created on construction");
        String header = readLines(tsvFile).get(0);
        Assert.assertTrue(header.contains("vrs_allele_id"));
        Assert.assertTrue(header.contains("vrs_location_id"));
        Assert.assertTrue(header.contains("refget_accession"));
        Assert.assertTrue(header.contains("state_type"));
    }

    @Test
    public void testTsvLiteralRecord() throws IOException {
                VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);

        VrsAlleleRecord rec = new VrsAlleleRecord(
                ALLELE_ID, LOCATION_ID, REFGET, 100L, 101L, "T");
        creator.apply(List.of(rec));
        creator.closeTool();

        File tsvFile = expectedTsvFile();
        List<String> lines = readLines(tsvFile);
        Assert.assertEquals(lines.size(), 2, "Expected header + 1 data line");

        String dataLine = lines.get(1);
        Assert.assertTrue(dataLine.contains(ALLELE_ID));
        Assert.assertTrue(dataLine.contains(LOCATION_ID));
        Assert.assertTrue(dataLine.contains(REFGET));
        Assert.assertTrue(dataLine.contains("LiteralSequenceExpression"));
        Assert.assertTrue(dataLine.contains("T"));
    }

    @Test
    public void testTsvRleRecord() throws IOException {
                VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);

        VrsAlleleRecord rec = new VrsAlleleRecord(
                ALLELE_ID, LOCATION_ID, REFGET, 200L, 206L, 3, 3);
        creator.apply(List.of(rec));
        creator.closeTool();

        File tsvFile = expectedTsvFile();
        List<String> lines = readLines(tsvFile);
        Assert.assertEquals(lines.size(), 2, "Expected header + 1 data line");

        String dataLine = lines.get(1);
        Assert.assertTrue(dataLine.contains(ALLELE_ID));
        Assert.assertTrue(dataLine.contains("ReferenceLengthExpression"));
        Assert.assertTrue(dataLine.contains("3"));
    }

    @Test
    public void testTsvMultipleRecords() throws IOException {
                VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);

        List<VrsAlleleRecord> records = Arrays.asList(
                new VrsAlleleRecord("ga4gh:VA.id1", LOCATION_ID, REFGET, 100L, 101L, "A"),
                new VrsAlleleRecord("ga4gh:VA.id2", LOCATION_ID, REFGET, 100L, 101L, "T"),
                new VrsAlleleRecord("ga4gh:VA.id3", LOCATION_ID, REFGET, 100L, 101L, "G")
        );
        creator.apply(records);
        creator.closeTool();

        File tsvFile = expectedTsvFile();
        List<String> lines = readLines(tsvFile);
        Assert.assertEquals(lines.size(), 4, "Expected header + 3 data lines");
    }

    // ------------------------------------------------------------------
    // Parquet file creation
    // ------------------------------------------------------------------

    @Test
    public void testParquetCreatesFile() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.apply(List.of(
                new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 100L, 101L, "T")));
        creator.closeTool();

        String expectedFilename = VrsAlleleCreator.getOutputFileName(
                TABLE_NUMBER, SAMPLE_ID, SAMPLE_IDENTIFIER, CommonCode.OutputType.PARQUET);
        File parquetFile = new File(outputDir, expectedFilename);
        Assert.assertTrue(parquetFile.exists(), "Parquet file should exist after writing");
        Assert.assertTrue(parquetFile.length() > 0, "Parquet file should not be empty");
    }

    @Test
    public void testParquetFileAlreadyExistsThrows() throws IOException {
        String expectedFilename = VrsAlleleCreator.getOutputFileName(
                TABLE_NUMBER, SAMPLE_ID, SAMPLE_IDENTIFIER, CommonCode.OutputType.PARQUET);
        File parquetFile = new File(outputDir, expectedFilename);
        Files.createFile(parquetFile.toPath());

        try {
                        newCreator(CommonCode.OutputType.PARQUET);
            Assert.fail("Expected RuntimeException wrapping FileAlreadyExistsException");
        } catch (RuntimeException e) {
            Assert.assertNotNull(e.getCause());
        } finally {
            Files.deleteIfExists(parquetFile.toPath());
        }
    }

    // ------------------------------------------------------------------
    // commitData / closeTool are safe to call with no records written
    // ------------------------------------------------------------------

    @Test
    public void testCommitAndCloseWithNoWrites() throws IOException {
                VrsAlleleCreator creator = newCreator(CommonCode.OutputType.TSV);
        creator.commitData();
        creator.closeTool();  // must not throw
    }
}
