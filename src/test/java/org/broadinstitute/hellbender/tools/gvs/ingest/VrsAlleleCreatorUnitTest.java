package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.example.data.Group;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.hadoop.example.GroupReadSupport;
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
import java.util.ArrayList;
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

        private File expectedParquetFile() {
                return new File(outputDir, VrsAlleleCreator.getOutputFileName(
                                TABLE_NUMBER, SAMPLE_ID, SAMPLE_IDENTIFIER, CommonCode.OutputType.PARQUET));
        }

        private List<Group> readRows(File file) throws IOException {
                final List<Group> rows = new ArrayList<>();
                try (ParquetReader<Group> reader = ParquetReader
                                .builder(new GroupReadSupport(), new Path(file.toURI()))
                                .withConf(new Configuration())
                                .build()) {
                        Group row;
                        while ((row = reader.read()) != null) {
                                rows.add(row);
                        }
                }
                return rows;
        }

        private static boolean hasValue(Group row, String field) {
                return row.getFieldRepetitionCount(row.getType().getFieldIndex(field)) > 0;
        }

    // ------------------------------------------------------------------
    // Output file naming
    // ------------------------------------------------------------------

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
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.apply(null);  // must not throw
        creator.closeTool();
        // The Parquet file is created on construction; it should hold no rows
        Assert.assertEquals(readRows(expectedParquetFile()).size(), 0, "Expected no rows when no records written");
    }

    @Test
    public void testApplyWithEmptyListIsNoop() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.apply(Collections.emptyList());  // must not throw
        creator.closeTool();
        Assert.assertEquals(readRows(expectedParquetFile()).size(), 0, "Expected no rows when no records written");
    }

    // ------------------------------------------------------------------
    // Parquet content
    // ------------------------------------------------------------------

    @Test
    public void testParquetLiteralRecord() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.apply(List.of(new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 100L, 101L, "T")));
        creator.closeTool();

        List<Group> rows = readRows(expectedParquetFile());
        Assert.assertEquals(rows.size(), 1);
        Group row = rows.get(0);
        Assert.assertEquals(row.getString("vrs_allele_id", 0), ALLELE_ID);
        Assert.assertEquals(row.getString("vrs_location_id", 0), LOCATION_ID);
        Assert.assertEquals(row.getString("refget_accession", 0), REFGET);
        Assert.assertEquals(row.getLong("ref_genome_coord_start", 0), 100L);
        Assert.assertEquals(row.getLong("ref_genome_coord_end", 0), 101L);
        Assert.assertEquals(row.getString("state_type", 0), "LiteralSequenceExpression");
        Assert.assertEquals(row.getString("state_sequence", 0), "T");
        Assert.assertFalse(hasValue(row, "state_length"));
        Assert.assertFalse(hasValue(row, "state_repeat_subunit_length"));
    }

    @Test
    public void testParquetReferenceLengthRecord() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.apply(List.of(new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 200L, 206L, 3, 2)));
        creator.closeTool();

        List<Group> rows = readRows(expectedParquetFile());
        Assert.assertEquals(rows.size(), 1);
        Group row = rows.get(0);
        Assert.assertEquals(row.getString("vrs_allele_id", 0), ALLELE_ID);
        Assert.assertEquals(row.getString("state_type", 0), "ReferenceLengthExpression");
        Assert.assertFalse(hasValue(row, "state_sequence"));
        Assert.assertEquals(row.getLong("state_length", 0), 3L);
        Assert.assertEquals(row.getLong("state_repeat_subunit_length", 0), 2L);
    }

    @Test
    public void testParquetMultipleRecords() throws IOException {
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.apply(Arrays.asList(
                new VrsAlleleRecord("ga4gh:VA.id1", LOCATION_ID, REFGET, 100L, 101L, "A"),
                new VrsAlleleRecord("ga4gh:VA.id2", LOCATION_ID, REFGET, 100L, 101L, "T"),
                new VrsAlleleRecord("ga4gh:VA.id3", LOCATION_ID, REFGET, 100L, 101L, "G")
        ));
        creator.closeTool();

        List<Group> rows = readRows(expectedParquetFile());
        Assert.assertEquals(rows.size(), 3);
        Assert.assertEquals(rows.get(0).getString("vrs_allele_id", 0), "ga4gh:VA.id1");
        Assert.assertEquals(rows.get(2).getString("vrs_allele_id", 0), "ga4gh:VA.id3");
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
        VrsAlleleCreator creator = newCreator(CommonCode.OutputType.PARQUET);
        creator.commitData();
        creator.closeTool();  // must not throw
    }
}
