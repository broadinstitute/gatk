package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.example.data.Group;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.hadoop.example.GroupReadSupport;
import org.apache.parquet.schema.GroupType;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.MessageTypeParser;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.HeaderParquetFileWriter;
import org.broadinstitute.hellbender.utils.Utils;
import org.json.JSONObject;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Unit tests for VS-1803 header-loading-on-Parquet: the {@code writeJson} record shape and the
 * end-to-end NAIVE Parquet write path (schema + {@code is_expected_unique} BOOLEAN support).
 */
public class VcfHeaderLineScratchCreatorUnitTest {

    private static final long SAMPLE_ID = 100;

    // Must match CreateVariantIngestFiles#headersRowSchema and the vcf_header_lines_scratch BQ table.
    private static final MessageType HEADERS_ROW_SCHEMA = MessageTypeParser.parseMessageType("""
            message HeaderRow {
                required int64 sample_id;
                optional binary vcf_header_lines (UTF8);
                required binary vcf_header_lines_hash (UTF8);
                required boolean is_expected_unique;
            }
            """);

    @Test
    public void testWriteJsonFullRecord() {
        final JSONObject record = HeaderParquetFileWriter.writeJson(5L, "##contig=<ID=chr1>", "hash1", true);
        Assert.assertEquals(record.getLong("sample_id"), 5L);
        Assert.assertEquals(record.getString("vcf_header_lines"), "##contig=<ID=chr1>");
        Assert.assertEquals(record.getString("vcf_header_lines_hash"), "hash1");
        Assert.assertTrue(record.getBoolean("is_expected_unique"));
    }

    @Test
    public void testWriteJsonNullChunkOmitsText() {
        // A null chunk yields an association-only row: the (nullable) text column must be omitted so
        // the Parquet writer skips it, matching the BQ path's NULL-text dedup rows.
        final JSONObject record = HeaderParquetFileWriter.writeJson(5L, null, "hash1", false);
        Assert.assertFalse(record.has("vcf_header_lines"), "text column should be absent for a null chunk");
        Assert.assertEquals(record.getLong("sample_id"), 5L);
        Assert.assertEquals(record.getString("vcf_header_lines_hash"), "hash1");
        Assert.assertFalse(record.getBoolean("is_expected_unique"));
    }

    /**
     * End-to-end NAIVE strategy: write header chunks via VcfHeaderLineScratchCreator to a real Parquet
     * file, read it back, and assert every column round-trips — including the boolean is_expected_unique,
     * which exercises the BOOLEAN case added to ParquetWriteSupport.
     */
    @Test
    public void testNaiveParquetRoundTrip() throws IOException {
        final File outputDir = Files.createTempDirectory("vcf-header-parquet-test").toFile();
        try {
            // A shared blob (is_expected_unique=false) plus a per-sample command line (true).
            final Map<String, Boolean> headers = new LinkedHashMap<>();
            headers.put("##contig=<ID=chr1>\n##contig=<ID=chr2>\n##INFO=<ID=AF>", false);
            headers.put("##GATKCommandLine=<ID=HaplotypeCaller,sampleA>", true);

            final VcfHeaderLineScratchCreator creator = new VcfHeaderLineScratchCreator(
                    SAMPLE_ID, "test", "test", outputDir, CommonCode.OutputType.PARQUET, HEADERS_ROW_SCHEMA);
            creator.apply(headers);
            creator.commitData();

            final File parquetFile = new File(outputDir, "vcf_header_lines_scratch_" + SAMPLE_ID + ".parquet");
            Assert.assertTrue(parquetFile.exists(), "expected vcf_header_lines_scratch_<sampleId>.parquet to be written");

            // hash -> (text, is_expected_unique) as read back from Parquet.
            final Map<String, String> readText = new HashMap<>();
            final Map<String, Boolean> readUnique = new HashMap<>();
            int rowCount = 0;
            try (ParquetReader<Group> reader = ParquetReader
                    .builder(new GroupReadSupport(), new Path(parquetFile.toURI()))
                    .withConf(new Configuration())
                    .build()) {
                Group row;
                while ((row = reader.read()) != null) {
                    rowCount++;
                    final GroupType schema = row.getType();
                    Assert.assertEquals(row.getLong("sample_id", 0), SAMPLE_ID);
                    final String hash = row.getString("vcf_header_lines_hash", 0);
                    // NAIVE writes the full text for every chunk.
                    Assert.assertEquals(row.getFieldRepetitionCount(schema.getFieldIndex("vcf_header_lines")), 1,
                            "NAIVE strategy should always write the text column");
                    final String text = row.getString("vcf_header_lines", 0);
                    final boolean unique = row.getBoolean("is_expected_unique", 0);
                    readText.put(hash, text);
                    readUnique.put(hash, unique);
                }
            }

            Assert.assertEquals(rowCount, headers.size(), "one Parquet row per header chunk");
            for (final Map.Entry<String, Boolean> chunk : headers.entrySet()) {
                final String expectedHash = Utils.calcMD5(chunk.getKey());
                Assert.assertEquals(readText.get(expectedHash), chunk.getKey(), "text should round-trip by hash");
                Assert.assertEquals(readUnique.get(expectedHash), chunk.getValue(), "is_expected_unique should round-trip");
            }
        } finally {
            deleteRecursively(outputDir);
        }
    }

    private static void deleteRecursively(final File file) {
        final File[] children = file.listFiles();
        if (children != null) {
            for (final File child : children) {
                deleteRecursively(child);
            }
        }
        // best-effort cleanup
        //noinspection ResultOfMethodCallIgnored
        file.delete();
    }
}
