package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Parquet writer for VCF header information.
 * <p>
 * This writer stores sample ID and header line hash associations in Parquet format.
 */
public class HeaderParquetFileWriter extends AbstractParquetFileWriter {

    public HeaderParquetFileWriter(Path file, MessageType schema, CompressionCodecName codecName) throws IOException {
        super(file, schema, codecName);
    }

    /**
     * Writes a JSON object to the Parquet file.
     * 
     * @param object the JSON object to write
     * @throws IOException if an error occurs during writing
     */
    public void write(JSONObject object) throws IOException {
        super.write(object);
    }

    /**
     * Creates a JSON object representing a header record.
     * <p>
     * Column names intentionally match the {@code vcf_header_lines_scratch} BigQuery table so the
     * resulting Parquet file can be loaded into that table by name (see the VS-1968 design doc).
     *
     * @param sampleId the sample ID
     * @param headerChunk the header text chunk; may be {@code null} to write an association-only row
     *                    (sample_id -> hash) without the text, mirroring the BQ path's dedup behavior
     * @param headerLineHash the MD5 hash of the header chunk
     * @param isExpectedUnique whether this chunk is expected to differ per sample
     * @return a JSON object containing the header information
     */
    public static JSONObject writeJson(Long sampleId, String headerChunk, String headerLineHash, Boolean isExpectedUnique) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        // Omit the (nullable) text column entirely when null so the Parquet writer skips it.
        if (headerChunk != null) {
            record.put("vcf_header_lines", headerChunk);
        }
        record.put("vcf_header_lines_hash", headerLineHash);
        record.put("is_expected_unique", isExpectedUnique);
        return record;
    }
}
