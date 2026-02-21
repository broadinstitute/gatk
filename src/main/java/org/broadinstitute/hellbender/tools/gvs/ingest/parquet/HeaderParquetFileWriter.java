package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Parquet writer for VCF header information.
 * 
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
     * 
     * @param sampleId the sample ID
     * @param headerLineHash the hash of the header line
     * @return a JSON object containing the header information
     */
    public static JSONObject writeJson(Long sampleId, String headerLineHash) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("headerLineHash", headerLineHash);
        return record;
    }
}
