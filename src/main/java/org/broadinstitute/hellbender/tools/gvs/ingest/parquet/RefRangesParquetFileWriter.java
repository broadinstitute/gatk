package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.ingest.RefRangesWriter;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Parquet writer for reference ranges data.
 * <p>
 * This writer extends AbstractParquetFileWriter for common Parquet infrastructure
 * and implements RefRangesWriter for domain-specific functionality.
 */
public class RefRangesParquetFileWriter extends AbstractParquetFileWriter implements RefRangesWriter {

    public RefRangesParquetFileWriter(Path file, MessageType schema, CompressionCodecName codecName) throws IOException {
        super(file, schema, codecName);
    }

    @Override
    public void commitData() {
        RefRangesWriter.super.commitData();
        try {
            this.close();
        } catch (IOException e) {
            throw new GATKException("Error while closing Parquet writer", e);
        }
    }

    @Override
    public void write(long location, long sampleId, int length, String state) throws IOException {
        JSONObject record = new JSONObject();
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("length", length);
        record.put("state", state);
        super.write(record);
    }

    @Override
    public void writeCompressed(long packedData, long sampleId) throws IOException {
        JSONObject compressedRecord = new JSONObject();
        compressedRecord.put("packed_ref_data", packedData);
        compressedRecord.put("sample_id", sampleId);
        super.write(compressedRecord);
    }
}
