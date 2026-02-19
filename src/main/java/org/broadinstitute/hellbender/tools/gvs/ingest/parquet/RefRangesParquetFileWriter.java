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
 * 
 * This writer extends RefRangesWriter for domain-specific functionality
 * and uses AbstractParquetFileWriter via composition for common Parquet writing infrastructure.
 */
public class RefRangesParquetFileWriter implements RefRangesWriter {
    private final AbstractParquetFileWriter delegate;

    public RefRangesParquetFileWriter(
            Path file,
            MessageType schema,
            CompressionCodecName codecName
    ) throws IOException {
        this.delegate = new AbstractParquetFileWriter(file, schema, codecName);
    }

    @Override
    public void close() throws IOException {
        this.delegate.close();
    }

    @Override
    public void commitData() {
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
        this.delegate.write(record);
    }

    @Override
    public void writeCompressed(long packedData, long sampleId) throws IOException {
        JSONObject compressedRecord = new JSONObject();
        compressedRecord.put("packedData", packedData);
        compressedRecord.put("sample_id", sampleId);
        this.delegate.write(compressedRecord);
    }
}
