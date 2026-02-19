package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.ingest.SamplePloidyWriter;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Parquet writer for sample ploidy data.
 * 
 * This writer extends SamplePloidyWriter for domain-specific functionality
 * and uses AbstractParquetFileWriter via composition for common Parquet writing infrastructure.
 */
public class SamplePloidyParquetFileWriter implements SamplePloidyWriter {
    private final AbstractParquetFileWriter delegate;

    public SamplePloidyParquetFileWriter(
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
    public void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("chromosome", chromosome);
        record.put("ploidy", ploidy);
        this.delegate.write(record);
    }

    @Override
    public void commitData() {
        SamplePloidyWriter.super.commitData();
        try {
            this.close();
        } catch (IOException e) {
            throw new GATKException("Error while closing Parquet writer", e);
        }
    }
}
