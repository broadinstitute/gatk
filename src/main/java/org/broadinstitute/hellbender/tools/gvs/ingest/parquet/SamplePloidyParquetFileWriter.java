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
 * This writer extends AbstractParquetFileWriter for common Parquet infrastructure
 * and implements SamplePloidyWriter for domain-specific functionality.
 */
public class SamplePloidyParquetFileWriter extends AbstractParquetFileWriter implements SamplePloidyWriter {

    public SamplePloidyParquetFileWriter(
            Path file,
            MessageType schema,
            CompressionCodecName codecName
    ) throws IOException {
        super(file, schema, codecName);
    }

    @Override
    public void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("chromosome", chromosome);
        record.put("ploidy", ploidy);
        super.write(record);
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
