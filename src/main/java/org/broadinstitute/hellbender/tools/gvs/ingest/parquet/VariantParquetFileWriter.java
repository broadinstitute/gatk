package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Parquet writer for variant data.
 * 
 * This writer stores variant information in Parquet format.
 */
public class VariantParquetFileWriter extends AbstractParquetFileWriter {

    public VariantParquetFileWriter(
            Path file,
            MessageType schema,
            CompressionCodecName codecName
    ) throws IOException {
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
}
