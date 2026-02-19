package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.api.WriteSupport;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Delegate helper class for Parquet file writers that simplifies writing JSON objects to Parquet files.
 * 
 * This implementation uses the modern composition pattern where WriteSupport is created
 * independently and passed to the Builder, rather than creating it inside a deprecated
 * getWriteSupport(Configuration) method. This approach provides better testability,
 * flexibility, and follows dependency injection principles.
 * 
 * This class is designed to be used via composition by concrete writer implementations,
 * allowing them to extend domain-specific abstract classes while reusing common Parquet logic.
 */
public class AbstractParquetFileWriter {
    private final ParquetWriter<JSONObject> parquetWriterImpl;

    /**
     * Creates a Parquet file writer with the specified configuration.
     * 
     * @param file the output file path
     * @param schema the Parquet message type schema
     * @param codecName the compression codec to use
     * @throws IOException if an error occurs during writer initialization
     */
    public AbstractParquetFileWriter(
            Path file,
            MessageType schema,
            CompressionCodecName codecName
    ) throws IOException {
        // Use composition pattern: create WriteSupport independently
        WriteSupport<JSONObject> writeSupport = new ParquetWriteSupport(schema);
        
        // Use our modern Builder that accepts WriteSupport via constructor
        Builder builder = new Builder(file, writeSupport);
        this.parquetWriterImpl = builder.withCompressionCodec(codecName).build();
    }

    /**
     * Writes a JSON object to the Parquet file.
     * 
     * @param object the JSON object to write
     * @throws IOException if an error occurs during writing
     */
    public void write(JSONObject object) throws IOException {
        this.parquetWriterImpl.write(object);
    }

    /**
     * Closes the Parquet writer and releases resources.
     * 
     * @throws IOException if an error occurs during closing
     */
    public void close() throws IOException {
        this.parquetWriterImpl.close();
    }

    /**
     * Modern Builder implementation using composition pattern.
     * 
     * The WriteSupport is injected via constructor rather than created inside
     * getWriteSupport(), which provides better separation of concerns and testability.
     * 
     * Note: The getWriteSupport(Configuration) method is deprecated in Parquet 1.13.1
     * but there's no alternative method available yet.
     */
    private static class Builder extends ParquetWriter.Builder<JSONObject, Builder> {
        private final WriteSupport<JSONObject> writeSupport;

        /**
         * Creates a builder with injected WriteSupport (composition pattern).
         * 
         * @param file the output file path
         * @param writeSupport the WriteSupport instance to use (injected dependency)
         */
        public Builder(Path file, WriteSupport<JSONObject> writeSupport) {
            super(file);
            this.writeSupport = writeSupport;
        }

        @Override
        protected Builder self() {
            return this;
        }

        /**
         * Returns the WriteSupport instance provided via constructor.
         * 
         * This method is deprecated but still required to implement ParquetWriter.Builder.
         * The Configuration parameter is properly ignored (not used).
         * 
         * @param configuration ignored parameter (as intended by modern design)
         * @return the WriteSupport instance provided at construction time
         */
        @Deprecated
        @Override
        protected WriteSupport<JSONObject> getWriteSupport(Configuration configuration) {
            return this.writeSupport;
        }
    }
}
