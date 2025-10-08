package org.broadinstitute.hellbender.utils.gvs.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.api.WriteSupport;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.io.IOException;

/**
 * A wrapper around ParquetWriter that simplifies writing JSON objects to Parquet files.
 * 
 * This implementation uses the modern composition pattern where WriteSupport is created
 * independently and passed to the Builder, rather than creating it inside a deprecated
 * getWriteSupport(Configuration) method. This approach provides better testability,
 * flexibility, and follows dependency injection principles.
 */
public class GvsVariantParquetFileWriter {
    private ParquetWriter<JSONObject> parquetWriterImpl;

    public GvsVariantParquetFileWriter(
            Path file,
            MessageType schema,
            CompressionCodecName codecName
    ) throws IOException {
        // Use composition pattern: create WriteSupport independently
        WriteSupport<JSONObject> writeSupport = new GvsVariantWriteSupport(schema);
        
        // Use our modern Builder that accepts WriteSupport via constructor
        Builder builder = new Builder(file, writeSupport);
        this.parquetWriterImpl = builder.withCompressionCodec(codecName).build();
    }

    public void write(JSONObject object) throws IOException {
        this.parquetWriterImpl.write(object);
    }

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
     * but there's no alternative method available yet. The @SuppressWarnings annotation
     * is used because:
     * 1. This is the only way to implement ParquetWriter.Builder in this version
     * 2. We're using best practices (composition) despite the deprecated signature
     * 3. The Configuration parameter is ignored as intended by the new design
     */
    public static class Builder extends ParquetWriter.Builder<JSONObject, Builder> {
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
         * We suppress the deprecation warning because:
         * - This is the only way to extend ParquetWriter.Builder in Parquet 1.13.1
         * - We follow modern composition patterns by injecting WriteSupport via constructor
         * - The Configuration parameter is properly ignored (not used)
         * - This approach is more testable and flexible than creating WriteSupport here
         * 
         * @param configuration ignored parameter (as intended by modern design)
         * @return the WriteSupport instance provided at construction time
         */
        @Override
        @SuppressWarnings("deprecation")
        protected WriteSupport<JSONObject> getWriteSupport(Configuration configuration) {
            return this.writeSupport;
        }
    }
}
