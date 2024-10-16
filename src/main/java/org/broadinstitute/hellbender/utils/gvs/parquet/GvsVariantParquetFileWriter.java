package org.broadinstitute.hellbender.utils.gvs.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.column.ParquetProperties;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.api.WriteSupport;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.io.OutputFile;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;


import java.io.IOException;

public class GvsVariantParquetFileWriter {
    /**
     * Let the top-level class simply create and delegate to the actual writer, so users of it
     * don't have to worry about details like builders and obscure arguments
     */
    private ParquetWriter<JSONObject> parquetWriterImpl;

    public GvsVariantParquetFileWriter(
            Path file,
            MessageType schema,
            CompressionCodecName codecName
    ) throws IOException {
        GvsVariantParquetFileWriter.Builder builder = new Builder(file);
        this.parquetWriterImpl = builder.withType(schema)
                .withCompressionCodec(codecName).build();
    }

    public void write(JSONObject object) throws IOException {
        this.parquetWriterImpl.write(object);
    }

    public void close() throws IOException {
        this.parquetWriterImpl.close();
    }

    public static class Builder extends ParquetWriter.Builder<JSONObject, Builder> {
        private MessageType schema = null;

        private Builder(Path file) {
            super(file);
        }

        public Builder withType(MessageType type) {
            this.schema = type;
            return this;
        }
        @Override
        protected Builder self() {
            return this;
        }

        @Override
        protected WriteSupport<JSONObject> getWriteSupport(Configuration configuration) {
            return new GvsVariantWriteSupport(schema);
        }
    }
}
