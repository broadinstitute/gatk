package org.broadinstitute.hellbender.utils.gvs.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileAlreadyExistsException;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.column.ParquetProperties;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.io.OutputFile;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.io.IOException;

public class GvsReferenceParquetFileWriter extends ParquetWriter<JSONObject> {

    /**
     * This is very deprecated, and we'll need to figure out how to do this from a builder once it works!
     * @param file
     * @param schema
     * @param enableDictionary
     * @param codecName
     * @throws IOException
     */
    public GvsReferenceParquetFileWriter(
            Path file,
            MessageType schema,
            boolean enableDictionary,
            CompressionCodecName codecName
    ) throws FileAlreadyExistsException, IOException {
        super(file, new GvsReferenceWriteSupport(schema), codecName, DEFAULT_BLOCK_SIZE, DEFAULT_PAGE_SIZE, enableDictionary, false);
    }

    GvsReferenceParquetFileWriter(
            Path file,
            GvsVariantWriteSupport writeSupport,
            CompressionCodecName compressionCodecName,
            int blockSize,
            int pageSize,
            boolean enableDictionary,
            boolean enableValidation,
            ParquetProperties.WriterVersion writerVersion,
            Configuration conf)
            throws IOException {
        super(
                file,
                writeSupport,
                compressionCodecName,
                blockSize,
                pageSize,
                pageSize,
                enableDictionary,
                enableValidation,
                writerVersion,
                conf);
    }

    public static JSONObject writeJson(long location, Long sampleId, int length, String state) {
        JSONObject record = new JSONObject();
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("length", length);
        record.put("state", state);
        return record;
    }

    public static JSONObject writeCompressed(long packedData, long sampleId) {
        JSONObject compressedRecord = new JSONObject();
        compressedRecord.put("packedData", packedData);
        compressedRecord.put("sample_id", sampleId);
        return compressedRecord;
    }


    public static class Builder extends ParquetWriter.Builder<JSONObject, Builder> {
        private MessageType schema = null;

        private Builder(Path file) {
            super(file);
        }

        private Builder(OutputFile file) {
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
        protected GvsReferenceWriteSupport getWriteSupport(Configuration conf) {
            return new GvsReferenceWriteSupport(schema);
        }
    }

}
