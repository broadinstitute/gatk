package org.broadinstitute.hellbender.utils.gvs.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.column.ParquetProperties;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.io.OutputFile;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;


import java.io.IOException;

public class GvsVariantParquetFileWriter extends ParquetWriter<JSONObject> {

//    GvsVariantParquetFileWriter(
//            OutputFile file,
//            ParquetFileWriter.Mode mode,
//            GVSVariantWriteSupport writeSupport,
//            CompressionCodecName compressionCodecName,
//            long rowGroupSize,
//            boolean validating,
//            Configuration conf,
//            int maxPaddingSize,
//            ParquetProperties encodingProps,
//            FileEncryptionProperties encryptionProperties) throws IOException {
//        super(file, mode, writeSupport, compressionCodecName, rowGroupSize, validating, conf, maxPaddingSize, encodingProps, encryptionProperties);
//    }

    /**
     * This is very deprecated, and we'll need to figure out how to do this from a builder once it works!
     * @param file
     * @param schema
     * @param enableDictionary
     * @param codecName
     * @throws IOException
     */
    public GvsVariantParquetFileWriter(
            Path file,
            MessageType schema,
            boolean enableDictionary,
            CompressionCodecName codecName
    ) throws IOException {
        super(file, new GvsVariantWriteSupport(schema), codecName, DEFAULT_BLOCK_SIZE, DEFAULT_PAGE_SIZE, enableDictionary, false);
    }

    GvsVariantParquetFileWriter(
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

//        @Override
//        protected GVSVariantWriteSupport getWriteSupport(Configuration conf) {
//            return getWriteSupport((Configuration) null);
//        }

        @Override
        protected GvsVariantWriteSupport getWriteSupport(Configuration conf) {
            return new GvsVariantWriteSupport(schema);
        }

//        @Override
//        public Builder withExtraMetaData(Map<String, String> extraMetaData) {
//            return super.withExtraMetaData(extraMetaData);
//        }



    }

}

/*
public class GvsVariantParquetFileWriter extends ParquetWriter<JSONObject> {
    @SuppressWarnings("deprecation")
    public GvsVariantParquetFileWriter(
            Path file,
            MessageType schema,
            boolean enableDictionary,
            CompressionCodecName codecName
    ) throws IOException {
        // update this to the new constructor soon
        super(file, new GVSVariantWriteSupport(schema), codecName, DEFAULT_BLOCK_SIZE, DEFAULT_PAGE_SIZE, enableDictionary, false);
    }
}*/
