package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import java.io.Closeable;
import java.io.IOException;

import org.apache.hadoop.fs.Path;
import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.generic.GenericData;
import org.apache.hadoop.conf.Configuration;
import org.apache.parquet.avro.AvroParquetWriter;
import org.apache.parquet.column.ParquetProperties.WriterVersion;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.ParquetFileWriter.Mode;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;

public class PetParquetWriter implements Closeable {
    private ParquetWriter<GenericData.Record> writer;
    private Schema schema;

    public PetParquetWriter(String outputFile) throws IOException{
        this.schema = SchemaBuilder.record("pet")
                                    .namespace("org.broadinstitute.dsp")
                                    .fields().requiredInt("location").requiredInt("sample_id").requiredString("state")
                                    .endRecord();
                                    
        writer = AvroParquetWriter.
                            <GenericData.Record>builder(new Path(outputFile))
                            .withRowGroupSize(ParquetWriter.DEFAULT_BLOCK_SIZE)
                            .withPageSize(ParquetWriter.DEFAULT_PAGE_SIZE)
                            .withSchema(schema)
                            .withWriteMode(Mode.OVERWRITE)
                            .withConf(new Configuration())
                            .withCompressionCodec(CompressionCodecName.SNAPPY)
                            .withValidation(false)
                            .withDictionaryEncoding(true)
                            .withWriterVersion(WriterVersion.PARQUET_1_0)
                            .build();
        
    }

    public void addRow(long location, long sampleId, String state) throws IOException {
        GenericData.Record record = new GenericData.Record(schema);
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("state", state);
        writer.write(record);
    }

    public void close() throws IOException {
        writer.close();
    }
}
