package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.io.DatumWriter;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

public class RefRangesAvroWriter implements Closeable {
    private DataFileWriter<GenericData.Record> writer;
    private Schema schema;

    public RefRangesAvroWriter(String outputFile) throws IOException{
        this.schema = SchemaBuilder.record("ref_ranges")
                                    .namespace("org.broadinstitute.dsp")
                                    .fields().requiredLong("location").requiredInt("sample_id").requiredInt("length").requiredString("state")
                                    .endRecord();

        
        DatumWriter<GenericData.Record> datumWriter = new GenericDatumWriter<GenericData.Record>(schema);
		
        writer = new DataFileWriter<GenericData.Record>(datumWriter);
        writer.setCodec(CodecFactory.deflateCodec(5));

        try {
            writer.create(schema, new File(outputFile));
        } catch (IOException ioe) {
            throw new GATKException("Unable to create AvroWriter", ioe);
        }
	
    }

    public void write(long location, long sampleId, long length, String state) throws IOException {
        GenericData.Record record = new GenericData.Record(schema);
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("length", length);
        record.put("state", state);
        writer.append(record);
    }

    public void close() throws IOException {
        writer.close();
    }
}
