package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import java.io.Closeable;
import java.io.IOException;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;

import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.io.DatumWriter;

public class PetAvroWriter implements Closeable {
    private DataFileWriter<GenericData.Record> writer;
    private Schema schema;

    public PetAvroWriter(String outputFile) throws IOException{
        this.schema = SchemaBuilder.record("pet")
                                    .namespace("org.broadinstitute.dsp")
                                    .fields().requiredLong("location").requiredInt("sample_id").requiredString("state")
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

    public void addRow(long location, long sampleId, String state) throws IOException {
        GenericData.Record record = new GenericData.Record(schema);
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("state", state);
        writer.append(record);
    }

    public void close() throws IOException {
        writer.close();
    }
}
