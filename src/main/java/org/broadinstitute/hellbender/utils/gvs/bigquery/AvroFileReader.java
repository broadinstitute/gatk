package org.broadinstitute.hellbender.utils.gvs.bigquery;

import org.apache.avro.Schema;
import org.apache.avro.file.DataFileStream;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumReader;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bigquery.GATKAvroReader;

import java.io.IOException;
import java.util.Iterator;

public class AvroFileReader implements GATKAvroReader {
    private final DataFileStream<GenericRecord> dataFileStream;
    private final org.apache.avro.Schema schema;

    public AvroFileReader(final GATKPath avroFilePath ) {
        try {
            DatumReader<GenericRecord> datumReader = new GenericDatumReader<>();
            dataFileStream = new DataFileStream<>(avroFilePath.getInputStream(), datumReader);
            schema = dataFileStream.getSchema();
        } catch (IOException e ) {
            throw new GATKException("I/O Error", e);
        }
    }

    @Override
    public Schema getSchema() {
        return schema;
    }

    @Override
    public boolean hasNext() {
        return dataFileStream.hasNext();
    }

    @Override
    public GenericRecord next() {
        return dataFileStream.next();
    }

    @Override
    public void close() {
        try {
            dataFileStream.close();
        } catch (IOException e ) {
            throw new GATKException("Error closing stream", e);
        }
    }

    @Override
    public Iterator<GenericRecord> iterator() {
        return this;
    }

}
