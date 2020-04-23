package org.broadinstitute.hellbender.utils.localsort;

import org.apache.avro.Schema;
import org.apache.avro.file.DataFileStream;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class AvroSortingCollectionCodec implements EvoquerSortingCollection.Codec<GenericRecord> {

    private Schema schema;
    private DataFileWriter<GenericRecord> outputWriter;
    private DataFileStream<GenericRecord> inputStream;

    public AvroSortingCollectionCodec(final Schema schema) {
        this.schema = schema;
    }

    @Override
    public void setOutputStream( OutputStream os ) {
        this.outputWriter = new DataFileWriter<>(new GenericDatumWriter<>(schema));
        try {
            outputWriter.create(schema, os);
        } catch ( IOException e ) {
            throw new GATKException("Error initializing Avro output stream in SortingCollection", e);
        }
    }

    @Override
    public void flushOutput() {
        if ( outputWriter == null ) {
            throw new IllegalStateException("flushOutput() called when there is no OutputWriter");
        }

        try {
            outputWriter.close(); // This does a flush()
            outputWriter = null;
        } catch ( IOException e ) {
            throw new GATKException("Error closing Avro output writer in SortingCollection", e);
        }
    }

    @Override
    public void setInputStream( InputStream is ) {
        try {
            this.inputStream = new DataFileStream<>(is, new GenericDatumReader<>(schema));
        }
        catch ( IOException e ) {
            throw new GATKException("Error initializing Avro input stream in SortingCollection", e);
        }
    }

    @Override
    public void encode( GenericRecord val ) {
        try {
            outputWriter.append(val);
        } catch ( IOException e ) {
            throw new GATKException("Error writing to Avro output stream in SortingCollection", e);
        }
    }

    @Override
    public GenericRecord decode() {
        if ( inputStream == null ) {
            throw new IllegalStateException("decode() called without an inputStream");
        }

        if ( ! inputStream.hasNext() ) {
            return null;
        }

        return inputStream.next();
    }

    @Override
    public EvoquerSortingCollection.Codec<GenericRecord> clone() {
        return new AvroSortingCollectionCodec(schema);
    }
}
