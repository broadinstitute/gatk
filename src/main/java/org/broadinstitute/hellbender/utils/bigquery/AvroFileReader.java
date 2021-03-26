package org.broadinstitute.hellbender.utils.bigquery;

import org.apache.avro.file.DataFileStream;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.BinaryDecoder;
import org.apache.avro.io.DatumReader;
import org.apache.avro.io.DecoderFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bigquery.GATKAvroReader;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.apache.avro.Schema;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;

public class AvroFileReader implements GATKAvroReader {
    private DatumReader<GenericRecord> datumReader;
    private DataFileStream<GenericRecord> dataFileStream;
    private org.apache.avro.Schema schema;

    // Decoder object will be reused to avoid re-allocation and too much garbage collection.
    private BinaryDecoder decoder = null;

    // GenericRecord object will be reused.
    private GenericRecord nextRow = null;

    public AvroFileReader(final String avroFileURI ) {
        try {
            final Path avroFilePath = IOUtils.getPath(avroFileURI);

            datumReader = new GenericDatumReader<>();
            dataFileStream = new DataFileStream<>(Files.newInputStream(avroFilePath), datumReader);
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
