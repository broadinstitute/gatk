package org.broadinstitute.hellbender.tools.evoquer;

import org.apache.avro.file.DataFileStream;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumReader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;

public class GCSAvroReader implements GATKAvroReader {

    private DatumReader<GenericRecord> datumReader;
    private DataFileStream<GenericRecord> dataFileStream;
    private org.apache.avro.Schema schema;

    public GCSAvroReader( final String avroFileURI ) {
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
    public org.apache.avro.Schema getSchema() {
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
