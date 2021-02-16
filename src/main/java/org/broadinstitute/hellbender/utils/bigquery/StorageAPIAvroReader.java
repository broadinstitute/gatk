package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.gax.rpc.ServerStream;
import com.google.cloud.bigquery.BigQueryOptions;
import com.google.cloud.bigquery.storage.v1beta1.*;
import com.google.cloud.bigquery.storage.v1beta1.ReadOptions.TableReadOptions.Builder;
import com.google.common.base.Preconditions;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.BinaryDecoder;
import org.apache.avro.io.DatumReader;
import org.apache.avro.io.DecoderFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class StorageAPIAvroReader implements GATKAvroReader {

    private static final Logger logger = LogManager.getLogger(StorageAPIAvroReader.class);

    private static int rowCount = 0;

    private BigQueryStorageClient client;

    private Iterator<Storage.ReadRowsResponse> serverStream;

    private org.apache.avro.Schema schema;

    private DatumReader<GenericRecord> datumReader;

    // Decoder object will be reused to avoid re-allocation and too much garbage
    // collection.
    private BinaryDecoder decoder = null;

    private AvroProto.AvroRows currentAvroRows;

    // GenericRecord object will be reused.
    private GenericRecord nextRow = null;

    public StorageAPIAvroReader(final TableReference tableRef) {
        this(tableRef, null, null);
    }

    public StorageAPIAvroReader(final TableReference tableRef, String parentProjectId) {
        this(tableRef, null, parentProjectId);
    }

    public StorageAPIAvroReader(final TableReference tableRef, final String rowRestriction, String parentProjectId) {

        try {
            logger.info("Using Storage API from " + tableRef + " with '" + rowRestriction + "'");

            this.client = BigQueryStorageClient.create();

            final String parent = String.format("projects/%s", parentProjectId == null || parentProjectId.isEmpty() ? tableRef.tableProject : parentProjectId);

            final TableReferenceProto.TableReference tableReference = TableReferenceProto.TableReference.newBuilder()
                    .setProjectId(tableRef.tableProject)
                    .setDatasetId(tableRef.tableDataset)
                    .setTableId(tableRef.tableName)
                    .build();

            Builder readOptions = ReadOptions.TableReadOptions.newBuilder()
                    .addAllSelectedFields(tableRef.fields);

            if (rowRestriction != null) {
                readOptions.setRowRestriction(rowRestriction);
            }
            final ReadOptions.TableReadOptions tableReadOptions = readOptions.build();

            final Storage.CreateReadSessionRequest.Builder builder = Storage.CreateReadSessionRequest.newBuilder()
                    .setParent(parent)
                    .setTableReference(tableReference)
                    .setReadOptions(tableReadOptions)
                    .setRequestedStreams(1)
                    .setFormat(Storage.DataFormat.AVRO);

            final Storage.ReadSession session = client.createReadSession(builder.build());

            if (session.getStreamsCount() > 0) {

                this.schema = new org.apache.avro.Schema.Parser().parse(session.getAvroSchema().getSchema());

                this.datumReader = new GenericDatumReader<>(
                        new org.apache.avro.Schema.Parser().parse(session.getAvroSchema().getSchema()));

                // Use the first stream to perform reading.
                Storage.StreamPosition readPosition = Storage.StreamPosition.newBuilder()
                        .setStream(session.getStreams(0))
                        .build();

                Storage.ReadRowsRequest readRowsRequest = Storage.ReadRowsRequest.newBuilder()
                        .setReadPosition(readPosition)
                        .build();

                this.serverStream = client.readRowsCallable().call(readRowsRequest).iterator();

                loadNextRow();
            }
        } catch ( IOException e ) {
            throw new GATKException("I/O Error", e);
        }
    }

    private void loadNextRow() {
        try {
            if ( decoder != null && ! decoder.isEnd() ) {
                nextRow = datumReader.read(null, decoder);
            } else {
                fetchNextAvroRows();

                if ( decoder != null && ! decoder.isEnd() ) {
                    nextRow = datumReader.read(null, decoder);
                } else {
                    nextRow = null; // end of traversal
                }
            }
        } catch ( IOException e ) {
            throw new GATKException("I/O error", e);
        }
    }

    private void fetchNextAvroRows() {
        if ( serverStream.hasNext() ) {
            currentAvroRows = serverStream.next().getAvroRows();
            decoder = DecoderFactory.get()
                    .binaryDecoder(currentAvroRows.getSerializedBinaryRows().toByteArray(), decoder);
        } else {
            currentAvroRows = null;
            decoder = null;
        }
    }

    @Override
    public org.apache.avro.Schema getSchema() {
        return schema;
    }

    @Override
    public Iterator<GenericRecord> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return nextRow != null;
    }

    @Override
    public GenericRecord next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when ! hasNext()");
        }

        final GenericRecord recordToReturn = nextRow;
        loadNextRow();
        return recordToReturn;
    }

    @Override
    public void close() {
        client.shutdownNow();
        /* TODO: do we need to wait for termination here?
        boolean terminated = false;
        while ( ! terminated ) {
            try {
                terminated = client.awaitTermination(100, TimeUnit.MILLISECONDS);
            }
            catch ( InterruptedException e ) {
                throw new GATKException("Interrupted during shutdown", e);
            }
        } */
        client.close();
    }
}
