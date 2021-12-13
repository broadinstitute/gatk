package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.storage.v1.AvroRows;
import com.google.cloud.bigquery.storage.v1.BigQueryReadClient;
import com.google.cloud.bigquery.storage.v1.CreateReadSessionRequest;
import com.google.cloud.bigquery.storage.v1.DataFormat;
import com.google.cloud.bigquery.storage.v1.ReadRowsRequest;
import com.google.cloud.bigquery.storage.v1.ReadSession;
import com.google.cloud.bigquery.storage.v1.ReadSession.TableReadOptions;
import com.google.cloud.bigquery.storage.v1.ReadSession.TableReadOptions.Builder;

import com.google.cloud.bigquery.storage.v1.ReadRowsResponse;
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

    private BigQueryReadClient client;

    private Iterator<ReadRowsResponse> serverStream;

    private org.apache.avro.Schema schema;

    private DatumReader<GenericRecord> datumReader;

    // Decoder object will be reused to avoid re-allocation and too much garbage
    // collection.
    private BinaryDecoder decoder = null;

    private AvroRows currentAvroRows;

    // GenericRecord object will be reused.
    private GenericRecord nextRow = null;

    private final ReadSession session;

    public StorageAPIAvroReader(final TableReference tableRef) {
        this(tableRef, null, null);
    }

    public StorageAPIAvroReader(final TableReference tableRef, String parentProjectId) {
        this(tableRef, null, parentProjectId);
    }

    public StorageAPIAvroReader(final TableReference tableRef, final String rowRestriction, String parentProjectId) {

        try {
            logger.info("Using Storage API from " + tableRef.getFQTableName() + " with '" + rowRestriction + "'");

            this.client = BigQueryReadClient.create();

            final String parent = String.format("projects/%s", parentProjectId == null || parentProjectId.isEmpty() ? tableRef.tableProject : parentProjectId);

            final String srcTable =
                    String.format(
                            "projects/%s/datasets/%s/tables/%s",
                            tableRef.tableProject, tableRef.tableDataset, tableRef.tableName);

            Builder readOptions =
                    ReadSession.TableReadOptions.newBuilder()
                            .addAllSelectedFields(tableRef.fields);

            if (rowRestriction != null) {
                readOptions.setRowRestriction(rowRestriction);
            }
            final TableReadOptions tableReadOptions = readOptions.build();

            // Start specifying the read session we want created.
            ReadSession.Builder sessionBuilder =
                    ReadSession.newBuilder()
                            .setTable(srcTable)
                            .setDataFormat(DataFormat.AVRO)
                            .setReadOptions(tableReadOptions);

            // Begin building the session creation request.
            CreateReadSessionRequest.Builder builder =
                    CreateReadSessionRequest.newBuilder()
                            .setParent(parent)
                            .setReadSession(sessionBuilder)
                            .setMaxStreamCount(1);

            this.session = client.createReadSession(builder.build());
            if (this.session.getStreamsCount() > 0) {

                this.schema = new org.apache.avro.Schema.Parser().parse(session.getAvroSchema().getSchema());

                this.datumReader = new GenericDatumReader<>(
                        new org.apache.avro.Schema.Parser().parse(session.getAvroSchema().getSchema()));

                logger.info("Storage API Session ID: " + session.getName());
                logger.info("Storage API Estimated Bytes Scanned for " + tableRef.getFQTableName() + " :" + getEstimatedTotalBytesScanned());

                // Use the first stream to perform reading.
                String streamName = session.getStreams(0).getName();

                ReadRowsRequest readRowsRequest = ReadRowsRequest.newBuilder()
                        .setReadStream(streamName)
                        .build();

                this.serverStream = client.readRowsCallable().call(readRowsRequest).iterator();

                loadNextRow();
            }
        } catch ( IOException e ) {
            throw new GATKException("I/O Error", e);
        }
    }

    public long getEstimatedTotalBytesScanned() {
        return session.getEstimatedTotalBytesScanned();
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
