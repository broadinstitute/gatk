package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.concurrent.ExecutionException;

public class PendingBQWriter {
    static final Logger logger = LogManager.getLogger(PendingBQWriter.class);

    private BigQueryWriteClient bqWriteClient;
    private WriteStream writeStream;
    protected JsonStreamWriter writer;
    private TableName parentTable;

    public PendingBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Exception {
        this.bqWriteClient = bqWriteClient;
        WriteStream writeStreamConfig = WriteStream.newBuilder().setType(WriteStream.Type.PENDING).build();
        parentTable = TableName.of(projectId, datasetName, tableName);
        CreateWriteStreamRequest createWriteStreamRequest =
                CreateWriteStreamRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .setWriteStream(writeStreamConfig)
                        .build();
        writeStream = bqWriteClient.createWriteStream(createWriteStreamRequest);
        writer = JsonStreamWriter.newBuilder(writeStream.getName(), writeStream.getTableSchema()).build();
    }

    public void commitWriteStreams() {
        FinalizeWriteStreamResponse finalizeResponse =
                bqWriteClient.finalizeWriteStream(writeStream.getName());
        logger.info("Rows written: " + finalizeResponse.getRowCount());

        BatchCommitWriteStreamsRequest commitRequest =
                BatchCommitWriteStreamsRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .addWriteStreams(writeStream.getName())
                        .build();
        BatchCommitWriteStreamsResponse commitResponse =
                bqWriteClient.batchCommitWriteStreams(commitRequest);
        // If the response does not have a commit time, it means the commit operation failed.
        if (commitResponse.hasCommitTime() == false) {
            for (StorageError err : commitResponse.getStreamErrorsList()) {
                logger.error(err.getErrorMessage());
            }
            throw new RuntimeException("Error committing the streams");
        }
        logger.info("Appended and committed records successfully.");
    }

    public void addJsonRow(JSONObject row) throws InterruptedException, ExecutionException {
        JSONArray jsonArr = new JSONArray();
        jsonArr.put(row);
        // TODO figure out what to use as the offset
        ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
        AppendRowsResponse response = future.get();
    }

    public void close() {
        writer.close();
    }
}
