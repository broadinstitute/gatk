package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.concurrent.ExecutionException;

public class PendingBQWriter {
    private BigQueryWriteClient bqWriteClient;
    private WriteStream writeStream;
    protected JsonStreamWriter writer;
    private TableName parentTable;

    public PendingBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Exception {
        this.bqWriteClient = bqWriteClient;
        writeStream = WriteStream.newBuilder().setType(WriteStream.Type.PENDING).build();
        parentTable = TableName.of(projectId, datasetName, tableName);
        CreateWriteStreamRequest createWriteStreamRequest =
                CreateWriteStreamRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .setWriteStream(writeStream)
                        .build();
        WriteStream writeStream = bqWriteClient.createWriteStream(createWriteStreamRequest);
        writer = JsonStreamWriter.newBuilder(writeStream.getName(), writeStream.getTableSchema()).build();
    }

    public void commitWriteStreams() {
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
                System.out.println(err.getErrorMessage());
            }
            throw new RuntimeException("Error committing the streams");
        }
        System.out.println("Appended and committed records successfully.");
    }

    public void addJsonRow(JSONObject row) throws InterruptedException, ExecutionException {
        JSONArray jsonArr = new JSONArray();
        jsonArr.put(row);
        ApiFuture<AppendRowsResponse> future = writer.append(jsonArr, row.getLong(VetFieldEnum.location.toString()));
        AppendRowsResponse response = future.get();
    }

    public void close() {
        writer.close();
    }
}
