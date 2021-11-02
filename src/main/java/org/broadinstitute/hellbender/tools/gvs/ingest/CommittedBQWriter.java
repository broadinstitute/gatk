package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import com.google.protobuf.Descriptors;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;
import scala.collection.immutable.Range;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class CommittedBQWriter {
    static final Logger logger = LogManager.getLogger(CommittedBQWriter.class);

    protected BigQueryWriteClient bqWriteClient;
    protected WriteStream writeStream;
    protected JsonStreamWriter writer;
    protected TableName parentTable;
    protected static int MAX_RETRIES = 3;

    public CommittedBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        this(bqWriteClient, projectId, datasetName, tableName, WriteStream.Type.COMMITTED);
    }

    protected CommittedBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName, WriteStream.Type type) throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        this.bqWriteClient = bqWriteClient;
        WriteStream writeStreamConfig = WriteStream.newBuilder().setType(type).build();
        parentTable = TableName.of(projectId, datasetName, tableName);
        CreateWriteStreamRequest createWriteStreamRequest =
                CreateWriteStreamRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .setWriteStream(writeStreamConfig)
                        .build();
        writeStream = bqWriteClient.createWriteStream(createWriteStreamRequest);
        writer = JsonStreamWriter.newBuilder(writeStream.getName(), writeStream.getTableSchema()).build();
    }

    public AppendRowsResponse addJsonRow(JSONObject row, int offset) throws InterruptedException, ExecutionException {
        JSONArray jsonArr = new JSONArray();
        jsonArr.put(row);

        return addJsonRows(jsonArr, offset);
    }

    protected AppendRowsResponse addJsonRows(JSONArray jsonArray, int offset) throws InterruptedException, ExecutionException {
        int retryCount = 0;
        AppendRowsResponse response = null;
        while (retryCount < MAX_RETRIES && response == null) {
            try {
                ApiFuture<AppendRowsResponse> future = writer.append(jsonArray, offset);
                response = future.get();
            } catch (StatusRuntimeException ex) {
                Status status = ex.getStatus();
                if (status == Status.ABORTED || status == Status.INTERNAL || status == Status.CANCELLED) {
                    logger.warn("Caught exception " + ex + "\n Retrying. " + (MAX_RETRIES - retryCount - 1) + " retries remaining.");
                }
            }
            retryCount ++;
        }
        return response;
    }


    public void close() {
        writer.close();
    }
}
