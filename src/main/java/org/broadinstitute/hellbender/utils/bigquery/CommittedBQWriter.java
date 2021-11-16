package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import com.google.protobuf.Descriptors;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class CommittedBQWriter {
    static final Logger logger = LogManager.getLogger(CommittedBQWriter.class);

    protected BigQueryWriteClient bqWriteClient;
    protected WriteStream writeStream;
    protected JsonStreamWriter writer;
    protected TableName parentTable;
    protected int MAX_RETRIES = 3;
    protected WriteStream.Type steamType;
    protected int BATCH_SIZE = 10000;
    protected JSONArray jsonArr = new JSONArray();


    public void setMaxRetries(int maxRetries) {
        MAX_RETRIES = maxRetries;
    }

    public void setBatchSize(int batchSize) {
        BATCH_SIZE = batchSize;
    }

    public CommittedBQWriter(String projectId, String datasetName, String tableName) {
        this(projectId, datasetName, tableName, WriteStream.Type.COMMITTED);
    }

    protected CommittedBQWriter(String projectId, String datasetName, String tableName, WriteStream.Type type) {
        this.parentTable = TableName.of(projectId, datasetName, tableName);
        this.steamType = type;
    }

    protected void createStream() throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        if (bqWriteClient == null) {
            bqWriteClient = BigQueryWriteClient.create();
        }
        WriteStream writeStreamConfig = WriteStream.newBuilder().setType(steamType).build();
        CreateWriteStreamRequest createWriteStreamRequest =
                CreateWriteStreamRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .setWriteStream(writeStreamConfig)
                        .build();
        writeStream = bqWriteClient.createWriteStream(createWriteStreamRequest);
        writer = JsonStreamWriter.newBuilder(writeStream.getName(), writeStream.getTableSchema()).build();
    }

    public AppendRowsResponse addJsonRow(JSONObject row) throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        AppendRowsResponse response = null;
        if (writer == null) {
            createStream();
        }
        jsonArr.put(row);

        if (jsonArr.length() >= BATCH_SIZE) {
            response = writeJsonArray();
        }
        return response;
    }

    protected AppendRowsResponse writeJsonArray() throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        return writeJsonArray(0);
    }

    protected AppendRowsResponse writeJsonArray(int retryCount) throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        AppendRowsResponse response = null;
        try {
            ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
            response = future.get();
            jsonArr = new JSONArray();
        } catch (StatusRuntimeException ex) {
            Status status = ex.getStatus();
            if (status == Status.ABORTED || status == Status.INTERNAL || status == Status.CANCELLED) {
                logger.warn("Caught exception " + ex + "\n Retrying. " + (MAX_RETRIES - retryCount - 1) + " retries remaining.");
            }
            if (retryCount < MAX_RETRIES) {
                createStream();
                response = writeJsonArray(retryCount + 1);
            } else {
                throw ex;
            }
        }
    return response;
    }


    public void close() {
        if (writer != null) {
            writer.close();
        }
        if (bqWriteClient != null) {
            bqWriteClient.close();
        }
    }
}
