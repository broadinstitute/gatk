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

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class CommittedBQWriter {
    static final Logger logger = LogManager.getLogger(CommittedBQWriter.class);

    protected BigQueryWriteClient bqWriteClient;
    protected WriteStream writeStream;
    protected JsonStreamWriter writer;
    protected TableName parentTable;
    protected static int MAX_RETRIES = 3;
    WriteStream.Type steamType;
    protected static int BATCH_SIZE = 1000;
    protected JSONArray jsonArr = new JSONArray();


    public CommittedBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        this(bqWriteClient, projectId, datasetName, tableName, WriteStream.Type.COMMITTED);
    }

    protected CommittedBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName, WriteStream.Type type) {
        this.bqWriteClient = bqWriteClient;
        this.parentTable = TableName.of(projectId, datasetName, tableName);
        this.steamType = type;
    }

    protected void createStream() throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
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
            response = writeJsonArray(0);
        }
        return response;
    }

    protected AppendRowsResponse writeJsonArray(int retryCount) throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        AppendRowsResponse response = null;
        try {
            ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
            logger.info("Wrote " + jsonArr.length() + " records");
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

        } catch (Exception ex) {
            logger.warn(ex);
            logger.info("creating new stream and retrying");
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
        if (jsonArr.length() > 0) {
            try {
                writeJsonArray(0);
            } catch (Exception ex) {
                logger.error("Caught exception writing last records on close", ex);
            }
        }
        writer.close();
    }
}
