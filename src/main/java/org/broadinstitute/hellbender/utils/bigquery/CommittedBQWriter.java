package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.client.util.ExponentialBackOff;
import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import com.google.common.collect.ImmutableSet;
import com.google.protobuf.Descriptors;
import io.grpc.Status.Code;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

import static io.grpc.Status.Code.*;


public class CommittedBQWriter implements AutoCloseable {
    static final Logger logger = LogManager.getLogger(CommittedBQWriter.class);

    protected BigQueryWriteClient bqWriteClient;
    protected WriteStream writeStream;
    protected JsonStreamWriter writer;
    protected TableName parentTable;
    protected WriteStream.Type streamType;
    protected int maxRetries = 3;
    protected int batchSize = 10000;
    protected JSONArray jsonArr = new JSONArray();

    private final ExponentialBackOff backoff = new ExponentialBackOff.Builder().
            setInitialIntervalMillis(2000).
            setMaxIntervalMillis(30000).
            setMultiplier(2).
            setRandomizationFactor(0.5).
            build();

    protected CommittedBQWriter(String projectId, String datasetName, String tableName, WriteStream.Type type) {
        this.parentTable = TableName.of(projectId, datasetName, tableName);
        this.streamType = type;
    }

    protected void createStream() throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        if (bqWriteClient == null) {
            bqWriteClient = BigQueryWriteClient.create();
        }
        WriteStream writeStreamConfig = WriteStream.newBuilder().setType(streamType).build();
        CreateWriteStreamRequest createWriteStreamRequest =
                CreateWriteStreamRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .setWriteStream(writeStreamConfig)
                        .build();
        writeStream = bqWriteClient.createWriteStream(createWriteStreamRequest);
        writer = JsonStreamWriter.newBuilder(writeStream.getName(), writeStream.getTableSchema()).build();
    }

    public void addJsonRow(JSONObject row) throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        if (writer == null) {
            createStream();
        }
        jsonArr.put(row);

        if (jsonArr.length() >= batchSize) {
            writeJsonArray();
        }
    }

    protected void writeJsonArray() throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        writeJsonArray(0);
    }

    private StatusRuntimeException findCausalStatusRuntimeException(Exception e) {
        StatusRuntimeException se = null;
        Throwable t = e;
        while (true) {
            if (t == null) break;
            if (t instanceof StatusRuntimeException) {
                se = (StatusRuntimeException) t;
                break;
            }
            t = t.getCause();
        }
        return se;
    }

    protected AppendRowsResponse writeJsonArray(int retryCount) throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        AppendRowsResponse response;
        try {
            ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
            response = future.get();
            jsonArr = new JSONArray();
            return response;
        } catch (Exception e) {
            StatusRuntimeException se = findCausalStatusRuntimeException(e);
            if (se == null) {
                throw e;
            }

            Code code = se.getStatus().getCode();
            if (ImmutableSet.of(ABORTED, CANCELLED, INTERNAL, UNAVAILABLE).contains(code)) {
                if (retryCount >= maxRetries) {
                    throw new GATKException("Caught exception writing to BigQuery and " + maxRetries + " write retries are exhausted", se);
                }

                // .error during validation as I am suspicious about not seeing output
                logger.error("Caught exception writing to BigQuery, " + (maxRetries - retryCount - 1) + " retries remaining.", se);
                long backOffMillis = backoff.nextBackOffMillis();
                Thread.sleep(backOffMillis);
                createStream();
                return writeJsonArray(retryCount + 1);
            } else {
                throw se;
            }
        }
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
