package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.client.util.ExponentialBackOff;
import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import com.google.protobuf.Descriptors;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

import static org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils.extractCausalStatusRuntimeExceptionOrThrow;


public class CommittedBQWriter implements AutoCloseable {
    static final Logger logger = LogManager.getLogger(CommittedBQWriter.class);
    protected BigQueryWriteClient bqWriteClient;
    protected WriteStream writeStream;
    protected JsonStreamWriter writer;
    protected TableName parentTable;
    protected int maxRetries = 3;
    protected WriteStream.Type streamType;
    protected int batchSize = 10000;
    protected JSONArray jsonArr = new JSONArray();

    private final ExponentialBackOff backoff = new ExponentialBackOff.Builder().
            setInitialIntervalMillis(2000).
            setMaxIntervalMillis(30000).
            // The default value is 900000 millis or 900 seconds or 15 minutes which is not long enough given how this
            // backoff is being used. The retry in this code is limited by number of attempts and not total time.
            setMaxElapsedTimeMillis(Integer.MAX_VALUE).
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

    protected void writeJsonArray() throws InterruptedException, IOException {
        writeJsonArray(0);
    }

    protected void writeJsonArray(int retryCount) throws InterruptedException, IOException {
        try {
            ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
            future.get();
            jsonArr = new JSONArray();
        } catch (Exception e) {
            StatusRuntimeException se = extractCausalStatusRuntimeExceptionOrThrow(e);
            // Google BigQuery write API error handling
            // https://cloud.google.com/bigquery/docs/write-api#error_handling
            // The comments in the case matching below are nearly all quotations from the documentation linked above.
            Status.Code code = se.getStatus().getCode();
            switch (code) {
                case ALREADY_EXISTS:
                    // ALREADY_EXISTS: The row was already written. This error can happen when you provide stream offsets.
                    // It indicates that a duplicate record was detected. You can safely ignore this error.

                    // We don't expect to see this, but to "ignore" we should not retry sending the same `jsonArr`.
                    jsonArr = new JSONArray();
                    return;
                case INVALID_ARGUMENT:
                case NOT_FOUND:
                case OUT_OF_RANGE:
                case PERMISSION_DENIED:
                    // Do not retry these:
                    //
                    // INVALID_ARGUMENT: Invalid argument. This error is an application error.
                    //
                    // NOT_FOUND: The stream or table was not found.
                    //
                    // OUT_OF_RANGE. The offset is beyond the current write offset. This error can happen if you provide
                    // stream offsets and a previous write operation failed. In that case, you can retry from the last
                    // successful write. This error can also happen if the application sets the wrong offset value.
                    //
                    // PERMISSION_DENIED. The application does not have permission to write to this table.
                    throw new GATKException("Caught StatusRuntimeException with unretryable status code " + code + ", throwing", e);
                default:
                    switch (code) {
                        case INTERNAL:
                        case CANCELLED:
                        case ABORTED:
                            // INTERNAL, CANCELLED, or ABORTED: The operation could not be completed. You can safely retry the operation.
                            if (retryCount >= maxRetries) {
                                throw new GATKException("Caught exception writing to BigQuery and " + maxRetries + " write retries are exhausted", e);
                            }
                            logger.warn("Caught exception writing to BigQuery, " + (maxRetries - retryCount - 1) + " retries remaining.", e);
                            long backOffMillis = backoff.nextBackOffMillis();
                            logger.warn("Sleeping for {} milliseconds before retrying.", backOffMillis);
                            Thread.sleep(backOffMillis);

                            writeJsonArray(retryCount + 1);
                            break;
                        default:
                            // If you receive an error that's not listed above, then try to open a new connection by closing the
                            // writer object and creating a new instance.
                            //
                            // In practice with a PENDING WriteStream.Type the advice above does not work out;
                            // everything that has been `append`ed prior to the writer being closed is lost. Instead,
                            // just throw and let `maxRetries` start up a subsequent attempt from the beginning.
                            throw new GATKException("Caught StatusRuntimeException with status code "  + code + " which is not known to be retryable, throwing", e);
                    }
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
