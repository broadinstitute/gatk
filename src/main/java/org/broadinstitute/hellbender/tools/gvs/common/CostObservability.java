package org.broadinstitute.hellbender.tools.gvs.common;

import com.google.api.client.util.ExponentialBackOff;
import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.util.Date;

import static org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils.extractCausalStatusRuntimeExceptionOrThrow;

public class CostObservability {
    static final Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.tools.gvs.common.CostObservability.class);

    private final TableName costObservabilityTable;

    public CostObservability(String projectID, String datasetName, String costObservabilityTableName) {
        this.costObservabilityTable = TableName.of(projectID, datasetName, costObservabilityTableName);
    }

    public TableSchema getCostObservabilityTableSchema() {
        TableSchema.Builder builder = TableSchema.newBuilder();
        builder.addFields(
                TableFieldSchema.newBuilder().setName("call_set_identifier").setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("step").setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("call").setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.NULLABLE).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("shard_identifier").setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.NULLABLE).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("call_start_timestamp").setType(TableFieldSchema.Type.TIMESTAMP).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("event_timestamp").setType(TableFieldSchema.Type.TIMESTAMP).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("event_key").setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName("event_bytes").setType(TableFieldSchema.Type.INT64).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        return builder.build();
    }

    public void writeCostObservability(String callSetIdentifier, String step, String call, String shardIdentifier,
                                       Date callStartTimestamp, Date eventTimestamp, String eventKey, long eventBytes) {
        final ExponentialBackOff backoff = new ExponentialBackOff.Builder().
                setInitialIntervalMillis(2000).
                setMaxIntervalMillis(30000).
                setMultiplier(2).
                setRandomizationFactor(0.5).
                build();
        int retryCount = 0;
        int maxRetries = 3;

        while (true) {
            // This uses the _default stream since it (a) commits immediately and (b) doesn't count
            // towards the CreateStreamWriter quota
            try (JsonStreamWriter writer =
                         JsonStreamWriter.newBuilder(costObservabilityTable.toString(), getCostObservabilityTableSchema()).build()) {

                // Create a JSON object that is compatible with the table schema.
                JSONArray jsonArr = new JSONArray();
                JSONObject jsonObject = new JSONObject();
                jsonObject.put("call_set_identifier", callSetIdentifier);
                jsonObject.put("step", step);
                jsonObject.put("call", call);
                jsonObject.put("shard_identifier", shardIdentifier);
                jsonObject.put("call_start_timestamp", callStartTimestamp.getTime() * 1000); // google wants this in microseconds since epoch...
                jsonObject.put("event_timestamp", eventTimestamp.getTime() * 1000); // google wants this in microseconds since epoch...
                jsonObject.put("event_key", eventKey);
                jsonObject.put("event_bytes", eventBytes);
                jsonArr.put(jsonObject);

                ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
                future.get();

                logger.info("Cost Observability for " + callSetIdentifier + "." + step + " appended successfully");
                break;
            } catch (Exception e) {
                @SuppressWarnings("ThrowableNotThrown")
                StatusRuntimeException se = extractCausalStatusRuntimeExceptionOrThrow(e);

                if (retryCount >= maxRetries) {
                    throw new GATKException("Caught exception writing to BigQuery and " + maxRetries + " write retries are exhausted", e);
                }

                switch (se.getStatus().getCode()) {
                    case ALREADY_EXISTS:
                        // This is okay, no need to retry
                        break;
                    case INVALID_ARGUMENT:
                    case NOT_FOUND:
                    case OUT_OF_RANGE:
                    case PERMISSION_DENIED:
                        throw new GATKException("Caught non-retryable StatusRuntimeException based exception", e);
                    default:
                        try {
                            logger.warn("Caught exception writing to BigQuery, " + (maxRetries - retryCount - 1) + " retries remaining.", e);
                            long backOffMillis = backoff.nextBackOffMillis();
                            //noinspection BusyWait
                            Thread.sleep(backOffMillis);
                            retryCount++;
                        } catch (final IOException | InterruptedException ie) {
                            throw new GATKException("Error attempting to sleep between retry attempts", ie);
                        }
                }
            }
        }
    }
}
