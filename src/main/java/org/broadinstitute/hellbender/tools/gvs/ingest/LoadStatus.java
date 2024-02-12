package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.client.util.ExponentialBackOff;
import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1.*;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryResultAndStatistics;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.util.Set;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils.extractCausalStatusRuntimeExceptionOrThrow;

public class LoadStatus {
    static final Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.tools.gvs.ingest.LoadStatus.class);

    private enum LoadStatusValue { STARTED, HEADERS_LOADED, FINISHED }
    public static class LoadState {
        private final Set<LoadStatusValue> statusValues;

        public LoadState(Set<LoadStatusValue> statusValues) {
            this.statusValues = statusValues;
        }

        public boolean isComplete() {
            return statusValues.contains(LoadStatusValue.FINISHED);
        }

        public boolean areHeadersLoaded() {
            return statusValues.contains(LoadStatusValue.HEADERS_LOADED);
        }

        public boolean isStarted() {
            return statusValues.contains(LoadStatusValue.STARTED);
        }
    }

    private final String projectID;
    private final String datasetName;
    private final String loadStatusTableName;
    private final TableName loadStatusTable;

    public LoadStatus(String projectID, String datasetName, String loadStatusTableName) {
        this.projectID = projectID;
        this.datasetName = datasetName;
        this.loadStatusTableName = loadStatusTableName;
        this.loadStatusTable = TableName.of(projectID, datasetName, loadStatusTableName);
    }

    public TableSchema getLoadStatusTableSchema() {
        TableSchema.Builder builder = TableSchema.newBuilder();
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_ID_FIELD_NAME).setType(TableFieldSchema.Type.INT64).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.LOAD_STATUS_FIELD_NAME).setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.LOAD_STATUS_EVENT_TIMESTAMP_NAME).setType(TableFieldSchema.Type.TIMESTAMP).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        return builder.build();
    }

    public LoadState getSampleLoadState(long sampleId) {
        String query = "SELECT " + SchemaUtils.LOAD_STATUS_FIELD_NAME +
                " FROM `" + projectID + "." + datasetName + "." + loadStatusTableName + "` " +
                " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId +
                " ORDER BY " + SchemaUtils.LOAD_STATUS_EVENT_TIMESTAMP_NAME;

        BigQueryResultAndStatistics results = BigQueryUtils.executeQuery(projectID, query, true, null);

        Set<LoadStatusValue> loadStatusValues =
                results.result.streamAll().map(v -> LoadStatusValue.valueOf(v.get(0).getStringValue())).collect(Collectors.toSet());

        if (loadStatusValues.contains(LoadStatusValue.FINISHED) && !loadStatusValues.contains(LoadStatusValue.STARTED)) {
            // Should this be an error?
            logger.warn("Found Load Status 'FINISHED' where no previous Load Status: 'STARTED' for sample " + sampleId);
        }
        return new LoadState(loadStatusValues);
    }

    public void writeLoadStatusStarted(long sampleId) {
        writeLoadStatus(LoadStatusValue.STARTED, sampleId);
    }

    public void writeVCFHeadersLoaded(long sampleId) {
        writeLoadStatus(LoadStatusValue.HEADERS_LOADED, sampleId);
    }

    public void writeLoadStatusFinished(long sampleId) {
        writeLoadStatus(LoadStatusValue.FINISHED, sampleId);
    }

    protected void writeLoadStatus(LoadStatusValue status, long sampleId) {
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
                         JsonStreamWriter.newBuilder(loadStatusTable.toString(), getLoadStatusTableSchema()).build()) {

                // Create a JSON object that is compatible with the table schema.
                JSONArray jsonArr = new JSONArray();
                JSONObject jsonObject = new JSONObject();
                jsonObject.put("sample_id", sampleId);
                jsonObject.put("status", status.toString());
                jsonObject.put("event_timestamp", System.currentTimeMillis() * 1000L); // google wants this in microseconds since epoch...
                jsonArr.put(jsonObject);

                ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
                future.get();

                logger.info("Load status " + status + " appended successfully");
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
