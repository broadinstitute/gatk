package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.TableResult;
import com.google.cloud.bigquery.storage.v1beta2.AppendRowsResponse;
import com.google.cloud.bigquery.storage.v1beta2.JsonStreamWriter;
import com.google.cloud.bigquery.storage.v1beta2.TableFieldSchema;
import com.google.cloud.bigquery.storage.v1beta2.TableName;
import com.google.cloud.bigquery.storage.v1beta2.TableSchema;
import com.google.protobuf.Descriptors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.time.Instant;
import java.util.concurrent.ExecutionException;

public class LoadStatus {
    static final Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.tools.gvs.ingest.LoadStatus.class);

    private enum LoadStatusValues { STARTED, FINISHED };
    public enum LoadState { NONE, PARTIAL, COMPLETE };

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
                " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId;

        TableResult result = BigQueryUtils.executeQuery(projectID, query, true, null);

        int startedCount = 0;
        int finishedCount = 0;
        for ( final FieldValueList row : result.iterateAll() ) {
            final String status = row.get(0).getStringValue();
            if (LoadStatusValues.STARTED.toString().equals(status)) {
                startedCount++;
            } else if (LoadStatusValues.FINISHED.toString().equals(status)) {
                finishedCount++;
            }

        }

        // if fully loaded, exit successfully!
        if (startedCount == 1 && finishedCount == 1) {
            return LoadState.COMPLETE;
        }

        if (startedCount == 0 && finishedCount == 0) {
            return LoadState.NONE;
        }

        // otherwise if there are any records, return partial
        return LoadState.PARTIAL;
    }

    public void writeLoadStatusStarted(long sampleId) {
        writeLoadStatus(LoadStatusValues.STARTED, sampleId);
    }

    public void writeLoadStatusFinished(long sampleId) {
        writeLoadStatus(LoadStatusValues.FINISHED, sampleId);
    }

    protected void writeLoadStatus(LoadStatusValues status, long sampleId) {

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
            AppendRowsResponse response = future.get();

            logger.info("Load status " + status + " appended successfully");
        } catch (ExecutionException | InterruptedException | Descriptors.DescriptorValidationException | IOException e) {
            // If the wrapped exception is a StatusRuntimeException, check the state of the operation.
            // If the state is INTERNAL, CANCELLED, or ABORTED, you can retry. For more information, see:
            // https://grpc.github.io/grpc-java/javadoc/io/grpc/StatusRuntimeException.html
            throw new GATKException(e.getMessage());
        }
    }
}
