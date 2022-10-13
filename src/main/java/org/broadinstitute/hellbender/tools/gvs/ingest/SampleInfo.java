package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.client.util.ExponentialBackOff;
import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.storage.v1beta2.*;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryResultAndStatistics;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;

import static org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils.extractCausalStatusRuntimeExceptionOrThrow;

public class SampleInfo {
    static final Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.tools.gvs.ingest.SampleInfo.class);

    private final String projectID;
    private final String datasetName;
    private final String sampleInfoTableName;
    private final TableName sampleInfoTable;


    public SampleInfo(String projectID, String datasetName, String sampleInfoTableName) {
        this.projectID = projectID;
        this.datasetName = datasetName;
        this.sampleInfoTableName = sampleInfoTableName;
        this.sampleInfoTable = TableName.of(projectID, datasetName, sampleInfoTableName);
    }

    public TableSchema getSampleInfoTableSchema() {
        TableSchema.Builder builder = TableSchema.newBuilder();
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_NAME_FIELD_NAME).setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_ID_FIELD_NAME).setType(TableFieldSchema.Type.INT64).setMode(TableFieldSchema.Mode.NULLABLE).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME).setType(TableFieldSchema.Type.BOOL).setMode(TableFieldSchema.Mode.NULLABLE).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_INFO_IS_CONTROL_FIELD_NAME).setType(TableFieldSchema.Type.BOOL).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_INFO_WITHDRAWN_FIELD_NAME).setType(TableFieldSchema.Type.TIMESTAMP).setMode(TableFieldSchema.Mode.NULLABLE).build()
        );
        return builder.build();
    }

    public void setSampleInfoIsLoaded(long sampleId) {
        final String query = "UPDATE `" + projectID + "." + datasetName + "." + sampleInfoTableName + "` " +
                " SET " + SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME + " = TRUE " +
                " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId;

        BigQueryUtils.executeQuery(projectID, query, false, null);
    }

    public boolean isLoaded(long sampleId) {
        final String query = "SELECT " + SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME +
                " FROM `" + projectID + "." + datasetName + "." + sampleInfoTableName + "` " +
                " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId;

        boolean isLoaded = false;
        BigQueryResultAndStatistics results = BigQueryUtils.executeQuery(projectID, query, false, null);
        if (results.result.getTotalRows() > 1) {
            // TODO - is it overkill to throw an exception here?
            throw new GATKException("Found more than 1 row in '" + projectID + "." + datasetName + "." + sampleInfoTableName + "' for sampleId: " + sampleId);
        }
        for ( final FieldValueList row : results.result.iterateAll() ) {
            isLoaded = !row.get(0).isNull() && row.get(0).getBooleanValue();
        }

        return isLoaded;
    }

    public void writeSampleInfoUsingDML(long sampleId, String sampleName, Boolean is_loaded, boolean is_control, Long withdrawnDate) {
        withdrawnDate = withdrawnDate == null ? null : withdrawnDate * 1000L;
        final String query = "INSERT INTO `" + projectID + "." + datasetName + "." + sampleInfoTableName + "` " +
                " (" + SchemaUtils.SAMPLE_NAME_FIELD_NAME + ", " +
                SchemaUtils.SAMPLE_ID_FIELD_NAME + ", " +
                SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME + ", " +
                SchemaUtils.SAMPLE_INFO_IS_CONTROL_FIELD_NAME + ", " +
                SchemaUtils.SAMPLE_INFO_WITHDRAWN_FIELD_NAME + ") VALUES (" +
                sampleName + ", " +
                sampleId + ", " +
                is_loaded + ", " +
                is_control + ", " +
                withdrawnDate + ") ";
        BigQueryUtils.executeQuery(projectID, query, false, null);
    }

    public void writeSampleInfo(long sampleId, String sampleName, Boolean is_loaded, boolean is_control, Long withdrawnDate) {
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
                         JsonStreamWriter.newBuilder(sampleInfoTable.toString(), getSampleInfoTableSchema()).build()) {

                // Create a JSON object that is compatible with the table schema.
                JSONArray jsonArr = new JSONArray();
                JSONObject jsonObject = new JSONObject();
                jsonObject.put(SchemaUtils.SAMPLE_NAME_FIELD_NAME, sampleName);
                jsonObject.put(SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
                jsonObject.put(SchemaUtils.SAMPLE_INFO_IS_LOADED_FIELD_NAME, is_loaded);
                jsonObject.put(SchemaUtils.SAMPLE_INFO_IS_CONTROL_FIELD_NAME, is_control);
                jsonObject.put(SchemaUtils.SAMPLE_INFO_WITHDRAWN_FIELD_NAME, withdrawnDate == null ? null : withdrawnDate * 1000L); // google wants this in microseconds since epoch...
                jsonArr.put(jsonObject);

                ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
                future.get();

                logger.info("Sample Info for sample name: " + sampleName + " loaded successfully");
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
