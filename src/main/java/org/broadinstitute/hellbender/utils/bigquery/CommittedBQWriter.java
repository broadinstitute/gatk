package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import com.google.protobuf.Descriptors;
import org.json.JSONArray;
import org.json.JSONObject;
import org.threeten.bp.Duration;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class CommittedBQWriter implements AutoCloseable {
    protected BigQueryWriteClient bqWriteClient;
    protected WriteStream writeStream;
    protected JsonStreamWriter writer;
    protected TableName parentTable;
    protected WriteStream.Type steamType;
    protected int BATCH_SIZE = 10000;
    protected JSONArray jsonArr = new JSONArray();

    protected CommittedBQWriter(String projectId, String datasetName, String tableName, WriteStream.Type type) {
        this.parentTable = TableName.of(projectId, datasetName, tableName);
        this.steamType = type;
    }

    protected void createStream() throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        if (bqWriteClient == null) {
            BigQueryWriteSettings.Builder bigQueryWriteSettingsBuilder = BigQueryWriteSettings.newBuilder();
            bigQueryWriteSettingsBuilder
                    .createWriteStreamSettings()
                    .setRetrySettings(
                            bigQueryWriteSettingsBuilder
                                    .createWriteStreamSettings()
                                    .getRetrySettings()
                                    .toBuilder()
                                    .setInitialRetryDelay(Duration.ofSeconds(5))
                                    .setTotalTimeout(Duration.ofSeconds(30))
                                    .setRetryDelayMultiplier(1.5)
                                    .setMaxAttempts(3)
                                    .build());

            BigQueryWriteSettings bigQueryWriteSettings = bigQueryWriteSettingsBuilder.build();
            bqWriteClient = BigQueryWriteClient.create(bigQueryWriteSettings);
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

    public void addJsonRow(JSONObject row) throws Descriptors.DescriptorValidationException, ExecutionException, InterruptedException, IOException {
        if (writer == null) {
            createStream();
        }
        jsonArr.put(row);

        if (jsonArr.length() >= BATCH_SIZE) {
            writeJsonArray();
        }
    }

    protected void writeJsonArray() throws ExecutionException, InterruptedException {
        ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
        future.get();
        jsonArr = new JSONArray();
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
