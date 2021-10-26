package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.ExecutionException;

public class PetJsonWriter {
    private JsonStreamWriter writer;
    private JSONArray jsonArr = new JSONArray();
    private final static char SEPARATOR = IngestConstants.SEPARATOR;

    public PetJsonWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Exception {
        WriteStream stream = WriteStream.newBuilder().setType(WriteStream.Type.COMMITTED).build();
        TableName parentTable = TableName.of(projectId, datasetName, tableName);
        CreateWriteStreamRequest createWriteStreamRequest =
                CreateWriteStreamRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .setWriteStream(stream)
                        .build();
        WriteStream writeStream = bqWriteClient.createWriteStream(createWriteStreamRequest);
        JsonStreamWriter writer =
                JsonStreamWriter.newBuilder(writeStream.getName(), writeStream.getTableSchema())
                        .build();
    }

    public void addRow(long location, long sampleId, String state) throws IOException {
        JSONObject record = new JSONObject();

        record.put(SchemaUtils.LOCATION_FIELD_NAME, location);
        record.put(SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
        record.put(SchemaUtils.STATE_FIELD_NAME, state);

        jsonArr.put(record);
    }

    public void close() throws IOException {
        try {
            ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
            AppendRowsResponse response = future.get();
        } catch (Exception ex) {
            // for now wrap in IOException
            throw new IOException(ex);
        }
        writer.close();
    }}
