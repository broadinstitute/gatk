package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.concurrent.ExecutionException;

public class PetBQWriter extends PendingBQWriter {

    public PetBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Exception {
        super(bqWriteClient, projectId, datasetName, tableName);
    }

    public void addRow(long location, long sampleId, String state) throws InterruptedException, ExecutionException {
        JSONObject record = new JSONObject();
        JSONArray jsonArr = new JSONArray();

        record.put(SchemaUtils.LOCATION_FIELD_NAME, location);
        record.put(SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
        record.put(SchemaUtils.STATE_FIELD_NAME, state);

        jsonArr.put(record);

        ApiFuture<AppendRowsResponse> future = writer.append(jsonArr, location);
        AppendRowsResponse response = future.get();
    }

}
