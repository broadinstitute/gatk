package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import com.google.protobuf.Descriptors;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class AbstractBQWriter {
    protected PendingBQWriter pendingBQWriter;

    AbstractBQWriter(String projectId, String datasetName, String tableName) {
        this.pendingBQWriter = new PendingBQWriter(projectId, datasetName, tableName);
    }

    protected void write(JSONObject object) throws IOException {
        try {
            pendingBQWriter.addJsonRow(object);
        } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException e) {
            throw new IOException(e);
        }
    }

    public void commitData() {
        pendingBQWriter.flushBuffer();
        pendingBQWriter.commitWriteStreams();
    }

    public void close() throws IOException {
        pendingBQWriter.close();
    }
}
