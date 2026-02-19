package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class SamplePloidyBQWriter extends SamplePloidyWriter {

    private PendingBQWriter pendingBQWriter = null;

    public SamplePloidyBQWriter(String projectId, String datasetName, String samplePloidyTableName) {
        pendingBQWriter = new PendingBQWriter(projectId, datasetName, samplePloidyTableName);
    }

    @Override
    public void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException {
        JSONObject jsonRow = createJson(sampleId, chromosome, ploidy);
        try {
            this.pendingBQWriter.addJsonRow(jsonRow);
        } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException e) {
            throw new IOException(e);
        }
    }

    @Override
    public void commitData() {
        super.commitData();
        if (pendingBQWriter != null) {
            pendingBQWriter.flushBuffer();
            pendingBQWriter.commitWriteStreams();
        }
    }

    @Override
    public void close() throws IOException {
        if (pendingBQWriter != null) {
            pendingBQWriter.close();
        }
    }

    private JSONObject createJson(Long sampleId, Long chromosome, Integer ploidy) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("chromosome", chromosome);
        record.put("ploidy", ploidy);
        return record;
    }
}
