package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import org.broadinstitute.hellbender.tools.gvs.ingest.RefRangesWriter;
import org.json.JSONObject;

import java.io.IOException;

public class RefRangesBQWriter extends AbstractBQWriter implements RefRangesWriter {

    public RefRangesBQWriter(String projectId, String datasetName, String tableName) throws IOException {
        super(projectId, datasetName, tableName);
    }

    @Override
    public void write(long location, long sampleId, int length, String state) throws IOException {
        write(createJsonRow(location, sampleId, length, state));
    }

    @Override
    public void writeCompressed(long packedData, long sampleId) throws IOException {
        write(createCompressedJsonRow(packedData, sampleId));
    }

    private JSONObject createJsonRow(long location, long sampleId, int length, String state) {
        JSONObject record = new JSONObject();
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("length", length);
        record.put("state", state);
        return record;
    }

    private JSONObject createCompressedJsonRow(long packedData, long sampleId) {
        JSONObject record = new JSONObject();
        record.put("packed_ref_data", packedData);
        record.put("sample_id", sampleId);
        return record;
    }
}
