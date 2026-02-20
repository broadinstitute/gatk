package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import org.broadinstitute.hellbender.tools.gvs.ingest.SamplePloidyWriter;
import org.json.JSONObject;

import java.io.IOException;

public class SamplePloidyBQWriter extends AbstractBQWriter implements SamplePloidyWriter {

    public SamplePloidyBQWriter(String projectId, String datasetName, String samplePloidyTableName) {
        super(projectId, datasetName, samplePloidyTableName);
    }

    @Override
    public void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException {
        write(createJson(sampleId, chromosome, ploidy));
    }

    private JSONObject createJson(Long sampleId, Long chromosome, Integer ploidy) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("chromosome", chromosome);
        record.put("ploidy", ploidy);
        return record;
    }
}
