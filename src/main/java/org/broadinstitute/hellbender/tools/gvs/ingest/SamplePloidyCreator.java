package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.BitSet;
import java.util.Map;
import java.util.concurrent.ExecutionException;

public class SamplePloidyCreator {
    private static final Logger logger = LogManager.getLogger(SamplePloidyCreator.class);

    private PendingBQWriter samplePloidyBQJsonWriter = null;
    private static final String SAMPLE_PLOIDY_TABLE_NAME = "sample_chromosome_ploidy";

    private final Long sampleId;
    private final String projectId;
    private final String datasetName;


    public SamplePloidyCreator(Long sampleId, String projectId, String datasetName) {
        try {
            this.sampleId = sampleId;
            this.projectId = projectId;
            this.datasetName = datasetName;

            if (projectId == null || datasetName == null) {
                throw new UserException("Must specify project-id and dataset-name.");
            }
            samplePloidyBQJsonWriter = new PendingBQWriter(projectId, datasetName, SAMPLE_PLOIDY_TABLE_NAME);
        } catch (Exception e) {
            throw new UserException("Could not create VCF Header Scratch Table Writer", e);
        }
    }

    public void apply(Map<String, BitSet> ploidyData) throws IOException {
        for (final Map.Entry<String, BitSet> ploidyLine : ploidyData.entrySet()) {
            try {
                // make sure we don't have mixed ploidy in a single chromosome
                if (ploidyLine.getValue().length() > 1) {
                    throw new UserException("Detected mixed ploidy in sample "+this.sampleId+" on chromosome "+ploidyLine.getKey());
                }

                samplePloidyBQJsonWriter.addJsonRow(createJson(this.sampleId, SchemaUtils.encodeLocation(ploidyLine.getKey(),0), ploidyLine.getValue().nextSetBit(0)));
            } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                throw new IOException("BQ exception", ex);
            }
        }
    }

    public JSONObject createJson(Long sampleId, Long chromosome, Integer ploidy) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);
        record.put("chromosome", chromosome);
        record.put("ploidy", ploidy);
        return record;
    }


    public void commitData() {
        if (samplePloidyBQJsonWriter != null) {
            samplePloidyBQJsonWriter.flushBuffer();
            samplePloidyBQJsonWriter.commitWriteStreams();
        }
    }

    public void closeTool() {
        if (samplePloidyBQJsonWriter != null) {
            try {
                samplePloidyBQJsonWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close Sample Ploidy Table writer", e);
            }
        }
        if (samplePloidyBQJsonWriter != null) {
            samplePloidyBQJsonWriter.close();
        }
    }

}
