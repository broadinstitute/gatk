package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ExecutionException;

public class VcfHeaderLineScratchCreator {
    private final Long sampleId;
    private final String projectId;
    private final String datasetName;

    private PendingBQWriter vcfHeaderBQJsonWriter = null;
    private static final String NON_SCRATCH_TABLE_NAME = "vcf_header_lines";
    private static final String SCRATCH_TABLE_NAME = "vcf_header_lines_scratch";

    private static boolean doScratchRowsExistFor(String projectId, String datasetName, String headerLineHash) {
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, SCRATCH_TABLE_NAME, "vcf_header_lines_hash", headerLineHash);
    }

    private static boolean doNonScratchRowsExistFor(String projectId, String datasetName, String headerLineHash) {
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, NON_SCRATCH_TABLE_NAME, "vcf_header_lines_hash", headerLineHash);
    }

    public VcfHeaderLineScratchCreator(Long sampleId, String projectId, String datasetName) {
        try {
            this.sampleId = sampleId;
            this.projectId = projectId;
            this.datasetName = datasetName;

            if (projectId == null || datasetName == null) {
                throw new UserException("Must specify project-id and dataset-name.");
            }
            vcfHeaderBQJsonWriter = new PendingBQWriter(projectId, datasetName, SCRATCH_TABLE_NAME);
        }
        catch (Exception e) {
            throw new UserException("Could not create VCF Header Scratch Table Writer", e);
        }

    }

    public void apply(Map<String, Boolean>  allLineHeaders) throws IOException {
        for (final Map.Entry<String, Boolean> headerChunk : allLineHeaders.entrySet()) {
            try {
                // if this header chunk has already been added to the scratch table, only add an association between the
                // sample_id and the hash, no need to rewrite the header chunk to the DB
                String chunkHash = Utils.calcMD5(headerChunk.getKey());
                Boolean isExpectedUnique = headerChunk.getValue();
                boolean vcfScratchHeaderRowsExist = doScratchRowsExistFor(this.projectId, this.datasetName, chunkHash);
                boolean vcfNonScratchHeaderRowsExist = doNonScratchRowsExistFor(this.projectId, this.datasetName, chunkHash);
                if (vcfScratchHeaderRowsExist || vcfNonScratchHeaderRowsExist) {
                    vcfHeaderBQJsonWriter.addJsonRow(createJson(this.sampleId, null, chunkHash, isExpectedUnique));
                }
                else {
                    vcfHeaderBQJsonWriter.addJsonRow(createJson(this.sampleId, headerChunk.getKey(), chunkHash, isExpectedUnique));
                }
            } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                throw new IOException("BQ exception", ex);
            }
        }
    }

    public JSONObject createJson(Long sampleId, String headerChunk, String headerHash, Boolean isExpectedUnique) {
        JSONObject record = new JSONObject();
        record.put("sample_id", sampleId);

        if (headerChunk != null) {
            record.put("vcf_header_lines", headerChunk);
        }
        record.put("vcf_header_lines_hash", headerHash);
        record.put("is_expected_unique", isExpectedUnique);
        return record;
    }

    public void commitData() {
        if (vcfHeaderBQJsonWriter != null) {
            vcfHeaderBQJsonWriter.flushBuffer();
            vcfHeaderBQJsonWriter.commitWriteStreams();
        }
    }

    public void closeTool() {
        if (vcfHeaderBQJsonWriter != null) {
            try {
                vcfHeaderBQJsonWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VCF Header Line writer", e);
            }
        }
        if (vcfHeaderBQJsonWriter != null) {
            vcfHeaderBQJsonWriter.close();
        }
    }

}
