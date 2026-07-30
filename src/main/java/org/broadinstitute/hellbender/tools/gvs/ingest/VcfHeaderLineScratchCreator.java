package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.apache.hadoop.fs.Path;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.HeaderParquetFileWriter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ExecutionException;


public class VcfHeaderLineScratchCreator {
    private static final Logger logger = LogManager.getLogger(VcfHeaderLineScratchCreator.class);

    /**
     * How the Parquet ingest path writes VCF header data to the scratch file. Most header text is a large
     * blob that is identical across all samples in a callset (only a few command-line chunks differ per
     * sample), so the strategies trade off write simplicity against how many redundant copies of that
     * shared blob get written. See the VS-1968 design doc for the full analysis.
     */
    public enum HeaderParquetStrategy {
        /**
         * Write the full header text for every chunk of every sample, redundant copies of the shared blob
         * included. Deduplication happens downstream in the scratch&rarr;final promotion, so the scratch
         * file is larger but the write path is trivial. This is the only implemented strategy.
         */
        NAIVE,

        /**
         * Deduplicate the shared blob at write time by routing on {@code is_expected_unique}: per-sample
         * unique (command-line) chunks are written inline, while the shared blob is content-addressed to
         * GCS keyed by its hash so it is stored exactly once for the whole callset. Cheaper at scale, but
         * the GCS content-addressing half is not yet implemented &mdash; selecting this today throws
         * (see {@link #apply}).
         */
        HYBRID
    }

    // TODO(VS-1803): promote this to a real CLI/WDL parameter. Hardcoded for now.
    private static final HeaderParquetStrategy PARQUET_STRATEGY = HeaderParquetStrategy.NAIVE;

    private final CommonCode.OutputType outputType;
    private final Long sampleId;
    private final String projectId;
    private final String datasetName;

    private PendingBQWriter vcfHeaderBQJsonWriter = null;
    private HeaderParquetFileWriter vcfHeaderParquetFileWriter = null;
    private static final String NON_SCRATCH_TABLE_NAME = "vcf_header_lines";
    private static final String SCRATCH_TABLE_NAME = "vcf_header_lines_scratch";


    public static boolean doScratchRowsExistFor(String projectId, String datasetName, Long sampleId) {
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, "vcf_header_lines_scratch", "sample_id", sampleId);
    }

    public static boolean doNonScratchRowsExistFor(String projectId, String datasetName, Long sampleId) {
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, "sample_vcf_header", "sample_id", sampleId);
    }

    private static boolean doScratchRowsExistFor(String projectId, String datasetName, String headerLineHash) {
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, SCRATCH_TABLE_NAME, "vcf_header_lines_hash", headerLineHash);
    }

    private static boolean doNonScratchRowsExistFor(String projectId, String datasetName, String headerLineHash) {
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, NON_SCRATCH_TABLE_NAME, "vcf_header_lines_hash", headerLineHash);
    }

    public VcfHeaderLineScratchCreator(Long sampleId, String projectId, String datasetName, File outputDirectory, CommonCode.OutputType outputType, MessageType headersRowSchema) {
        try {
            this.sampleId = sampleId;
            this.outputType = outputType;
            this.projectId = projectId;
            this.datasetName = datasetName;

            // String PREFIX_SEPARATOR = "_"; // TODO should this be moved to a common place?

            if (projectId == null || datasetName == null) {
                throw new UserException("Must specify project-id and dataset-name.");
            }

            switch (outputType) {

                case BQ:
                    // null-check for projectId/datasetName already performed above
                    vcfHeaderBQJsonWriter = new PendingBQWriter(projectId, datasetName, SCRATCH_TABLE_NAME);
                    break;
                case PARQUET:
                    // One header file per sample. Name it "<SCRATCH_TABLE_NAME>_<sampleId>.parquet" so the
                    // loader groups it under the vcf_header_lines_scratch table and parses the sample id
                    // (regular-table pattern: /<prefix>_<sampleId>...parquet). See the VS-1968 design doc.
                    final File parquetOutputFile = new File(outputDirectory, SCRATCH_TABLE_NAME + "_" + this.sampleId + ".parquet");
                    vcfHeaderParquetFileWriter = new HeaderParquetFileWriter(new Path(parquetOutputFile.toURI()), headersRowSchema, CompressionCodecName.SNAPPY);
                    break;

            }
        }
        catch (Exception e) {
            throw new UserException("Could not create VCF Header Scratch Table writer", e);
        }

    }

    public void apply(Map<String, Boolean>  allLineHeaders) throws IOException {
        for (final Map.Entry<String, Boolean> headerChunk : allLineHeaders.entrySet()) {
            switch (outputType) {
                case BQ:
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
                    break;
                case PARQUET: {
                    // The Parquet path does no write-time dedup (that would reintroduce the per-sample BQ
                    // query storm); dedup + idempotency happen in the scratch->final load. See the design doc.
                    final String chunkHash = Utils.calcMD5(headerChunk.getKey());
                    final Boolean isExpectedUnique = headerChunk.getValue();
                    switch (PARQUET_STRATEGY) {
                        case NAIVE:
                            // Write the full record for every chunk of every sample; dedup happens downstream.
                            vcfHeaderParquetFileWriter.write(
                                    HeaderParquetFileWriter.writeJson(this.sampleId, headerChunk.getKey(), chunkHash, isExpectedUnique));
                            break;
                        case HYBRID:
                            // The shared blob is content-addressed to GCS keyed by its hash (see the enum
                            // javadoc), and that half is not yet implemented -- so HYBRID cannot durably persist
                            // the shared blob text, and selecting it would silently drop header data.
                            // Fail loudly rather than write a partial, lossy scratch file.
                            // TODO(VS-1803): write the blob text to gs://.../headers/text/<hash>.parquet with an
                            //   ifGenerationMatch=0 precondition, load those objects into vcf_header_lines, then
                            //   restore the routing here and promote PARQUET_STRATEGY to a real CLI/WDL parameter.
                            throw new UnsupportedOperationException(
                                    "HYBRID header Parquet strategy is not yet implemented (shared-blob content-addressing " +
                                    "is incomplete; see TODO VS-1803). Use NAIVE until it lands.");
                    }
                    break;
                }
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
        if (outputType == CommonCode.OutputType.BQ && vcfHeaderBQJsonWriter != null) {
            vcfHeaderBQJsonWriter.flushBuffer();
            vcfHeaderBQJsonWriter.commitWriteStreams();
        } else if (outputType == CommonCode.OutputType.PARQUET && vcfHeaderParquetFileWriter != null) {
            try {
                vcfHeaderParquetFileWriter.close();
            } catch (IOException exception) {
                logger.error("Error closing Parquet VCF header file", exception);
            }
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
    }
}
