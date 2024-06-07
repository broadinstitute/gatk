package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.broadinstitute.hellbender.utils.gvs.parquet.GvsHeaderParquetFileWriter;
import org.broadinstitute.hellbender.utils.gvs.parquet.GvsReferenceParquetFileWriter;
import org.broadinstitute.hellbender.utils.gvs.parquet.GvsVariantParquetFileWriter;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;



public class VcfHeaderLineScratchCreator {
    private final CommonCode.OutputType outputType;
    private final Long sampleId;
    private final String projectId;
    private final String datasetName;

    private PendingBQWriter vcfHeaderBQJsonWriter = null;
    private GvsHeaderParquetFileWriter vcfHeaderParquetFileWriter = null;
    private static final String NON_SCRATCH_TABLE_NAME = "vcf_header_lines";
    private static final String SCRATCH_TABLE_NAME = "vcf_header_lines_scratch";

    private static final String HEADER_FILETYPE_PREFIX = "header_file";


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
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    vcfHeaderBQJsonWriter = new PendingBQWriter(projectId, datasetName, SCRATCH_TABLE_NAME);
                    break;
                case PARQUET:
                    // TODO ensure that there doesn't need to be a table_number or sampleIdentifierForOutputFileName--it's all tables/samples, yes?
                    final File parquetOutputFile = new File(outputDirectory, HEADER_FILETYPE_PREFIX + ".parquet");
                    vcfHeaderParquetFileWriter = new GvsHeaderParquetFileWriter(new Path(parquetOutputFile.toURI()), headersRowSchema, false, CompressionCodecName.SNAPPY);
                    break;

            }
        }
        catch (Exception e) {
            throw new UserException("Could not create VCF Header Scratch Table Writer", e);
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
                case PARQUET:
                    String chunkHash = Utils.calcMD5(headerChunk.getKey());
                    Boolean isExpectedUnique = headerChunk.getValue();
                    JSONObject record = vcfHeaderParquetFileWriter.writeJson(this.sampleId, chunkHash);
                    vcfHeaderParquetFileWriter.write(record);

                    //boolean vcfScratchHeaderRowsExist = doScratchRowsExistFor(this.projectId, this.datasetName, chunkHash);
                    //boolean vcfNonScratchHeaderRowsExist = doNonScratchRowsExistFor(this.projectId, this.datasetName, chunkHash);
                    // why is there no isExpectedUnique check here?
                    //if (vcfScratchHeaderRowsExist || vcfNonScratchHeaderRowsExist) {
                        //vcfHeaderParquetFileWriter.writeJson(this.sampleId, chunkHash);
                    //}
                    //else {
                        //vcfHeaderParquetFileWriter.writeJson(this.sampleId, chunkHash);
                    //}
                    break;
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
                System.out.println("ERROR CLOSING PARQUET FILE: ");
                exception.printStackTrace();
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
        if (vcfHeaderBQJsonWriter != null) {
            vcfHeaderBQJsonWriter.close();
        }
    }
}
