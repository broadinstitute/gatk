package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.ingest.bq.VetBQWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.VetParquetFileWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.tsv.VetTsvWriter;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;

import java.io.File;
import java.io.IOException;

public class VetCreator {
    private VetWriter vetWriter = null;
    private final Long sampleId;

    private static final String VET_FILETYPE_PREFIX = "vet_";

    public static boolean doRowsExistFor(CommonCode.OutputType outputType, String projectId, String datasetName, String tableNumber, Long sampleId) {
        if (outputType != CommonCode.OutputType.BQ) return false;
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, VET_FILETYPE_PREFIX + tableNumber, SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
    }

    public VetCreator(String sampleIdentifierForOutputFileName, Long sampleId, String tableNumber, final File outputDirectory, final CommonCode.OutputType outputType, final String projectId, final String datasetName, final boolean forceLoadingFromNonAlleleSpecific, final boolean skipLoadingVqsrFields, final MessageType parquetSchema) {
        this.sampleId = sampleId;

        try {
            String PREFIX_SEPARATOR = "_";
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    vetWriter = new VetBQWriter(projectId, datasetName, VET_FILETYPE_PREFIX + tableNumber, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific);
                    break;
                case TSV:
                    // If the vet directory inside it doesn't exist yet -- create it
                    final File vetOutputFile = new File(outputDirectory, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + IngestConstants.FILETYPE);
                    vetWriter = new VetTsvWriter(vetOutputFile, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific);
                    break;
                case PARQUET:
                    String[] sampleComponents = {tableNumber, sampleId.toString(), sampleIdentifierForOutputFileName};
                    String filename = VET_FILETYPE_PREFIX + String.join(PREFIX_SEPARATOR, sampleComponents) +
                            "." + outputType.toString().toLowerCase();
                    final File parquetOutputFile = new File(outputDirectory, filename);
                    vetWriter = new VetParquetFileWriter(new Path(parquetOutputFile.toURI()), skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific, parquetSchema, CompressionCodecName.SNAPPY);
                    break;
            }
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void apply(VariantContext variant) throws IOException {
        final int start = variant.getStart();
        final long location = SchemaUtils.encodeLocation(variant.getContig(), start);

        vetWriter.write(location, variant, sampleId);
    }


    public void commitData() {
        if (vetWriter != null) {
            vetWriter.commitData();
        }
    }

    public void closeTool() {
        try {
            if (vetWriter != null) vetWriter.close();
        } catch (final Exception e) {
            throw new IllegalArgumentException(e);
        }
    }
}
