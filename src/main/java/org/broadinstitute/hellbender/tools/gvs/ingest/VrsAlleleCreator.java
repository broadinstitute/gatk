package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.tools.gvs.ingest.bq.VrsAlleleBQWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.VrsAlleleParquetFileWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.tsv.VrsAlleleTsvWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Orchestrates writing of {@link VrsAlleleRecord} rows to the {@code vrs_allele} output.
 *
 * <p>One {@code VrsAlleleCreator} is created per sample, mirroring the {@link VetCreator}
 * pattern. It receives the list of fully-computed {@link VrsAlleleRecord} objects from
 * {@link VetCreator#apply} and writes them to the appropriate output format.
 *
 * <p>Duplicate rows (same {@code vrs_allele_id}) may be written across samples; deduplication
 * is deferred to query time via a BigQuery view or a post-load MERGE.
 */
public class VrsAlleleCreator {

    private VrsAlleleWriter vrsAlleleWriter = null;

    private static final String VRS_ALLELE_FILETYPE_PREFIX = "vrs_allele_";
    private static final String PREFIX_SEPARATOR = "_";

    public VrsAlleleCreator(String sampleIdentifierForOutputFileName,
                            Long sampleId,
                            String tableNumber,
                            File outputDirectory,
                            CommonCode.OutputType outputType,
                            String projectId,
                            String datasetName,
                            MessageType parquetSchema) {
        try {
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    // All samples write to a single vrs_allele table (no sharding needed)
                    vrsAlleleWriter = new VrsAlleleBQWriter(projectId, datasetName, "vrs_allele");
                    break;
                case TSV:
                    final File tsvOutputFile = new File(outputDirectory,
                            VRS_ALLELE_FILETYPE_PREFIX + tableNumber
                            + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName
                            + IngestConstants.FILETYPE);
                    vrsAlleleWriter = new VrsAlleleTsvWriter(tsvOutputFile);
                    break;
                case PARQUET:
                    final String filename = getOutputFileName(tableNumber, sampleId,
                            sampleIdentifierForOutputFileName, outputType);
                    final File parquetOutputFile = new File(outputDirectory, filename);
                    vrsAlleleWriter = new VrsAlleleParquetFileWriter(
                            new Path(parquetOutputFile.toURI()), parquetSchema, CompressionCodecName.SNAPPY);
                    break;
                default:
                    // Other output types not supported; writer stays null
                    break;
            }
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Write all VRS allele records produced for a single variant call.
     *
     * @param records list of records from {@link VetCreator#computeVrsAlleleRecords}
     * @throws IOException if writing fails
     */
    public void apply(List<VrsAlleleRecord> records) throws IOException {
        if (vrsAlleleWriter == null || records == null || records.isEmpty()) {
            return;
        }
        for (VrsAlleleRecord record : records) {
            vrsAlleleWriter.write(record);
        }
    }

    public static String getOutputFileName(String tableNumber, Long sampleId,
                                           String sampleIdentifierForOutputFileName,
                                           CommonCode.OutputType outputType) {
        String[] components = {tableNumber, sampleId.toString(), sampleIdentifierForOutputFileName};
        return VRS_ALLELE_FILETYPE_PREFIX + String.join(PREFIX_SEPARATOR, components)
                + "." + outputType.toString().toLowerCase();
    }

    public void commitData() {
        if (vrsAlleleWriter != null) {
            vrsAlleleWriter.commitData();
        }
    }

    public void closeTool() {
        try {
            if (vrsAlleleWriter != null) {
                vrsAlleleWriter.close();
            }
        } catch (final Exception e) {
            throw new IllegalArgumentException(e);
        }
    }
}
