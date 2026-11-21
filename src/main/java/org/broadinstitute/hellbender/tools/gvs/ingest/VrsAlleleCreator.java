package org.broadinstitute.hellbender.tools.gvs.ingest;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.tools.gvs.ingest.bq.VrsAlleleBQWriter;
import org.broadinstitute.hellbender.tools.gvs.ingest.parquet.VrsAlleleParquetFileWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Orchestrates writing of {@link VrsAlleleRecord} rows to the {@code vrs_allele_scratch} output.
 *
 * <p>One {@code VrsAlleleCreator} is created per sample, mirroring the {@link VetCreator}
 * pattern. It receives the list of fully-computed {@link VrsAlleleRecord} objects from
 * {@link VetCreator#apply} and writes them, tagged with the sample's id, to the appropriate output format.
 *
 * <p>Each row carries {@code sample_id} so the Parquet ingest verification, which is keyed on
 * {@code (table, sample_id)}, covers this table like any other. The same {@code vrs_allele_id} is
 * written once per sample that carries the allele; {@code ProcessVrsAlleles} merges the distinct rows
 * into the {@code vrs_allele} lookup table and then empties the scratch table.
 */
public class VrsAlleleCreator {

    private VrsAlleleWriter vrsAlleleWriter = null;
    private final Long sampleId;

    public static final String SCRATCH_TABLE_NAME = "vrs_allele_scratch";
    private static final String PREFIX_SEPARATOR = "_";

    public VrsAlleleCreator(String inputVcfFileName,
                            Long sampleId,
                            File outputDirectory,
                            CommonCode.OutputType outputType,
                            String projectId,
                            String datasetName,
                            MessageType parquetSchema) {
        this.sampleId = sampleId;
        try {
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    // All samples write to a single scratch table (no sharding needed)
                    vrsAlleleWriter = new VrsAlleleBQWriter(projectId, datasetName, SCRATCH_TABLE_NAME);
                    break;
                case PARQUET:
                    final String filename = getOutputFileName(sampleId, inputVcfFileName, outputType);
                    final File parquetOutputFile = new File(outputDirectory, filename);
                    vrsAlleleWriter = new VrsAlleleParquetFileWriter(
                            new Path(parquetOutputFile.toURI()), parquetSchema, CompressionCodecName.SNAPPY);
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
            vrsAlleleWriter.write(sampleId, record);
        }
    }

    /**
     * Named {@code vrs_allele_scratch_<sampleId>_<inputVcfFileName>} so the Parquet loader groups it
     * under the {@code vrs_allele_scratch} table and parses the sample id (regular-table pattern
     * {@code /<prefix>_<sampleId>...parquet}).
     */
    public static String getOutputFileName(Long sampleId,
                                           String inputVcfFileName,
                                           CommonCode.OutputType outputType) {
        String[] components = {SCRATCH_TABLE_NAME, sampleId.toString(), inputVcfFileName};
        return String.join(PREFIX_SEPARATOR, components) + "." + outputType.toString().toLowerCase();
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
