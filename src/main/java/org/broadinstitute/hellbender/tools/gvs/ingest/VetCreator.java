package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.fs.FileAlreadyExistsException;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.PendingBQWriter;
import org.broadinstitute.hellbender.utils.gvs.parquet.GvsVariantParquetFileWriter;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

public class VetCreator {
    private final CommonCode.OutputType outputType;

    private SimpleXSVWriter vetWriter = null;
    private final Long sampleId;
    private PendingBQWriter vetBQJsonWriter = null;
    private GvsVariantParquetFileWriter vetParquetFileWriter = null;
    private final boolean forceLoadingFromNonAlleleSpecific;
    private final boolean skipLoadingVqsrFields;

    private static final String VET_FILETYPE_PREFIX = "vet_";

    public static boolean doRowsExistFor(CommonCode.OutputType outputType, String projectId, String datasetName, String tableNumber, Long sampleId) {
        if (outputType != CommonCode.OutputType.BQ) return false;
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, VET_FILETYPE_PREFIX + tableNumber, SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
    }

    public VetCreator(String sampleIdentifierForOutputFileName, Long sampleId, String tableNumber, final File outputDirectory, final CommonCode.OutputType outputType, final String projectId, final String datasetName, final boolean forceLoadingFromNonAlleleSpecific, final boolean skipLoadingVqsrFields, final MessageType parquetSchema) {
        this.sampleId = sampleId;
        this.outputType = outputType;
        this.forceLoadingFromNonAlleleSpecific = forceLoadingFromNonAlleleSpecific;
        this.skipLoadingVqsrFields = skipLoadingVqsrFields;

        try {
            String PREFIX_SEPARATOR = "_";
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    vetBQJsonWriter = new PendingBQWriter(projectId, datasetName, VET_FILETYPE_PREFIX + tableNumber);

                    break;
                case TSV:
                    // If the vet directory inside it doesn't exist yet -- create it
                    final File vetOutputFile = new File(outputDirectory, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + IngestConstants.FILETYPE);
                    vetWriter = new SimpleXSVWriter(vetOutputFile.toPath(), IngestConstants.SEPARATOR);
                    vetWriter.setHeaderLine(getHeaders());
                    break;
                case PARQUET:
                    final File parquetOutputFile = new File(outputDirectory, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + ".parquet");
                    vetParquetFileWriter = new GvsVariantParquetFileWriter(new Path(parquetOutputFile.toURI()), parquetSchema, false, CompressionCodecName.SNAPPY);
                    break;
            }
        } catch (final FileAlreadyExistsException fs) {
            throw new UserException("This variants parquet file already exists", fs);
        } catch (final IOException ioex) {
            throw new UserException("Could not create vet outputs", ioex);
        }
    }

    public void apply(VariantContext variant) throws IOException {
        final int start = variant.getStart();
        final long location = SchemaUtils.encodeLocation(variant.getContig(), start);

        switch(outputType) {
            case BQ:
                try {
                    vetBQJsonWriter.addJsonRow(createJson(location, variant, sampleId));
                } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                    throw new IOException("BQ exception", ex);
                }
                break;
            case TSV:
                // write the variant to the XSV
                final List<String> row = createRow(
                        location,
                        variant,
                        String.valueOf(sampleId)
                );
                vetWriter.getNewLineBuilder().setRow(row).write();
                break;
            case PARQUET:
                vetParquetFileWriter.write(createJson(location, variant, sampleId));
                break;

        }
    }

   public List<String> createRow(final long location, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        for ( final VetFieldEnum fieldEnum : VetFieldEnum.values() ) {
            if (fieldEnum.equals(VetFieldEnum.location)) {
                row.add(String.valueOf(location));
            } else if (fieldEnum.equals(VetFieldEnum.sample_id)) {
                row.add(sampleId);
            } else if (!(skipLoadingVqsrFields && fieldEnum.isVqsrSpecificField())) {
                row.add(fieldEnum.getColumnValue(variant, forceLoadingFromNonAlleleSpecific));
            }
        }
        return row;
    }

    // Similar to create row but with types for JSON object
    public JSONObject createJson(final long location, final VariantContext variant, final long sampleId) {
        JSONObject jsonObject = new JSONObject();
        for ( final VetFieldEnum fieldEnum : VetFieldEnum.values() ) {
            if (fieldEnum.equals(VetFieldEnum.location)) {
                jsonObject.put(VetFieldEnum.location.toString(), location);
            } else if (fieldEnum.equals(VetFieldEnum.sample_id)) {
                jsonObject.put(VetFieldEnum.sample_id.toString(), sampleId);
            } else if (fieldEnum.equals(VetFieldEnum.call_GQ)) {
                jsonObject.put(fieldEnum.toString(), Integer.valueOf(fieldEnum.getColumnValue(variant, forceLoadingFromNonAlleleSpecific)));
            } else {
                final String strVal = !(skipLoadingVqsrFields && fieldEnum.isVqsrSpecificField()) ? fieldEnum.getColumnValue(variant, forceLoadingFromNonAlleleSpecific) : "";
                jsonObject.put(fieldEnum.toString(), StringUtils.isEmpty(strVal) ? null : strVal);
            }
        }
        return jsonObject;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(VetFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void commitData() {
        if (outputType == CommonCode.OutputType.BQ && vetBQJsonWriter != null) {
            vetBQJsonWriter.flushBuffer();
            vetBQJsonWriter.commitWriteStreams();
        } else if (outputType == CommonCode.OutputType.PARQUET && vetParquetFileWriter != null) {
            try {
                vetParquetFileWriter.close();
            } catch (IOException exception) {
                System.out.println("ERROR CLOSING PARQUET FILE: ");
                exception.printStackTrace();
            }
        }
    }

    public void closeTool() {
        if (vetWriter != null) {
            try {
                vetWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }
        if (vetBQJsonWriter != null) {
            vetBQJsonWriter.close();
        }
    }
}
