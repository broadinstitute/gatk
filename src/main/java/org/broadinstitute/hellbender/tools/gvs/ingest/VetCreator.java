package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

import org.json.JSONObject;

public class VetCreator {

    static final Logger logger = LogManager.getLogger(VetCreator.class);

    private CommonCode.OutputType outputType;

    private SimpleXSVWriter vetWriter = null;
    private final String sampleId;
    private static String VET_FILETYPE_PREFIX = "vet_";
    private static String PREFIX_SEPARATOR = "_";
    private PendingBQWriter vetBQJsonWriter = null;


    public VetCreator(String sampleIdentifierForOutputFileName, String sampleId, String tableNumber, final File outputDirectory, final CommonCode.OutputType outputType, final String projectId, final String datasetName) {
        this.sampleId = sampleId;
        this.outputType = outputType;
        try {
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    vetBQJsonWriter = new PendingBQWriter(projectId, datasetName,VET_FILETYPE_PREFIX + tableNumber);

                    break;
                case TSV:
                    // If the vet directory inside it doesn't exist yet -- create it
                    final File vetOutputFile = new File(outputDirectory, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + IngestConstants.FILETYPE);
                    vetWriter = new SimpleXSVWriter(vetOutputFile.toPath(), IngestConstants.SEPARATOR);
                    vetWriter.setHeaderLine(getHeaders());
            }
        } catch (final IOException ioex) {
            throw new UserException("Could not create vet outputs", ioex);
        }
    }

    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) throws IOException {
        int start = variant.getStart();
        long location = SchemaUtils.encodeLocation(variant.getContig(), start);
        List<String> row = createRow(
                location,
                variant,
                sampleId
        );

        switch(outputType) {
            case BQ:
                try {
                    vetBQJsonWriter.addJsonRow(createJson(location, variant, Long.parseLong(sampleId)));
                } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                    throw new IOException("BQ exception", ex);
                }
                break;
            case TSV:
                // write the variant to the XSV
                vetWriter.getNewLineBuilder().setRow(row).write();
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
            } else {
                row.add(fieldEnum.getColumnValue(variant));
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
                jsonObject.put(fieldEnum.toString(), Integer.valueOf(fieldEnum.getColumnValue(variant)));
            } else {
                jsonObject.put(fieldEnum.toString(), fieldEnum.getColumnValue(variant));
            }
        }
        return jsonObject;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(VetFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void commitData() {
        vetBQJsonWriter.flushBuffer();
        if (outputType == CommonCode.OutputType.BQ && vetBQJsonWriter != null) {
            vetBQJsonWriter.commitWriteStreams();
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

    }
}
