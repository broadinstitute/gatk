package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class MetadataTsvCreator {

    private SimpleXSVWriter sampleMetadataWriter = null;

    /**
     * Expected headers for the Sample List Table
     */
    public enum HeaderFieldEnum {
        sample_name,
        sample_id,
        interval_list_blob,
        inferred_state,
    }

    private List<String> createSampleListRow(
            String sampleName,
            String sampleId,
            String intervalListBlob,
            PetTsvCreator.GQStateEnum inferredMissingState
    ) {

        List<String> row = new ArrayList<>();
        row.add(sampleName);
        row.add(sampleId);
        row.add(intervalListBlob);
        if (inferredMissingState == null) {
            row.add("");
        } else {
            row.add(inferredMissingState.getValue());
        }

        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(MetadataTsvCreator.HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

//    public void createRow(String sampleName, String sampleId) {
//        createRow(sampleName, sampleId, null, PetTsvCreator.GQStateEnum.SIXTY);
//
//    }
//
    public void createRow(String sampleName, String sampleId, String tableNumberPrefix, List<SimpleInterval> userIntervals, PetTsvCreator.GQStateEnum gqStateToIgnore) {
        // if the metadata tsvs don't exist yet -- create them
        try {
            // Create a metadata file to go into the metadata dir for _this_ sample
            // TODO--this should just be one file per sample set?
            final String sampleMetadataName = IngestConstants.sampleMetadataFilePrefix + tableNumberPrefix + sampleName + IngestConstants.FILETYPE;
            // write header to it
            List<String> sampleListHeader = MetadataTsvCreator.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(Paths.get(sampleMetadataName), IngestConstants.SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);
            String intervalListMd5 = "NA";

            if (userIntervals != null) {
                // write values
                List<String> intervalList = userIntervals.stream().map(interval -> interval.toString())
                        .collect(Collectors.toList());
                String intervalListBlob = StringUtils.join(intervalList, ", ");
                intervalListMd5 = Utils.calcMD5(intervalListBlob);
            }
            final List<String> TSVLineToCreateSampleMetadata = createSampleListRow(
                    sampleName,
                    sampleId,
                    intervalListMd5,
                    gqStateToIgnore);
            sampleMetadataWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleMetadata).write();

        } catch (final IOException e) {
            throw new UserException("Could not create sample metadata outputs", e);
        }

    }

    public void closeTool() {
        if (sampleMetadataWriter != null) {
            try {
                sampleMetadataWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }

    }
}
