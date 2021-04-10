package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class SampleInfoTsvCreator {

    private SimpleXSVWriter sampleInfoWriter = null;

    /**
     * Expected headers for the Sample List Table
     */
    public enum HeaderFieldEnum {
        sample_name,
        sample_id,
        inferred_state,
    }

    public SampleInfoTsvCreator(String sampleName, String sampleId, String tableNumberPrefix, final File outputDirectory) {
        try {
            final File sampleInfoFile = new File(outputDirectory, IngestConstants.sampleInfoFilePrefix + tableNumberPrefix + sampleName + IngestConstants.FILETYPE);
            // write header to it
            List<String> sampleListHeader = SampleInfoTsvCreator.getHeaders();
            sampleInfoWriter = new SimpleXSVWriter(sampleInfoFile.toPath(), IngestConstants.SEPARATOR);
            sampleInfoWriter.setHeaderLine(sampleListHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create sample info outputs", e);
        }

    }
    private List<String> createSampleListRow(
            String sampleName,
            String sampleId,
            PetTsvCreator.GQStateEnum inferredMissingState
    ) {

        List<String> row = new ArrayList<>();
        row.add(sampleName);
        row.add(sampleId);
        if (inferredMissingState == null) {
            row.add("");
        } else {
            row.add(inferredMissingState.getValue());
        }

        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(SampleInfoTsvCreator.HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void createRow(String sampleName, String sampleId, List<SimpleInterval> userIntervals, PetTsvCreator.GQStateEnum gqStateToIgnore) {
        // if the sample_info tsvs don't exist yet -- create them
        // Create a sample_info file to go into the sample_info dir for _this_ sample
        // TODO--this should just be one file per sample set?
        final List<String> TSVLineToCreateSampleInfo = createSampleListRow(
                sampleName,
                sampleId,
                gqStateToIgnore);
        sampleInfoWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleInfo).write();
    }

    public void closeTool() {
        if (sampleInfoWriter != null) {
            try {
                sampleInfoWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close SampleInfo writer", e);
            }
        }

    }
}
