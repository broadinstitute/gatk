package org.broadinstitute.hellbender.tools.variantdb.arrays;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.nextgen.PetTsvCreator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class ArrayMetadataTsvCreator {

    private SimpleXSVWriter sampleMetadataWriter = null;

    /**
     * Expected headers for the Sample List Table
     */
    public enum HeaderFieldEnum {
        sample_name,
        sample_id
    }

    public static List<String> getHeaders() {
        return Arrays.stream(ArrayMetadataTsvCreator.HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void createRow(String sampleName, String sampleId, String tableNumberPrefix, File outputDirectory) {
        // if the metadata tsvs don't exist yet -- create them
        try {
            // Create a metadata file to go into the metadata dir for _this_ sample
            // TODO--this should just be one file per sample set?
            final File sampleMetadataName = new File (outputDirectory, IngestConstants.metadataFilePrefix + tableNumberPrefix + sampleName + IngestConstants.FILETYPE);
            // write header to it
            List<String> sampleListHeader = ArrayMetadataTsvCreator.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(sampleMetadataName.toPath(), IngestConstants.SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);

            final List<String> TSVLineToCreateSampleMetadata = createSampleListRow(
                    sampleName,
                    sampleId);
            sampleMetadataWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleMetadata).write();

        } catch (final IOException e) {
            throw new UserException("Could not create sample metadata outputs", e);
        }

    }

    private List<String> createSampleListRow(String sampleName, String sampleId) {
        List<String> row = new ArrayList<>();
        row.add(sampleName);
        row.add(sampleId);
        return row;
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
