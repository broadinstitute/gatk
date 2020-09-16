package org.broadinstitute.hellbender.tools.variantdb.arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ArraySampleTsvCreator {

    private static final Logger logger = LogManager.getLogger(ArraySampleTsvCreator.class);

    private SimpleXSVWriter sampleMetadataWriter = null;
    private Map<String, String> metricsMap;


    public ArraySampleTsvCreator(String metricsFilepath) {
        if (metricsFilepath != null && !metricsFilepath.isEmpty()) {
            BufferedReader reader = null;
            try {
                String columns = null;
                String values = null;
                reader = new BufferedReader(new FileReader(metricsFilepath));
                String line = reader.readLine();
                while (line != null) {
                    if (!line.startsWith("#") && !line.trim().isEmpty()) {
                        if (columns == null) {
                            columns = line;
                        } else if (values == null) {
                            values = line;
                        } else {
                            // there are more lines than expected - output a warning
                            logger.warn("more lines than expected in metrics file: " + line);
                        }
                    }
                    line = reader.readLine();
                }

                List<String> colList = Arrays.asList(columns.split("\t"));
                List<String> valList = Arrays.asList(values.split("\t"));
                metricsMap = IntStream.range(0, colList.size()).boxed().collect(Collectors.toMap(colList::get, valList::get));

            } catch (IOException e) {
                throw new RuntimeException("could not read metrics file", e);
            }
        }
    }

    public static List<String> getHeaders() {
        return Arrays.stream(ArraySampleFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void createRow(String sampleName, String sampleId, String tableNumberPrefix, File outputDirectory) {
        // if the metadata tsvs don't exist yet -- create them
        try {
            // Create a metadata file to go into the metadata dir for _this_ sample
            // TODO--this should just be one file per sample set?
            final File sampleMetadataFileName = new File (outputDirectory, IngestConstants.sampleMetadataFilePrefix + tableNumberPrefix + sampleName + IngestConstants.FILETYPE);
            // write header to it
            List<String> sampleListHeader = ArraySampleTsvCreator.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(sampleMetadataFileName.toPath(), IngestConstants.SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);

            final List<String> TSVLineToCreateSampleMetadata = createSampleListRow(
                    sampleName,
                    sampleId);
            sampleMetadataWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleMetadata).write();

        } catch (final IOException e) {
            throw new UserException("Could not create sample outputs", e);
        }

    }

    private List<String> createSampleListRow(String sampleName, String sampleId) {
        List<String> row = new ArrayList<>();
        row.add(sampleName);
        row.add(sampleId);

        for (final ArraySampleFieldEnum fieldEnum : ArraySampleFieldEnum.values()) {
            if (fieldEnum != ArraySampleFieldEnum.sample_id && fieldEnum != ArraySampleFieldEnum.sample_name) {
                if (metricsMap == null) {
                    row.add("null");
                } else {
                    row.add(fieldEnum.getColumnValue(metricsMap));
                }
            }
        }
        return row;
    }

    public void closeTool() {
        if (sampleMetadataWriter != null) {
            try {
                sampleMetadataWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close array sample writer", e);
            }
        }

    }
}
