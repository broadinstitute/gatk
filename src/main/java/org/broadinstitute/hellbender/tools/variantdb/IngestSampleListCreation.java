package org.broadinstitute.hellbender.tools.variantdb;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class IngestSampleListCreation {

    /**
     * Expected headers for the Position Table (PET)
     */
    public enum HeaderFieldEnum {
        sample,
        sample_id,
        interval_list_blob,
        inferred_state,
    }

    public static List<String> createSampleListRow(
            String sampleName,
            String sampleId,
            String intervalListBlob,
            IngestPetCreation.GQStateEnum inferredMissingState
    ) {

        List<String> row = new ArrayList<>();
        row.add(sampleName);
        row.add(sampleId);
        row.add(intervalListBlob);
        if (inferredMissingState == null) {
            row.add("");
        } else {
            row.add(inferredMissingState.value);
        }

        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(IngestSampleListCreation.HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
