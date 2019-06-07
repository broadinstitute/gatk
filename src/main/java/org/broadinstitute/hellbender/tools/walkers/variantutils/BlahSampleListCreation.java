package org.broadinstitute.hellbender.tools.walkers.variantutils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class BlahSampleListCreation {

    /**
     * Expected headers for the Position Table (PET)
     */
    public enum HeaderFieldEnum {
        sample,
        interval_list,
        inferred_state,
    }

    public static List<String> createSampleListRow(String sampleName, String intervalListPath, BlahPetCreation.GQStateEnum inferredMissingState) {

        List<String> row = new ArrayList<>();
        row.add(sampleName);
        row.add(intervalListPath);
        if (inferredMissingState == null) {
            row.add("");
        } else {
            row.add(inferredMissingState.value);
        }

        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(BlahSampleListCreation.HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
