package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.apache.commons.lang.StringUtils;

import java.util.ArrayList;
import java.util.List;

public final class AnnotationUtils {
    private AnnotationUtils(){}

    /**
     * Helper function to parse the list into the annotation string
     * @param valueList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    protected static String encodeValueList(final List<Double> valueList, final String precisionFormat ) {
        List<String> outputList = new ArrayList<>();
        for (Double d : valueList) {
            outputList.add(String.format(precisionFormat, d));
        }
        return StringUtils.join(outputList, ",");
    }
}
