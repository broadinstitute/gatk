package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;
import java.util.*;

public final class AnnotationUtils {
    public static final String ALLELE_SPECIFIC_PRINT_DELIM = "|";
    public static final String ALLELE_SPECIFIC_REDUCED_DELIM = ",";
    public static final String ALLELE_SPECIFIC_SPLIT_REGEX = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String BRACKET_REGEX = "\\[|\\]";
    public static final String LIST_DELIMITER = ",";

    private AnnotationUtils(){}

    /**
     * Helper function to parse the list into the annotation string
     * @param valueList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    public static String encodeValueList(final List<Double> valueList, final String precisionFormat ) {
        List<String> outputList = new ArrayList<>();
        for (Double d : valueList) {
            outputList.add(String.format(precisionFormat, d));
        }
        return StringUtils.join(outputList, LIST_DELIMITER);
    }

    /**
     * Helper function to convert a List of Strings to a comma-separated String
     * @param stringList the ArrayList with String data
     * @return a comma-separated String
     */
    public static String encodeStringList( final List<String> stringList) {
        return StringUtils.join(stringList, LIST_DELIMITER);
    }

    /**
     * Helper function to convert a List of Strings to a pipe-separated String, as for raw annotations
     * @param somethingList the ArrayList with String data
     * @return a pipe-separated String
     */
    public static String encodeAnyASList( final List<?> somethingList) {
        return StringUtils.join(somethingList, ALLELE_SPECIFIC_PRINT_DELIM).replaceAll(BRACKET_REGEX, "");  //Who actually wants brackets at the ends of their string?  Who???
    }

    /**
     * Helper function to convert a comma-separated String (such as a vc.getAttrbute().toString() output) to a List of Strings
     * @param somethingList the allele-specific annotations string; may have brackets
     * @return a list of allele-specific annotation entries
     */
    public static List<String> decodeAnyASList( final String somethingList) {
        return Arrays.asList(StringUtils.splitByWholeSeparatorPreserveAllTokens(somethingList.replaceAll(BRACKET_REGEX, ""), ALLELE_SPECIFIC_PRINT_DELIM));
    }

    /**
     * Helper function to determine if an annotation is allele-specific
     * @param annotation the annotation to be tested
     * @return true if the annotation is expected to have values per-allele
     */
    public static boolean isAlleleSpecific(final InfoFieldAnnotation annotation) {
        return annotation instanceof AlleleSpecificAnnotation;
    }

    /**
     * Handles all the Java and htsjdk parsing shenanigans
     * @param rawDataString should not have surrounding brackets
     * @return
     */
    public static List<String> getAlleleLengthListOfString(String rawDataString) {
        if (rawDataString == null) {
            return Collections.emptyList();
        }
        if (rawDataString.startsWith("[")) {
            rawDataString = rawDataString.substring(1, rawDataString.length() - 1).replaceAll("\\s", "");
        }
        return Arrays.asList(rawDataString.split(ALLELE_SPECIFIC_SPLIT_REGEX));
    }
}
