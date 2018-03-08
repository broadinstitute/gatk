package org.broadinstitute.hellbender.tools.funcotator;

import java.util.*;

/**
 * Class to store argument definitions specific to {@link Funcotator}.
 * Created by jonn on 11/17/17.
 */
public class FuncotatorArgumentDefinitions {

    public static final char MAP_NAME_VALUE_DELIMITER = ':';

    // ------------------------------------------------------------
    // Definitions for required arguments:

    public static final String REFERENCE_VERSION_LONG_NAME = "ref-version";

    public static final String DATA_SOURCES_PATH_LONG_NAME = "data-sources-path";

    public static final String OUTPUT_FORMAT_LONG_NAME = "output-file-format";

    // ------------------------------------------------------------
    // Definitions for optional arguments:

    public static final String IGNORE_FILTERED_VARIANTS_LONG_NAME = "ignore-filtered-variants";

    public static final String TRANSCRIPT_SELECTION_MODE_LONG_NAME = "transcript-selection-mode";
    public static final TranscriptSelectionMode TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE = TranscriptSelectionMode.CANONICAL;

    public static final String TRANSCRIPT_LIST_LONG_NAME = "transcript-list";
    public static final Set<String> TRANSCRIPT_LIST_DEFAULT_VALUE = new HashSet<>();

    public static final String ANNOTATION_DEFAULTS_LONG_NAME = "annotation-default";
    public static final List<String> ANNOTATION_DEFAULTS_DEFAULT_VALUE = new ArrayList<>();

    public static final String ANNOTATION_OVERRIDES_LONG_NAME = "annotation-override";
    public static final List<String> ANNOTATION_OVERRIDES_DEFAULT_VALUE = new ArrayList<>();

    public static final String ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_LONG_NAME = "allow-hg19-gencode-b37-contig-matching";

    public static final String HG19_REFERENCE_VERSION_STRING = "hg19";
    public static final String HG38_REFERENCE_VERSION_STRING = "hg38";

    // ------------------------------------------------------------
    // Helper Types:

    /**
     * An enum to handle the different types of input files for data sources.
     */
    public enum DataSourceType {
        /**
         * This values indicates a simple arbitrary separated value (XSV) file that can be
         * annotated on records via matching by gene name or transcript ID.
         */
        SIMPLE_XSV("simpleXSV"),

        /**
         * This values indicates a simple arbitrary separated value (XSV) file that can be
         * annotated on records via matching by gene location.
         */
        LOCATABLE_XSV("locatableXSV"),

        /**
         * This values indicates a GENCODE GTF data file.
         */
        GENCODE("gencode"),

        /**
         * This values indicates a pre-processed COSMIC database file.
         * For more information on the pre-processing steps see {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory}
         * and the Funcotator Scripts directory.
         */
        COSMIC("cosmic");

        private final String serialized;

        DataSourceType(final String serializedValue) {
            serialized = serializedValue;
        }

        @Override
        public String toString() {
            return serialized;
        }

        public static DataSourceType getEnum(final String s) {
            for( final DataSourceType val : values() ) {
                if(val.serialized.equalsIgnoreCase(s)) {
                    return val;
                }
            }
            throw new IllegalArgumentException("Unexpected value: " + s);
        }
    }

    /**
     * The file format of the output file.
     */
    public enum OutputFormatType {
        VCF,
        MAF
    }
}
