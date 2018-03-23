package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;

import java.nio.file.Path;
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
         * This value indicates a simple arbitrary separated value (XSV) file that can be
         * annotated on records via matching by gene name or transcript ID.
         */
        SIMPLE_XSV("simpleXSV") {
            @Override
            public void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_DELIMITER, configFileProperties, configFilePath);
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_KEY, configFileProperties, configFilePath);
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_KEY_COLUMN, configFileProperties, configFilePath);
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_PERMISSIVE_COLS, configFileProperties, configFilePath);

                // Ensure typed values:
                DataSourceUtils.assertIntegerPropertiesField(configFileProperties, DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_KEY_COLUMN, configFilePath);
                DataSourceUtils.assertBooleanPropertiesField(configFileProperties, DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_PERMISSIVE_COLS, configFilePath);

                // Validate our xsv_key:
                final String stringXsvKey = configFileProperties.getProperty(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_KEY);
                try {
                    SimpleKeyXsvFuncotationFactory.XsvDataKeyType.valueOf(stringXsvKey);
                }
                catch (final IllegalArgumentException ex) {
                    throw new UserException.BadInput("ERROR in config file: " + configFilePath.toUri().toString() +
                            " - Invalid value in \"" + DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_KEY + "\" field: " + stringXsvKey, ex);
                }
            }
        },

        /**
         * This value indicates a simple arbitrary separated value (XSV) file that can be
         * annotated on records via matching by gene location.
         */
        LOCATABLE_XSV("locatableXSV") {
            @Override
            public void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_XSV_DELIMITER, configFileProperties, configFilePath);
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_CONTIG_COLUMN, configFileProperties, configFilePath);
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_START_COLUMN, configFileProperties, configFilePath);
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_END_COLUMN, configFileProperties, configFilePath);

                // Ensure typed values:
                DataSourceUtils.assertIntegerPropertiesField(configFileProperties, DataSourceUtils.CONFIG_FILE_FIELD_NAME_CONTIG_COLUMN, configFilePath);
                DataSourceUtils.assertIntegerPropertiesField(configFileProperties, DataSourceUtils.CONFIG_FILE_FIELD_NAME_START_COLUMN, configFilePath);
                DataSourceUtils.assertIntegerPropertiesField(configFileProperties, DataSourceUtils.CONFIG_FILE_FIELD_NAME_END_COLUMN, configFilePath);
            }
        },

        /**
         * This value indicates a Variant Context File (VCF)
         */
        VCF("vcf") {
            @Override
            public void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
                // There is no special check required for vcf files.
            }
        },

        /**
         * This value indicates a GENCODE GTF data file.
         */
        GENCODE("gencode") {
            @Override
            public void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
                DataSourceUtils.assertConfigPropertiesContainsKey(DataSourceUtils.CONFIG_FILE_FIELD_NAME_GENCODE_FASTA_PATH, configFileProperties, configFilePath);

                // Assert that the path is good:
                DataSourceUtils.assertPathFilePropertiesField(configFileProperties, DataSourceUtils.CONFIG_FILE_FIELD_NAME_GENCODE_FASTA_PATH, configFilePath);
            }
        },

        /**
         * This value indicates a pre-processed COSMIC database file.
         * For more information on the pre-processing steps see {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory}
         * and the Funcotator Scripts directory.
         */
        COSMIC("cosmic") {
            @Override
            public void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath) {
                // There is no special check required for cosmic.
            }
        };

        /**
         * Asserts that the given properties and corresponding config file path are valid for this {@link DataSourceType}.
         * @param configFileProperties {@link Properties} pulled from the given {@code configFilePath}.
         * @param configFilePath {@link Path} to the config file containing the {@code configFileProperties}.
         */
        abstract public void assertConfigFilePropertiesAreValid(final Properties configFileProperties, final Path configFilePath);

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
