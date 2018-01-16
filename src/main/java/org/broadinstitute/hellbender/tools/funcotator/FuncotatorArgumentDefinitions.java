package org.broadinstitute.hellbender.tools.funcotator;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
    public static final boolean IGNORE_FILTERED_VARIANTS_DEFFAULT_VALUE = false;

    public static final String TRANSCRIPT_SELECTION_MODE_LONG_NAME = "transcript-selection-mode";
    public static final TranscriptSelectionMode TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE = TranscriptSelectionMode.CANONICAL;

    public static final String TRANSCRIPT_LIST_LONG_NAME = "transcript-list";
    public static final Set<String> TRANSCRIPT_LIST_DEFAULT_VALUE = new HashSet<>();

    public static final String ANNOTATION_DEFAULTS_LONG_NAME = "annotation-default";
    public static final List<String> ANNOTATION_DEFAULTS_DEFAULT_VALUE = new ArrayList<>();

    public static final String ANNOTATION_OVERRIDES_LONG_NAME = "annotation-override";
    public static final List<String> ANNOTATION_OVERRIDES_DEFAULT_VALUE = new ArrayList<>();

    public static final String ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_LONG_NAME = "allow-hg19-gencode-b37-contig-matching";
    public static final boolean ALLOW_HG19_GENCODE_B37_CONTIG_MATCHING_DEFAULT_VALUE = false;

    public static final String HG19_REFERENCE_VERSION_STRING = "hg19";
    public static final String HG38_REFERENCE_VERSION_STRING = "hg38";

    // ------------------------------------------------------------
    // Helper Types:

    /**
     * The manner to select a single transcript from a set of transcripts to report as the "best" or main transcript.
     * That is, the transcript that has the most information on it in the output.
     */
    public enum TranscriptSelectionMode {

        /**
         * BEST_EFFECT
         *
         * Select a transcript to be reported with details with priority on effect according to the folowing list of selection criteria:
         *
         * Choose the transcript that is on the custom list specified by the user. If no list was specified, treat as if no transcripts were on the list (tie).
         * In case of tie, choose the transcript that yields the variant classification highest on the variant classification rank list (see below).
         * If still a tie, choose the transcript with highest level of curation. Note that this means lower number is better for level (see below).
         * If still a tie, choose the transcript with the best appris annotation (see below).
         * If still a tie, choose the transcript with the longest transcript sequence length.
         * If still a tie, choose the first transcript, alphabetically.
         *
         * Levels of Curation:
         *     1 (verified loci)
         *     2 manually annotated loci
         *     3 automatically annotated loci
         *
         * Variant Classification Scores (See {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} as well):
         *      De_novo_Start_OutOfFrame    0
         *      Nonsense_Mutation           0
         *      Nonstop_Mutation            0
         *      Missense_Mutation           1
         *      De_novo_Start_InFrame       1
         *      In_Frame_Del                1
         *      In_Frame_Ins                1
         *      Frame_Shift_Del             2
         *      Frame_Shift_Ins             2
         *      Frame_Shift_Sub             2
         *      Start_Codon_SNP             3
         *      Start_Codon_Del             3
         *      Start_Codon_Ins             3
         *      Start_Codon_DNP             3
         *      Start_Codon_TNP             3
         *      Start_Codon_ONP             3
         *      Stop_Codon_SNP              3
         *      Stop_Codon_Del              3
         *      Stop_Codon_Ins              3
         *      Stop_Codon_DNP              3
         *      Stop_Codon_TNP              3
         *      Stop_Codon_ONP              3
         *      Splice_Site                 4
         *      Splice_Site_SNP             4
         *      Splice_Site_Del             4
         *      Splice_Site_Ins             4
         *      Splice_Site_DNP             4
         *      Splice_Site_TNP             4
         *      Splice_Site_ONP             4
         *      Splice_Site                 4
         *      miRNA                       4
         *      Silent                      5
         *      3UTR                        6
         *      5UTR                        6
         *      Intron                      7
         *      5Flank                      8
         *      3Flank                      8
         *      Non-coding_Transcript       9
         *      IGR	                        20
         *      TX-REF-MISMATCH             100
         *
         * APPRIS Ranks (http://appris.bioinfo.cnio.es/):
         *
         *      appris_principal
         *      appris_candidate_highest_score
         *      appris_candidate_longest_ccds
         *      appris_candidate_ccds
         *      appris_candidate_longest_seq
         *      appris_candidate_longest
         *      appris_candidate
         *      no appris tag present
         */
        BEST_EFFECT,

        /**
         * CANONICAL
         *
         * Select a transcript to be reported with details with priority on canonical order according to the folowing list of selection criteria:
         *
         * Choose the transcript that is on the custom list specified by the user. If no list was specified, treat as if all transcripts were on the list (tie).
         * In case of tie, choose the transcript with highest level of curation. Note that this means lower number is better for level (see below).
         * If still a tie, choose the transcript that yields the variant classification highest on the variant classification rank list (see below).
         * If still a tie, choose the transcript with the best appris annotation (see below).
         * If still a tie, choose the transcript with the longest transcript sequence length.
         * If still a tie, choose the first transcript, alphabetically.
         *
         * Levels of Curation:
         *     1 (verified loci)
         *     2 manually annotated loci
         *     3 automatically annotated loci
         *
         * Variant Classification Scores (See {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation.VariantClassification} as well):
         *      De_novo_Start_OutOfFrame    0
         *      Nonsense_Mutation           0
         *      Nonstop_Mutation            0
         *      Missense_Mutation           1
         *      De_novo_Start_InFrame       1
         *      In_Frame_Del                1
         *      In_Frame_Ins                1
         *      Frame_Shift_Del             2
         *      Frame_Shift_Ins             2
         *      Frame_Shift_Sub             2
         *      Start_Codon_SNP             3
         *      Start_Codon_Del             3
         *      Start_Codon_Ins             3
         *      Start_Codon_DNP             3
         *      Start_Codon_TNP             3
         *      Start_Codon_ONP             3
         *      Stop_Codon_SNP              3
         *      Stop_Codon_Del              3
         *      Stop_Codon_Ins              3
         *      Stop_Codon_DNP              3
         *      Stop_Codon_TNP              3
         *      Stop_Codon_ONP              3
         *      Splice_Site                 4
         *      Splice_Site_SNP             4
         *      Splice_Site_Del             4
         *      Splice_Site_Ins             4
         *      Splice_Site_DNP             4
         *      Splice_Site_TNP             4
         *      Splice_Site_ONP             4
         *      Splice_Site                 4
         *      miRNA                       4
         *      Silent                      5
         *      3UTR                        6
         *      5UTR                        6
         *      Intron                      7
         *      5Flank                      8
         *      3Flank                      8
         *      Non-coding_Transcript       9
         *      IGR	                        20
         *      TX-REF-MISMATCH             100
         *
         * APPRIS Ranks (http://appris.bioinfo.cnio.es/):
         *
         *      appris_principal
         *      appris_candidate_highest_score
         *      appris_candidate_longest_ccds
         *      appris_candidate_ccds
         *      appris_candidate_longest_seq
         *      appris_candidate_longest
         *      appris_candidate
         *      no appris tag present
         */
        CANONICAL
    }

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
