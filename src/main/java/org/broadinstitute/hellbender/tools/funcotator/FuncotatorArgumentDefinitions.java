package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;

import java.util.*;

/**
 * Class to store argument definitions specific to {@link Funcotator}.
 * Created by jonn on 11/17/17.
 */
public class FuncotatorArgumentDefinitions {

    public static final char MAP_NAME_VALUE_DELIMITER = ':';

    // ------------------------------------------------------------
    // Definitions for required arguments:

    public static final String GTF_FILE_ARG_LONG_NAME = "gtfFile";
    public static final String GTF_FILE_ARG_SHORT_NAME = "gtf";

    public static final String GENCODE_FASTA_ARG_NAME = "fasta";

    // ------------------------------------------------------------
    // Definitions for optional arguments:

    public static final String LOCATABLE_XSV_IN_ARG_LONG_NAME = "xsv-locatable-input";
    public static final String LOCATABLE_XSV_IN_ARG_SHORT_NAME = "xsv-loc-in";

    public static final String XSV_INPUT_ARG_LONG_NAME = "xsv-input";
    public static final String XSV_INPUT_ARG_SHORT_NAME = "xsv";
    public static final List<String> XSV_INPUT_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String XSV_VERSION_ARG_LONG_NAME = "xsv-version";
    public static final String XSV_VERSION_ARG_SHORT_NAME = "xsvv";
    public static final List<String> XSV_VERSION_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String XSV_DELIMITER_ARG_LONG_NAME = "xsv-delimiter";
    public static final String XSV_DELIMITER_ARG_SHORT_NAME = "xsvd";
    public static final List<String> XSV_DELIMITER_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String XSV_KEY_COLUMN_ARG_LONG_NAME = "xsv-key-column";
    public static final String XSV_KEY_COLUMN_ARG_SHORT_NAME = "xsv-key-col";
    public static final List<Integer> XSV_KEY_COLUMN_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String XSV_FILE_TYPE_ARG_LONG_NAME = "xsv-file-type";
    public static final String XSV_FILE_TYPE_ARG_SHORT_NAME = "xsv-type";
    public static final List<SimpleKeyXsvFuncotationFactory.XsvDataKeyType> XSV_FILE_TYPE_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String XSV_NAME_ARG_LONG_NAME = "xsv-name";
    public static final String XSV_NAME_ARG_SHORT_NAME = "xsv-name";
    public static final List<String> XSV_NAME_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String XSV_PERMISSIVE_COLS_ARG_LONG_NAME = "xsv-permit-columns";
    public static final String XSV_PERMISSIVE_COLS_ARG_SHORT_NAME = "xsvpc";
    public static final List<Boolean> XSV_PERMISSIVE_COLS_ARG_DEFAULT_VALUE = new ArrayList<>();

    public static final String TRANSCRIPT_SELECTION_MODE_LONG_NAME = "transcript-selection-mode";
    public static final String TRANSCRIPT_SELECTION_MODE_SHORT_NAME = "tm";
    public static final TranscriptSelectionMode TRANSCRIPT_SELECTION_MODE_DEFAULT_VALUE = TranscriptSelectionMode.CANONICAL;

    public static final String TRANSCRIPT_LIST_LONG_NAME = "transcript-list";
    public static final String TRANSCRIPT_LIST_SHORT_NAME = "tl";
    public static final Set<String> TRANSCRIPT_LIST_DEFAULT_VALUE = new HashSet<>();

    public static final String ANNOTATION_DEFAULTS_LONG_NAME = "annotation-default";
    public static final String ANNOTATION_DEFAULTS_SHORT_NAME = "d";
    public static final List<String> ANNOTATION_DEFAULTS_DEFAULT_VALUE = new ArrayList<>();

    public static final String ANNOTATION_OVERRIDES_LONG_NAME = "annotation-override";
    public static final String ANNOTATION_OVERRIDES_SHORT_NAME = "a";
    public static final List<String> ANNOTATION_OVERRIDES_DEFAULT_VALUE = new ArrayList<>();

    // ------------------------------------------------------------
    // Helper Types:

    /**
     * The manner to select a single transcript from a set of transcripts to report as the "best" or main transcript.
     * That is, the transcript that has the most information on it in the output.
     */
    public enum TranscriptSelectionMode {

        /**
         * BEST_EFFECT
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
}
