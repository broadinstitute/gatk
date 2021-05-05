package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.variant.variantcontext.Allele;

public final class GATKSVVCFConstants {

    // todo: add these and the other standard SV info fields from the VCF spec to htsjdk VCFStandardHeaderLines
    // VCF standard keys reserved for sv
    public static final String SVTYPE = "SVTYPE";
    public static final String SVLEN = "SVLEN";
    public static final String IMPRECISE = "IMPRECISE";
    public static final String CIPOS = "CIPOS";
    public static final String CIEND = "CIEND";

    public static final String BREAKEND_STR = "BND";
    public static final String BND_MATEID_STR = "MATEID";

    // symbolic alt allele names
    public static final String SYMB_ALT_STRING_INV = "INV";
    public static final String SYMB_ALT_STRING_DEL = "DEL";
    public static final String SYMB_ALT_STRING_INS = "INS";
    public static final String SYMB_ALT_STRING_DUP = "DUP";
    public static final String SYMB_ALT_STRING_INVDUP = "DUP:INV";

    // symbolic alt alleles
    public static final Allele DEL_ALLELE = Allele.create("<DEL>", false);
    public static final Allele DUP_ALLELE = Allele.create("<DUP>", false);

    // GATK-SV specific header lines
    // TODO: 10/3/17 the following comment is a goal we are trying to achieve
    // applicable to all records all the time
    public static final String READ_PAIR_SUPPORT = "READ_PAIR_SUPPORT";
    public static final String SPLIT_READ_SUPPORT = "SPLIT_READ_SUPPORT";

    // applicable to all precise records all the time
    public static final String CONTIG_NAMES = "CTG_NAMES";
    public static final String TOTAL_MAPPINGS = "TOTAL_MAPPINGS";
    public static final String MAPPING_QUALITIES = "MAPPING_QUALITIES";
    public static final String HQ_MAPPINGS = "HQ_MAPPINGS";
    public static final String ALIGN_LENGTHS = "ALIGN_LENGTHS";
    public static final String MAX_ALIGN_LENGTH = "MAX_ALIGN_LENGTH";

    // applicable to all precise records when available
    public static final String SEQ_ALT_HAPLOTYPE = "SEQ_ALT_HAPLOTYPE";
    public static final String INSERTED_SEQUENCE = "INSSEQ";
    public static final String INSERTED_SEQUENCE_LENGTH = "INSLEN";
    public static final String INSERTED_SEQUENCE_MAPPINGS = "INSSEQ_MAP";
    public static final String HOMOLOGY = "HOMSEQ";
    public static final String HOMOLOGY_LENGTH = "HOMLEN";
    public static final String LINK = "LINK";
    public static final String EXTERNAL_CNV_CALLS = "EXTERNAL_CNV_CALLS";
    public static final String CTG_GOOD_NONCANONICAL_MAPPING = "CTG_GOOD_NONCANONICAL_MAPPING";

    // type specific: tandem duplication
    public static final String DUP_REPEAT_UNIT_REF_SPAN = "DUP_REPEAT_UNIT_REF_SPAN";
    public static final String DUP_SEQ_CIGARS = "DUP_SEQ_CIGARS";
    public static final String DUPLICATION_NUMBERS = "DUP_NUM";
    public static final String DUP_ANNOTATIONS_IMPRECISE = "DUP_ANNOTATIONS_IMPRECISE";
    public static final String DUP_IMPRECISE_AFFECTED_RANGE = "DUP_IMPRECISE_AFFECTED_RANGE";

    public static final String DUP_TAN_CONTRACTION_STRING = "CONTRACTION";
    public static final String DUP_TAN_EXPANSION_STRING = "EXPANSION";

    // type specific: inverted duplication
    public static final String DUP_ORIENTATIONS = "DUP_ORIENTATIONS";

    // type specific: inversion
    public static final String INV33 = "INV33";
    public static final String INV55 = "INV55";

    // type specific: CPX
    public static final String CPX_SV_SYB_ALT_ALLELE_STR = "CPX";
    public static final String CPX_EVENT_ALT_ARRANGEMENTS = "ALT_ARRANGEMENT";
    public static final String CPX_SV_REF_SEGMENTS = "SEGMENTS";
    public static final String CPX_EVENT_KEY = "CPX_EVENT";

    // not defined in output vcf header but used in internal id that is currently output in the ID column
    public static final String INTERVAL_VARIANT_ID_FIELD_SEPARATOR = "_";
    public static final String DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING = "DEL-DUPLICATION-TANDEM-CONTRACTION";
    public static final String DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING = "INS-DUPLICATION-TANDEM-EXPANSION";
    public static final String DUP_INV_INTERNAL_ID_START_STRING = "INS-DUPLICATION-INVERTED-EXPANSION";

    // for breakpoint segmentation
    public static final String ALGORITHMS_ATTRIBUTE = "ALGORITHMS";
    public static final String STRANDS_ATTRIBUTE = "STRANDS";
    public static final String DEPTH_ALGORITHM = "depth";
    public static final String CONTIG2_ATTRIBUTE = "CHR2";
    public static final String END2_ATTRIBUTE = "END2";

    // format block
    public static final String COPY_NUMBER_FORMAT = "CN";
    public static final String COPY_NUMBER_QUALITY_FORMAT = "CNQ";

    // filter block
    public static final String ASSEMBLY_BASED_VARIANT_MQ_FILTER_KEY = "LOW_MQ";
    public static final String ASSEMBLY_BASED_VARIANT_ALN_LENGTH_FILTER_KEY = "SHORT_ALN";
    public static final String LOW_QS_SCORE_FILTER_KEY = "LOW_QS";
    public static final String FREQUENCY_FILTER_KEY = "FREQ";

    // Clustering
    public static final String CLUSTER_MEMBER_IDS_KEY = "MEMBERS";

    // evidence metrics
    public static final String COPY_NUMBER_LOG_POSTERIORS_KEY = "CNLP";
    public static final String NEUTRAL_COPY_NUMBER_KEY = "NCN";
    public static final String DEPTH_OVERLAP_KEY = "RDOV";
    public static final String DEPTH_P_HARDY_WEINBERG_LOSS_FIELD = "PHW_L";
    public static final String DEPTH_P_HARDY_WEINBERG_GAIN_FIELD = "PHW_G";
    public static final String DEPTH_BACKGROUND_FIELD = "ERD";
    public static final String DEPTH_MEAN_BIAS_FIELD = "PHI_RD";
    public static final String START_SPLIT_READ_COUNT_ATTRIBUTE = "SR1";
    public static final String END_SPLIT_READ_COUNT_ATTRIBUTE = "SR2";
    public static final String DISCORDANT_PAIR_COUNT_ATTRIBUTE = "PE";

    // genotyping
    public static final String COPY_NUMBER_FIELD = "CN";
    public static final String PAIRED_END_PROB_FIELD = "PPE";
    public static final String FIRST_SPLIT_READ_PROB_FIELD = "PSR1";
    public static final String SECOND_SPLIT_READ_PROB_FIELD = "PSR2";
    public static final String PAIRED_END_BACKGROUND_FIELD = "EPE";
    public static final String FIRST_SPLIT_READ_BACKGROUND_FIELD = "ESR1";
    public static final String SECOND_SPLIT_READ_BACKGROUND_FIELD = "ESR2";
    public static final String PAIRED_END_MEDIAN_BIAS_FIELD = "PHI_PE";
    public static final String FIRST_SPLIT_READ_MEDIAN_BIAS_FIELD = "PHI_SR1";
    public static final String SECOND_SPLIT_READ_MEDIAN_BIAS_FIELD = "PHI_SR2";
    public static final String MEDIAN_HARDY_WEINBERG_Q_FIELD = "HWQ";
    public static final String MEDIAN_HARDY_WEINBERG_R_FIELD = "HWR";

    public static final String PAIRED_END_BACKGROUND_IQR_FIELD = "EPE_IQR";
    public static final String FIRST_SPLIT_READ_BACKGROUND_IQR_FIELD = "ESR1_IQR";
    public static final String SECOND_SPLIT_READ_BACKGROUND_IQR_FIELD = "ESR2_IQR";
    public static final String PAIRED_END_BIAS_IQR_FIELD = "PHI_PE_IQR";
    public static final String FIRST_SPLIT_READ_BIAS_IQR_FIELD = "PHI_SR1_IQR";
    public static final String SECOND_SPLIT_READ_BIAS_IQR_FIELD = "PHI_SR2_IQR";
    public static final String HARDY_WEINBERG_Q_IQR_FIELD = "HWQ_IQR";
    public static final String HARDY_WEINBERG_R_IQR_FIELD = "HWR_IQR";


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
