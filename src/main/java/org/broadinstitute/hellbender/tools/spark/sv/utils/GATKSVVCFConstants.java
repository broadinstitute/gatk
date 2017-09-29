package org.broadinstitute.hellbender.tools.spark.sv.utils;

public class GATKSVVCFConstants {

    // VCF standard keys reserved for sv
    // todo: add these and the other standard SV info fields from the VCF spec to htsjdk VCFStandardHeaderLines
    public static final String SVTYPE = "SVTYPE";
    public static final String SVLEN = "SVLEN";

    public static final String BREAKEND_STR = "BND";
    public static final String BND_MATEID_STR = "MATEID";

    // symbolic alt alleles
    public static final String INV33 = "INV33";
    public static final String INV55 = "INV55";
    public static final String SYMB_ALT_ALLELE_INV_IN_HEADER = "INV";
    public static final String SYMB_ALT_ALLELE_DEL_IN_HEADER = "DEL";
    public static final String SYMB_ALT_ALLELE_INS_IN_HEADER = "INS";
    public static final String SYMB_ALT_ALLELE_DUP_IN_HEADER = "DUP";

    // GATK-SV specific header lines
    // applicable to all precise records all the time
    public static final String CONTIG_NAMES = "CTG_NAMES";
    public static final String TOTAL_MAPPINGS = "TOTAL_MAPPINGS";
    public static final String MAPPING_QUALITIES = "MAPPING_QUALITIES";
    public static final String ALIGN_LENGTHS = "ALIGN_LENGTHS";
    public static final String HQ_MAPPINGS = "HQ_MAPPINGS";
    public static final String MAX_ALIGN_LENGTH = "MAX_ALIGN_LENGTH";
    // applicable to all precise records when available
    public static final String INSERTED_SEQUENCE = "INSERTED_SEQUENCE";
    public static final String INSERTED_SEQUENCE_MAPPINGS = "INSERTED_SEQUENCE_MAPPINGS";
    public static final String HOMOLOGY = "HOMOLOGY";
    public static final String HOMOLOGY_LENGTH = "HOMOLOGY_LENGTH";
    // type specific: tandem duplication
    public static final String DUP_REPEAT_UNIT_REF_SPAN = "DUP_REPEAT_UNIT_REF_SPAN";
    public static final String DUP_SEQ_CIGARS = "DUP_SEQ_CIGARS";
    public static final String DUPLICATION_NUMBERS = "DUP_NUM";
    public static final String DUP_ANNOTATIONS_IMPRECISE = "DUP_ANNOTATIONS_IMPRECISE";

    public static final String TANDUP_CONTRACTION_STRING = "CONTRACTION";
    public static final String TANDUP_EXPANSION_STRING = "EXPANSION";

    // not defined in output vcf header but used in internal id that is currently output in the ID column
    public static final String INTERVAL_VARIANT_ID_FIELD_SEPARATOR = "_";
    public static final String TANDUP_CONTRACTION_INTERNAL_ID_START_STRING = "DEL-DUPLICATION-TANDEM-CONTRACTION";
    public static final String TANDUP_EXPANSION_INTERNAL_ID_START_STRING = "INS-DUPLICATION-TANDEM-EXPANSION";
}
