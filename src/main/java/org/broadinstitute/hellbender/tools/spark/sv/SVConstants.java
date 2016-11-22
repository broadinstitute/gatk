package org.broadinstitute.hellbender.tools.spark.sv;

/**
 * Constants shared among SV tools.
 */
public final class SVConstants {
    private SVConstants() {}

    public static final int KMER_SIZE = 51;
    public static final double MIN_ENTROPY = 1.25;

    // -----------------------------------------------------------------------------------------------
    // below are a bunch of heuristics used in various places in the SV pipeline (most likely artifacts of prototyping)
    // (could have made them live closer to their use places but choose to make it easier to identify room for improvements)
    // -----------------------------------------------------------------------------------------------

    static final String INVERSIONS_OUTPUT_VCF = "inversions.vcf";

    static final String INS_DEL_OUTPUT_VCF = "del_ins.vcf";

    static final String VCF_ALT_ALLELE_STRING_INV = "<INV>";

    static final String VCF_ALT_ALLELE_STRING_INS = "<INS>";

    static final String VCF_ALT_ALLELE_STRING_DEL = "<DEL>";

    static final String VCF_ALT_ALLELE_STRING_INDEL = "<INS_DEL>";


    /**
     * Minimum flanking alignment length used in calling variants when going through contig alignments.
     */
    static final Integer DEFAULT_MIN_ALIGNMENT_LENGTH = 50;

    static final String VARIANT_ID_FIELD_SEPARATER = ";";

    static final int CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD = 60;

    static final int GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY = 50; // alignment with gap of size >= 50 will be break apart.
}
