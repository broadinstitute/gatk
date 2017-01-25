package org.broadinstitute.hellbender.tools.spark.sv;

/**
 * Constants shared among SV tools.
 */
public final class SVConstants {
    private SVConstants() {
    }

    public static final int KMER_SIZE = 51;
    public static final int MAX_DUST_SCORE = KMER_SIZE - 2;

    // -----------------------------------------------------------------------------------------------
    // below are a bunch of heuristics used in various places in the SV pipeline (most likely artifacts of prototyping)
    // (could have made them live closer to their use places but choose to make it easier to identify room for improvements)
    // -----------------------------------------------------------------------------------------------

    static final class CallingStepConstants {

        static final String INVERSIONS_OUTPUT_VCF = "inversions.vcf";
        static final String CURRENTLY_CAPABLE_VARIANTS_VCF = "inv_del_ins.vcf";
        static final String VCF_ALT_ALLELE_STRING_INV = "<INV>";
        static final String VCF_ALT_ALLELE_STRING_INS = "<INS>";
        static final String VCF_ALT_ALLELE_STRING_DEL = "<DEL>";
        static final String VCF_ALT_ALLELE_STRING_DUP = "<DUP>";
        static final String TANDUP_CONTRACTION_STRING = "CONTRACTION";
        static final String TANDUP_EXPANSION_STRING = "EXPANSION";


        static final Integer DEFAULT_MIN_ALIGNMENT_LENGTH = 50; // Minimum flanking alignment length filters used when going through contig alignments.
        static final String VARIANT_ID_FIELD_SEPARATOR = ";";
        static final int CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD = 60;
        static final int GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY = 50; // alignment with gap of size >= 50 will be broken apart.
        static final int MISSING_NM = Integer.MIN_VALUE;
        static final int ARTIFICIAL_MISMATCH = MISSING_NM;
    }
}
