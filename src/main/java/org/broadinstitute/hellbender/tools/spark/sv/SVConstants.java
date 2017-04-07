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

    public static final class DiscoveryStepConstants {

        public static final String VCF_ALT_ALLELE_STRING_INV = "<INV>";
        public static final String VCF_ALT_ALLELE_STRING_INS = "<INS>";
        public static final String VCF_ALT_ALLELE_STRING_DEL = "<DEL>";
        public static final String VCF_ALT_ALLELE_STRING_DUP = "<DUP>";
        public static final String TANDUP_CONTRACTION_STRING = "CONTRACTION";
        public static final String TANDUP_EXPANSION_STRING = "EXPANSION";


        public static final Integer DEFAULT_MIN_ALIGNMENT_LENGTH = 50; // Minimum flanking alignment length filters used when going through contig alignments.
        public static final String VARIANT_ID_FIELD_SEPARATOR = "_";
        public static final int CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD = 60;
        public static final int GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY = 50; // alignment with gap of size >= 50 will be broken apart.
        public static final int MISSING_NM = Integer.MIN_VALUE;
        public static final int ARTIFICIAL_MISMATCH = MISSING_NM;
    }
}
