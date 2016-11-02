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

    /**
     * Minimum flanking alignment length used in calling variants when going through contig alignments.
     */
    static final Integer DEFAULT_MIN_ALIGNMENT_LENGTH = 50;

    static final String VARIANT_ID_FIELD_SEPARATER = ";";

    static final int CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD = 60;
}
