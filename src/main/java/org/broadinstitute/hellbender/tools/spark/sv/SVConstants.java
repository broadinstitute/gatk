package org.broadinstitute.hellbender.tools.spark.sv;

/**
 * Constants shared among SV tools.
 */
public final class SVConstants {
    private SVConstants() {}

    public static final int KMER_SIZE = 51;
    public static final int MAX_DUST_SCORE = KMER_SIZE - 2;
}
