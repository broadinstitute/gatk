package org.broadinstitute.hellbender.tools.spark.sv;

/**
 * Constants shared among SV tools.
 */
public final class SVConstants {
    private SVConstants() {}

    public static final int KMER_SIZE = 51;
    public static final double MIN_ENTROPY = 1.25;

    public static final String FASTQ_OUT_PREFIX = "assembly";
}
