package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;

/**
 * Common interface for implementations of the core Smith Waterman pairwise alignment algorithm.
 * The API is designed to be used as follows.
 * First, {@link #initialize(SmithWatermanParameters)} is called once,
 * then for each pairs of sequences, {@link #align(byte[], byte[])} is called, followed by optional
 * {@link #getCigar()} and {@link #getAlignmentStart2wrt1()}.
 */
public interface CoreSmithWatermanAlgorithm {

    /**
     * Initialize the core algorithm. Returns true if this implementation supports
     * the given parameter values, false otherwise.
     * For example, not all implementations may support all values of {@link OverhangStrategy}.
     */
    public boolean initialize(SmithWatermanParameters params);

    /**
     * Aligns the alternate sequence to the reference sequence
     */
    public void align(final byte[] reference, final byte[] alternate);

    public Cigar getCigar();

    public int getAlignmentStart2wrt1();
}
