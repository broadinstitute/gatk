package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;

/**
 * Classes that provide a way to test kmers for set membership and keep track of the kmer size and mask
 */
public abstract class PSKmerCollection {

    abstract boolean contains(final SVKmerShort val);
    abstract int kmerSize();
    abstract SVKmerShort getMask();
    abstract double getFalsePositiveProbability();

    /**
     * Definition for the order of canonicalization and masking
     */
    public final static long canonicalizeAndMask(final SVKmerShort val, final int kmerSize, final SVKmerShort mask) {
        return val.canonical(kmerSize).mask(mask).getLong();
    }

}
