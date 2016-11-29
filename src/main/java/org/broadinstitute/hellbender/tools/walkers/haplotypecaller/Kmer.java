package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * Fast wrapper for byte[] kmers
 *
 * This objects has several important features that make it better than using a raw byte[] for a kmer:
 *
 * -- Can create kmer from a range of a larger byte[], allowing us to avoid Array.copyOfRange
 * -- Fast equals and hashcode methods
 * -- can get actual byte[] of the kmer, even if it's from a larger byte[], and this operation
 *    only does the work of that operation once, updating its internal state
 */
public final class Kmer {
    // this values may be updated in the course of interacting with this kmer
    private byte[] bases;
    private int start;

    // two constants
    private final int length;
    private final int hash;

    /**
     * Create a new kmer using all bases in kmer
     * @param kmer a non-null byte[]. The input array must not be modified by the caller.
     */
    public Kmer(final byte[] kmer) {
        this(kmer, 0, kmer.length);
    }

    /**
     * Create a new kmer based on the string kmer
     *
     * This is not a good method to use for performance
     *
     * @param kmer the bases as a string
     */
    public Kmer(final String kmer) {
        this(Utils.nonNull(kmer).getBytes());
    }

    /**
     * Create a new kmer backed by the bases in bases, spanning start -> start + length
     *
     * Under no circumstances can bases be modified anywhere in the client code.  This does not make a copy
     * of bases for performance reasons
     *
     * @param bases an array of bases
     * @param start the start of the kmer in bases, must be >= 0 and < bases.length
     * @param length the length of the kmer.  Must be >= 0 and start + length < bases.length
     */
    public Kmer(final byte[] bases, final int start, final int length) {
        this.bases = Utils.nonNull(bases, "bases cannot be null");
        Utils.validateArg(start >= 0, () -> "start must be >= 0 but got " + start);
        Utils.validateArg( length >= 0, () -> "length must be >= 0 but got " + length);
        Utils.validateArg(start + length <= bases.length, () -> "start + length " + (start + length) + " must be <= bases.length " + bases.length + " but got " + start + " with length " + length);
        this.start = start;
        this.length = length;
        hash = hashCode(bases, start, length);
    }

    /**
     *  Compute the hashcode for a KMer.
     *  Equivalent to <code>new String(bases, start, length).hashCode()</code>
     */
    private static int hashCode(final byte[] bases, final int start, final int length) {
        if (length == 0){
            return 0;
        }
        int h = 0;
        for (int i = start, stop = start + length; i < stop; i++) {
            h = 31 * h + bases[i];
        }
        return h;
    }

    /**
     * Create a derived shallow kmer that starts at newStart and has newLength bases
     * @param newStart the new start of kmer, where 0 means that start of the kmer, 1 means skip the first base
     * @param newLength the new length
     * @return a new kmer based on the data in this kmer.  Does not make a copy, so shares most of the data
     */
    public Kmer subKmer(final int newStart, final int newLength) {
        return new Kmer(bases, start + newStart, newLength);
    }

    /**
     * Get the bases of this kmer.  May create a copy of the bases, depending on how this kmer was constructed.
     *
     * Note that this function is efficient in that if it needs to copy the bases this only occurs once.
     *
     * @return a non-null byte[] containing length() bases of this kmer, regardless of how this kmer was created
     */
    public byte[] bases() {
        if ( start != 0 || bases.length != length ) {
            // update operation.  Rip out the exact byte[] and update start so we don't ever do this again
            bases = Arrays.copyOfRange(bases, start, start + length);
            start = 0;
        }

        return bases;
    }

    /**
     * The length of this kmer
     * @return an integer >= 0
     */
    public int length() {
        return length;
    }

    /**
     * Gets a set of differing positions and bases from another k-mer, limiting up to a max distance.
     * For example, if this = "ACATT" and other = "ACGGT":
     * - if maxDistance < 2 then -1 will be returned, since distance between kmers is 2.
     * - If maxDistance >=2, then 2 will be returned, and arrays will be filled as follows:
     * differingIndeces = {2,3}
     * differingBases = {'G','G'}
     * @param other                 Other k-mer to test
     * @param maxDistance           Maximum distance to search. If this and other k-mers are beyond this Hamming distance,
     *                              search is aborted and -1 is returned
     * @param differingIndeces      Array with indices of differing bytes in array
     * @param differingBases        Actual differing bases
     * @return                      Set of mappings of form (int->byte), where each elements represents index
     *                              of k-mer array where bases mismatch, and the byte is the base from other kmer.
     *                              If both k-mers differ by more than maxDistance, returns null
     */
    public int getDifferingPositions(final Kmer other,
                                                   final int maxDistance,
                                                   final int[] differingIndeces,
                                                   final byte[] differingBases) {

        Utils.nonNull(other);
        Utils.nonNull(differingIndeces);
        Utils.nonNull(differingBases);
        Utils.validateArg(maxDistance > 0, "maxDistance must be positive but was " + maxDistance);
        int dist = 0;
        if (length == other.length()) {
            final byte[] f2 = other.bases;
            for (int i=0; i < length; i++) {
                if (bases[start + i] != f2[i]) {
                    differingIndeces[dist] = i;
                    differingBases[dist++] = f2[i];
                    if (dist > maxDistance) {
                        return -1;
                    }
                }
            }

        }
        return dist;
    }

    @Override
    public String toString() {
        return "Kmer{" + new String(bases,start,length) + '}';
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof Kmer)) {
            return false;
        }

        final Kmer kmer = (Kmer) o;

        // very fast test.  If hash aren't equal you are done, otherwise compare the bases
        if ( hash != kmer.hash || length != kmer.length ) {
            return false;
        }

        return Utils.equalRange(bases, start, kmer.bases, kmer.start, length);
    }

    @Override
    public int hashCode() {
        return hash;
    }

    @VisibleForTesting
    int start() {
        return start;
    }
}
