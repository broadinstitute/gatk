package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Set;

/**
 * Represents a trivial k-best sub haplotype finder with no solutions.
 *
 * <p>To be used at vertices that do not have any valid path to the requested sink vertices</p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class DeadEndKBestSubHaplotypeFinder implements KBestSubHaplotypeFinder {

    /**
     * Sole instance of this class.
     */
    public static final DeadEndKBestSubHaplotypeFinder INSTANCE = new DeadEndKBestSubHaplotypeFinder();

    /**
     * Prevents instantiation of more than one instance; please use {@link #INSTANCE}.
     */
    private DeadEndKBestSubHaplotypeFinder() {
    }

    @Override
    public String id() {
        return "<DEAD>";
    }

    @Override
    public String label() {
        return "&lt;DEAD&gt;";
    }

    @Override
    public Set<Pair<? extends KBestSubHaplotypeFinder, String>> subFinderLabels() {
        return Collections.emptySet();
    }

    @Override
    public int getCount() {
        return 0;
    }

    @Override
    public KBestHaplotype getKBest(final int k) {
        if (k < 0) {
            throw new IllegalArgumentException("k cannot be negative");
        } else {
            throw new IllegalArgumentException("k cannot be equal or greater to the haplotype count");
        }
    }

    @Override
    public boolean isReference() {
        return false;
    }

    @Override
    public double score(final byte[] bases, final int offset, final int length) {
        Utils.nonNull(bases, "bases cannot be null");
        if (offset < 0) {
            throw new IllegalArgumentException("the offset cannot be negative");
        }
        if (length < 0) {
            throw new IllegalArgumentException("the length cannot be negative");
        }
        if (offset + length > bases.length) {
            throw new IllegalArgumentException("the offset and length go beyond the array size");
        }
        return Double.NaN;
    }
}
