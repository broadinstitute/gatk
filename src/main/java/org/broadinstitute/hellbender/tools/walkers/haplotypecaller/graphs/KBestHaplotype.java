package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

/**
 * Represents a result from a K-best haplotype search.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class KBestHaplotype implements Comparable<KBestHaplotype> {

    /**
     * Returns the original graph searched.
     *
     * @return never {@code null}
     */
    public abstract SeqGraph graph();

    /**
     * Returns the result haplotype score.
     *
     * <p>Currently, the score is the multiplicity total sum of edges along the haplotype path</p>
     *
     * @return 0 or greater.
     */
    public abstract double score();

    /**
     * Indicates whether this result is the reference haplotype.
     *
     * @return {@code true} if it is the reference haplotype, {@code false} otherwise.
     */
    public abstract boolean isReference();

    /**
     * The rank of this solution within the list of solutions that resulted from the same search.
     *
     * <p>0 would correspond to the best solution, 1 with the second best and so on</p>
     *
     * @return 0 or greater.
     */
    public abstract int rank();

    private byte[] bases;

    private Haplotype haplotype;

    private Path<SeqVertex,BaseEdge> path;

    /**
     * Returns the result haplotype bases.
     *
     * @return never {@code null}.
     */
    public final byte[] bases() {
        if (bases != null) {
            return bases;
        }
        final KBestHaplotype tail = tail();
        final SeqVertex head = head();
        if (tail == null) {
            bases = head.getSequence();
        } else {
            final byte[] tailBases = tail.bases();
            final byte[] headBases = head.getSequence();
            final int length = tailBases.length + headBases.length;
            bases = new byte[length];
            System.arraycopy(headBases, 0, bases, 0, headBases.length);
            System.arraycopy(tailBases, 0, bases, headBases.length, tailBases.length);
        }
        return bases;
    }

    /**
     * Returns the solution haplotype.
     *
     * @return never {@code null}.
     */
    public final Haplotype haplotype() {
        if (haplotype != null) {
            return haplotype;
        }
        haplotype = new Haplotype(bases(),isReference());
        if (score() > 0) {
            throw new IllegalStateException("score cannot be greater than 0: " + score());
        }
        haplotype.setScore(score());
        return haplotype;
    }

    /**
     * Returns the path across the original graph that correspond to the solution haplotype.
     *
     * @return never {@code null}, although perhaps a zero-length path (only one vertex).
     */
    public final Path<SeqVertex,BaseEdge> path() {
        if (path != null) {
            return path;
        }
        final KBestHaplotype tail = tail();
        if (tail == null) {
            path = new Path<>(head(), graph());
        } else {
            final Path<SeqVertex,BaseEdge> tailPath = tail.path();
            path = new Path<>(graph().getEdge(head(),tailPath.getFirstVertex()),tailPath);
        }
        return path;
    }

    /**
     * Compares k-best haplotypes based on the score where the one with larger score comes first (descending orther).
     *
     * @param other the other haplotype to compare to.
     * @return {@code -1} if the current score is larger than {@code other}'s, {@code 0} if they are the same, {@code 1}
     * if {@code other}'s score is larger.
     */
    @Override
    public final int compareTo(final KBestHaplotype other) {
        Utils.nonNull(other, "the other object cannot be null");
        return -1 * Double.compare(score(), other.score());
    }

    @Override
    public final int hashCode() {
        return haplotype().hashCode();
    }

    @Override
    public final boolean equals(final Object other) {
        if (!(other instanceof KBestHaplotype)) {
            return false;
        }
        return equals((KBestHaplotype) other);
    }

    @Override
    public final String toString() {
        return haplotype() + " Score = "  + score();
    }

    /**
     * Checks whether both solutions are equal.
     * <p>
     *     Both solutions are considered equal when the underlying haplotypes are equal. The path on the respective
     *     graph might deffer though.
     * </p>
     *
     * @return {@code true} iff both haplotypes are the same (considering the ref state).
     */
    private boolean equals(final KBestHaplotype other) {
       return haplotype().equals(other.haplotype(),false);
    }

    /**
     * The first vertex on the haplotype path.
     *
     * @return never {@code null}.
     */
    protected abstract SeqVertex head();

    /**
     * Returns the sub-haplotype from the second vertex involved in the haplotype until the end.
     *
     * @return {@code null} if there are no more vertices in the solution path a part from the one returned by {@link #head}.
     */
    protected abstract KBestHaplotype tail();
}
