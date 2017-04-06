package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;

import java.io.Serializable;
import java.util.Collections;
import java.util.Comparator;
import java.util.Objects;

/**
 * Parent class for {@link Locatable locatable} {@code Comparator comparators} that adheres to the standard position order
 * within a contig.
 *
 * <p>
 *     All descendant classes must respect the following order of {@link Locatable}s within a contig:
 *     <ul>
 *         <li>instances with a lower {@code start} position come first</li>
 *         <li>then if two instances have the same {@code start} location the ones with a smaller {@code end} location come first.</li>
 *     </ul>
 * </p>
 * <p>
 *     Also, the {@code null} contig is considered the last contig and so all instances with a non-null contig would come before
 *     any instances with a {@code null} contig.
 * </p>
 * <p>
 *     The relative order between contig is what is left to be determined by the implementing class by defining {@link #compareContigs(String,String)}.
 * </p>
 */
public abstract class LocatableComparator implements Comparator<Locatable>, Serializable {

    /**
     * Comparator that
     */
    public static final LocatableComparator LEXICOGRAPHIC = new LocatableComparator() {

        private static final long serialVersionUID = 1L;

        @Override
        protected int compareContigs(String leftContig, String rightContig) {
            return Utils.nonNull(leftContig).compareTo(Utils.nonNull(rightContig));
        }
    };

    private static final long serialVersionUID = 1L;

    @Override
    public final int compare(final Locatable left, final Locatable right) {
        final String leftContig = left.getContig();
        final String rightContig = right.getContig();
        if (Objects.equals(leftContig, rightContig)) {
            final int startCmp = Integer.compare(left.getStart(), right.getStart());
            return startCmp == 0 ? Integer.compare(left.getEnd(), right.getEnd()) : startCmp;
        } else if (leftContig == null) {
            return 1;
        } else if (rightContig == null) {
            return -1;
        } else {
            return compareContigs(leftContig, rightContig);
        }
    }

    public <L extends Locatable> Comparator<L> narrow(final Class<L> clazz) {
        final LocatableComparator that = this;
        return new Comparator<L>() {
            @Override
            public int compare(final L o1, final L o2) {
                Utils.nonNull(o1);
                Utils.nonNull(o2);
                if (!clazz.isAssignableFrom(o1.getClass())) {
                    throw new IllegalArgumentException("bad left object class " + o1.getClass());
                } else if (!clazz.isAssignableFrom(o2.getClass())) {
                    throw new IllegalArgumentException("bad right object class " + o1.getClass());
                }
                return that.compare(o1, o2);
            }
        };
    }

    /**
     * Determines the relative order of contigs. It must follow the common properties of
     * any ordering function such as transitivity and reciprocativity.
     * <p>
     *     The behavior of this method if any of the arguments in {@code null} is undetermined.
     * </p>
     *
     * @param leftContig the left contig to compare.
     * @param rightContig the right contig to compare.
     *
     * @return {@code 0} if left and right are the same, {@code -1} if left comes first and {@code 1} if right come first.
     */
    protected abstract int compareContigs(final String leftContig, final String rightContig);
}
