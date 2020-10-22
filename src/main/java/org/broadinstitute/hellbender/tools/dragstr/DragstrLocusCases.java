package org.broadinstitute.hellbender.tools.dragstr;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Collection of Dragstr Locus cases constraint to a particular period and (minimum) repeat-length
 */
public final class DragstrLocusCases extends AbstractList<DragstrLocusCase> {

    private final List<DragstrLocusCase> cases;

    private final int period;

    private final int repeatLength;


    /**
     * Creates a new collection with default initial capacity.
     * @param period the collection cases period.
     * @param repeatLength the collection cases (minimum) repeat-length.
     */
    public DragstrLocusCases(final int period, final int repeatLength) {
        this(100, period, repeatLength);
    }

    @Override
    public DragstrLocusCase get(final int index) {
        return cases.get(index);
    }

    /**
     * Creates a new dragstr cases collection providing it period and repeat-length constraints.
     * @param initialCapacity the expected maximum number of members in the collection.
     * @param period the member case period constraint. Only cases with this period are allowed.
     * @param repeatLength the member case repeat-length constraint. Only cases with a repeat-length
     *                     as large are allowed.
     */
    public DragstrLocusCases(final int initialCapacity, final int period, final int repeatLength) {
        Utils.validateArg( initialCapacity >= 0, "the initial capacity must be 0 or greater");
        Utils.validateArg(period >= 1, "the period cannot be less than 1");
        Utils.validateArg(repeatLength >= 1, "the repeat-length cannot be less than 1");
        this.period = period;
        this.repeatLength = repeatLength;
        cases = new ArrayList<>(initialCapacity);
    }

    private DragstrLocusCases(final int period, final int repeatLength, final List<DragstrLocusCase> cases) {
        this.period = period;
        this.repeatLength = repeatLength;
        this.cases = cases;
    }

    /**
     * Adds a case to the collection.
     * <p>
     *     This method makes sure that you don't add cases with a different period. The input repeat-length might be
     *     larger than the collection's to accomodate the special case of the "maximum" analyzed repeat-length.
     * </p>
     * @param caze the case to add.
     * @throws IllegalArgumentException if the input is {@code null} or has the wrong period or repeat-length.
     * @return always {@code true}.
     */
    public boolean add(final DragstrLocusCase caze) {
        Utils.validateArg(caze.getPeriod() == period, "the input case has the wrong period");
        // We allow for larger repeat lengths since that might be the case for the "maximum" analyzed repeat-length.
        Utils.validateArg(caze.getRepeatLength() >= repeatLength, "the input case has a repeat-length smaller than this collections");
        return cases.add(caze);
    }

    public boolean addAll(final Collection<? extends DragstrLocusCase> more) {
        for (final DragstrLocusCase caze : more) {
            add(caze);
        }
        return true;
    }

    public boolean addAll(final DragstrLocusCases other) {
        return addAll(other.cases);
    }

    public int getPeriod() {
        return period;
    }

    public int getRepeatLength() {
        return repeatLength;
    }

    public int size() {
        return cases.size();
    }

    /**
     * Returns a new collection that contains only those cases are qualify for parameter estimation given
     * the required constraints.
     * @param minDepth minimum depth admissible.
     * @param minMQ minimum mapping quality admissible.
     * @param maxSup maximum number of supplementary alignments admissible.
     * @return never {@code null}.
     */
    public DragstrLocusCases qualifyingOnly(final int minDepth, final int minMQ, final int maxSup) {
       final List<DragstrLocusCase> newCases = cases.stream()
                    .filter(caze -> caze.getDepth() >= minDepth && caze.getMinMQ() >= minMQ && caze.getNSup() <= maxSup)
                    .collect(Collectors.toList());
        return new DragstrLocusCases(period, repeatLength, newCases);
    }
}
