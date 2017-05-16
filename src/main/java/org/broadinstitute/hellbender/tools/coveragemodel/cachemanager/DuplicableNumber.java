package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

/**
 * A {@link Duplicable} wrapper around a {@link Number}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class DuplicableNumber<NUMBER extends Number> implements Duplicable {

    private final NUMBER value;

    public DuplicableNumber() {
        value = null;
    }

    public DuplicableNumber(final NUMBER value) {
        this.value = value;
    }

    @Override
    public DuplicableNumber<NUMBER> duplicate() {
        return this;
    }

    @Override
    public boolean hasValue() {
        return value == null;
    }

    public NUMBER value() {
        return value;
    }

    public static double of(final Duplicable obj) {
        if (obj instanceof DuplicableNumber) {
            return ((DuplicableNumber)obj).value().doubleValue();
        } else {
            throw new ClassCastException("Can not cast " + obj + " to a number.");
        }
    }

    @Override
    public String toString() {
        return value == null ? "null" : value.toString();
    }
}