package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

/**
 * A Duplicable wrapper around Numbers.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class DuplicableNumber<N extends Number> implements Duplicable {

    final private N value;

    public DuplicableNumber() {
        value = null;
    }

    public DuplicableNumber(final N value) {
        this.value = value;
    }

    @Override
    public DuplicableNumber<N> deepCopy() {
        return this;
    }

    @Override
    public boolean isNull() {
        return value == null;
    }
    public N value() {
        return value;
    }

    public static double of(final Duplicable obj) {
        if (obj instanceof DuplicableNumber) {
            return ((DuplicableNumber)obj).value().doubleValue();
        } else {
            throw new ClassCastException("Can not cast " + obj + " to an INDArray.");
        }
    }

    @Override
    public String toString() {
        if (value == null) {
            return "null";
        } else {
            return value.toString();
        }
    }

}