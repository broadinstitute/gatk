package org.broadinstitute.hellbender.utils.diffengine;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * An interface that must be implemented to allow us to calculate differences
 * between structured objects
 */
class DiffValue {
    private DiffElement binding = null;
    final private Object value;

    public DiffValue(final Object value) {
        this.value = value;
    }

    public DiffValue(final DiffElement binding, final Object value) {
        this.binding = binding;
        this.value = value;
    }

    public DiffValue(final DiffValue parent, final Object value) {
        this(parent.getBinding(), value);
    }

    public DiffValue(final String name, final DiffElement parent, final Object value) {
        this.binding = new DiffElement(name, parent, this);
        this.value = value;
    }

    public DiffValue(final String name, final DiffValue parent, final Object value) {
        this(name, parent.getBinding(), value);
    }

    public DiffElement getBinding() {
        return binding;
    }

    protected void setBinding(final DiffElement binding) {
        this.binding = binding;
    }

    public Object getValue() {
        return value;
    }

    @Override
    public String toString() {
        return getValue().toString();
    }

    public String toString(final int offset) {
        return toString();
    }

    public String toOneLineString() {
        return getValue().toString();
    }

    public boolean isAtomic() { return true; }

    public final boolean isCompound() { return ! isAtomic(); }

    public int size() { return 1; }

    public List<Difference> diffValues(final DiffValue test) {
        if ( this.getValue().equals(test.getValue()) ) {
            return Collections.emptyList();
        } else {
            return Arrays.asList(new Difference(this.getBinding(), test.getBinding()));
        }
    }
}
