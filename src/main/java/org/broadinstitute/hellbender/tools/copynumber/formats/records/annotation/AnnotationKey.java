package org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Function;

/**
 * Represents a key for a named, typed annotation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AnnotationKey<T> {
    private final String name;
    private final Class<T> clazz;
    private final Function<T, Boolean> validateValue;

    public AnnotationKey(final String name,
                         final Class<T> clazz,
                         final Function<T, Boolean> validateValue) {
        this.name = Utils.nonEmpty(name);
        this.clazz = Utils.nonNull(clazz);
        this.validateValue = Utils.nonNull(validateValue);
    }

    public String getName() {
        return name;
    }

    public Class<T> getType() {
        return clazz;
    }

    public T validate(final T value) {
        Utils.validateArg(validateValue.apply(value),
                String.format("Invalid value %s for annotation %s.", value, name));
        return value;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AnnotationKey<?> that = (AnnotationKey<?>) o;
        return name.equals(that.name) && clazz.equals(that.clazz);
    }

    @Override
    public int hashCode() {
        int result = name.hashCode();
        result = 31 * result + clazz.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "AnnotationKey{" +
                "name='" + name + '\'' +
                ", class=" + clazz +
                '}';
    }
}
