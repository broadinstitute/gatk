package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents a named parameter.
 * @param <T>   type of the parameter
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Parameter<T> {
    private final String name;
    private final T value;
    private final Class cls;

    /**
     * Constructs a Parameter given a name and value.
     * @param name  parameter name
     * @param value parameter value
     * @throws IllegalArgumentException if {@code name} or {@code value} is null
     */
    public Parameter(final String name, final T value) {
        Utils.nonNull(name, "The parameter name cannot be null.");
        Utils.nonNull(value, "The parameter value cannot be null.");
        this.name = name;
        this.value = value;
        this.cls = value.getClass();
    }

    /**
     * Returns the parameter name.
     * @return  parameter name
     */
    public String name() {
        return name;
    }

    /**
     * Returns the parameter value.
     * @return  parameter value
     */
    public T value() {
        return value;
    }

    /**
     * Returns the parameter class.
     * @return  parameter class
     */
    public Class cls() {return cls;}
}
