package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents a parameter value with a named {@link ParameterEnum} key.
 * @param <T>   type of the {@link ParameterEnum} key name
 * @param <U>   type of the value
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Parameter<T extends Enum<T> & ParameterEnum, U> {
    private final T name;
    private final U value;

    /**
     * Constructs a {@link Parameter} given a {@link ParameterEnum} key name and value.
     * @param name  parameter key name
     * @param value parameter value
     * @throws IllegalArgumentException if {@code name} or {@code value} is null
     */
    public Parameter(final T name, final U value) {
        this.name = Utils.nonNull(name, "The parameter name cannot be null.");
        this.value = Utils.nonNull(value, "The parameter value cannot be null.");
    }

    public T getName() {
        return name;
    }

    public U getValue() {
        return value;
    }
}
