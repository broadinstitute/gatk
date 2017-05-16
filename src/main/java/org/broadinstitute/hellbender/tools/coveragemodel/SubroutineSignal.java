package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * This class is used for communicating compound exit signals from computational subroutines.
 * It is essentially a wrapper around a map.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SubroutineSignal implements Serializable {

    private static final long serialVersionUID = -8990591574704241500L;

    private final Map<String, Object> result;

    public static final SubroutineSignal EMPTY_SIGNAL = new SubroutineSignal(new HashMap<>());

    /**
     * Private constructor
     * @param result a key-value map for exit signals
     */
    private SubroutineSignal(final Map<String, Object> result) {
        this.result = result;
    }

    @SuppressWarnings("unchecked")
    public <TYPE> TYPE get(final String key) {
        Utils.validateArg(result.containsKey(Utils.nonNull(key, "The subroutine signal key must be non-null")),
                "No exit signal is available for \"" + key + "\"");
        try {
            return (TYPE) result.get(key);
        } catch (final ClassCastException ex) {
            throw new UnsupportedOperationException("Can not cast the value associated to \"" + key + "\" to the" +
                    " inferred generic type");
        }
    }

    public <TYPE> TYPE get(@Nonnull final StandardSubroutineSignals key) {
        return get(Utils.nonNull(key, "The subroutine signal key must be non-null").name());
    }

    /**
     * Creates an instance of {@link SubroutineSignalBuilder}
     */
    public static SubroutineSignalBuilder builder() {
        return new SubroutineSignalBuilder();
    }

    /**
     * Static builder
     */
    public static final class SubroutineSignalBuilder implements Serializable {

        public static final long serialVersionUID = 1190387682156562190L;

        private final Map<String, Object> result;

        SubroutineSignalBuilder() {
            result = new HashMap<>();
        }

        public SubroutineSignalBuilder put(final String key, final Object value) {
            result.put(key, value);
            return this;
        }

        public SubroutineSignalBuilder put(final StandardSubroutineSignals key, final Object value) {
            result.put(key.name(), value);
            return this;
        }

        public SubroutineSignal build() {
            return new SubroutineSignal(result);
        }
    }
}