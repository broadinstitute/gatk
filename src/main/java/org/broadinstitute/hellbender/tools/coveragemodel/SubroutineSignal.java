package org.broadinstitute.hellbender.tools.coveragemodel;

import org.nd4j.linalg.api.ndarray.INDArray;

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

    /**
     * Private constructor
     * @param result a key-value map for exit signals
     */
    private SubroutineSignal(final Map<String, Object> result) {
        this.result = result;
    }

    /**
     * Fetch a double from exit signals
     *
     * @param key string identifier
     * @return a double
     */
    public double getDouble(final String key) {
        if (result.containsKey(key)) {
            return (double)result.get(key);
        } else {
            throw new IllegalArgumentException("No exit signal is available for \"" + key + "\"");
        }
    }

    /**
     * Fetch an integer from exit signals
     *
     * @param key string identifier
     * @return an integer
     */
    public int getInteger(final String key) {
        if (result.containsKey(key)) {
            return (int)result.get(key);
        } else {
            throw new IllegalArgumentException("No exit signal is available for \"" + key + "\"");
        }
    }

    /**
     * Fetch a String from exit signals
     *
     * @param key string identifier
     * @return a String
     */
    public String getString(final String key) {
        if (result.containsKey(key)) {
            return (String)result.get(key);
        } else {
            throw new IllegalArgumentException("No exit signal is available for \"" + key + "\"");
        }
    }

    /**
     * Fetch an INDArray from exit signals
     *
     * @param key string identifier
     * @return an {@link INDArray}
     */
    public INDArray getINDArray(final String key) {
        if (result.containsKey(key)) {
            return (INDArray)result.get(key);
        } else {
            throw new IllegalArgumentException("No exit signal is available for \"" + key + "\"");
        }
    }

    /**
     * Fetch an Object from exit signals
     *
     * @param key string identifier
     * @return an Object
     */
    public Object getObject(final String key) {
        if (result.containsKey(key)) {
            return result.get(key);
        } else {
            throw new IllegalArgumentException("No exit signal is available for \"" + key + "\"");
        }
    }

    /**
     * Creates an instance of {@link SubroutineSignalBuilder}
     *
     * @return an instance of {@link SubroutineSignalBuilder}
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

        public SubroutineSignal build() {
            return new SubroutineSignal(result);
        }
    }
}