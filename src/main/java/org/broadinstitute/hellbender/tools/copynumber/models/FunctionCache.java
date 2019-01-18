package org.broadinstitute.hellbender.tools.copynumber.models;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class FunctionCache<DATA> extends LinkedHashMap<DATA, Double> {
    private static final long serialVersionUID = 19841647L;
    private static final int MAX_SIZE = 100_000;

    private final Function<DATA, Double> mappingFunction;

    FunctionCache(final Function<DATA, Double> mappingFunction) {
        this.mappingFunction = mappingFunction;
    }

    Double computeIfAbsent(final DATA key) {
        return super.computeIfAbsent(key, mappingFunction);
    }

    @Override
    protected boolean removeEldestEntry(final Map.Entry<DATA, Double> eldest) {
        return size() >= MAX_SIZE;
    }
}
