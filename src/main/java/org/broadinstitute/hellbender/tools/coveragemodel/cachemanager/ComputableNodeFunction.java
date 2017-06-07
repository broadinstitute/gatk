package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.io.Serializable;
import java.util.Map;

/**
 * A functional interface for computation functions of a {@link ComputableCacheNode}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface ComputableNodeFunction {

    Duplicable apply(final Map<CacheNode.NodeKey, Duplicable> parents) throws ParentValueNotFoundException;

    /**
     * Fetches a parent node value from a given map
     *
     * @param key parent key
     * @param parents parent key-value map
     *
     * @return an instance of {@link Duplicable} by reference
     * @throws ParentValueNotFoundException if the parent key is not in the map
     */
    default Duplicable fetch(final CacheNode.NodeKey key, final Map<CacheNode.NodeKey, Duplicable> parents)
            throws ParentValueNotFoundException {
        if (!parents.containsKey(key)) {
            throw new ParentValueNotFoundException(key);
        }
        return parents.get(key);
    }

    /**
     * Fetches a parent node value from a given map and casts it to an INDArray
     *
     * @param key parent key
     * @param parents parent key-value map
     * @throws ParentValueNotFoundException if the parent key is not in the map
     */
    default INDArray fetchINDArray(final CacheNode.NodeKey key, final Map<CacheNode.NodeKey, Duplicable> parents)
            throws ParentValueNotFoundException, ClassCastException {
        return ((DuplicableNDArray)fetch(key, parents)).value();
    }

    /**
     * Fetches a parent node value from a given map and casts it to a double
     *
     * @param key parent key
     * @param parents parent key-value map
     * @throws ParentValueNotFoundException if the parent key is not in the map
     */
    @SuppressWarnings("unchecked")
    default double fetchDouble(final CacheNode.NodeKey key, final Map<CacheNode.NodeKey, Duplicable> parents)
            throws ParentValueNotFoundException, ClassCastException {
        return ((DuplicableNumber<Double>)fetch(key, parents)).value();
    }

    /**
     * This exception will be thrown if a required parent node is not in the map
     * supplied to {@link #apply(Map)}
     */
    final class ParentValueNotFoundException extends RuntimeException implements Serializable {
        private static final long serialVersionUID = -4557250891066141519L;

        private ParentValueNotFoundException(final CacheNode.NodeKey nodeKey) {
            super(String.format("The value of node \"%s\" is required for computation but it is not available", nodeKey));
        }
    }
}
