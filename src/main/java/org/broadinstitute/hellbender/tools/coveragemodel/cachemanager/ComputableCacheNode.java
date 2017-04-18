package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.Collection;
import java.util.Map;

/**
 * This class represents an automatically computable cache node.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ComputableCacheNode extends CacheNode {

    private final boolean cacheEvals;
    private final ComputableNodeFunction func;
    private Duplicable cachedValue = null;
    private boolean isCacheCurrent;

    /**
     * Public constructor
     *
     * @param key the key of the node
     * @param parents immediate parents of the node
     * @param func a function from a map that (at least) contains parents data to the computed value of this node
     * @param cacheEvals does it store the value or not
     */
    public ComputableCacheNode(@Nonnull final String key,
                               @Nonnull final Collection<String> tags,
                               @Nonnull final Collection<String> parents,
                               @Nullable final ComputableNodeFunction func,
                               final boolean cacheEvals) {
        super(key, tags, parents);
        this.func = func;
        this.cacheEvals = cacheEvals;
        Utils.validateArg(func != null || cacheEvals, "A computable node with null evaluation function is externally" +
                " mutable and must cache its values");
        isCacheCurrent = false;
    }

    private ComputableCacheNode(@Nonnull final String key,
                                @Nonnull final Collection<String> tags,
                                @Nonnull final Collection<String> parents,
                                @Nullable final ComputableNodeFunction func,
                                final boolean cacheEvals,
                                final Duplicable cachedValue,
                                final boolean isCacheCurrent) {
        super(key, tags, parents);
        this.func = func;
        this.cacheEvals = cacheEvals;
        this.isCacheCurrent = isCacheCurrent;
        this.cachedValue = cachedValue;
    }

    @Override
    public boolean isPrimitive() { return false; }

    @Override
    public boolean isExternallyComputable() { return func == null; }

    public boolean isCacheCurrent() { return isCacheCurrent; }

    public boolean doesCacheEvaluations() { return cacheEvals; }

    /**
     * Available means (1) the node caches its value, and (2) a value is already cached (though may
     * not be up-to-date)
     *
     * @return a boolean
     */
    @Override
    public boolean isStoredValueAvailable() {
        return cacheEvals && cachedValue != null && !cachedValue.hasValue();
    }

    /**
     * In addition to being available, this methods checks if the cached value is up-to-date
     *
     * @return a boolean
     */
    public boolean isStoredValueAvailableAndCurrent() {
        return isStoredValueAvailable() && isCacheCurrent();
    }

    @Override
    public void set(@Nullable final Duplicable val) {
        if (isExternallyComputable()) {
            cachedValue = val;
            isCacheCurrent = true;
        } else {
            throw new UnsupportedOperationException("Can not explicitly set the value of a computable cache node with" +
                    " non-null function.");
        }
    }

    /**
     * Returns the value of the node, either from cache or by computation
     *
     * @param parentsValues a lookup map for parents values
     * @return value
     * @throws ExternallyComputableNodeValueUnavailableException if the node is externally mutable but its value
     *   is out of date
     * @throws ComputableNodeFunction.ParentValueNotFoundException if a required parent value is not given
     */
    @Override
    public Duplicable get(@Nonnull final Map<String, Duplicable> parentsValues)
            throws ComputableNodeFunction.ParentValueNotFoundException, ExternallyComputableNodeValueUnavailableException {
        if (isStoredValueAvailableAndCurrent()) {
            return cachedValue;
        } else if (!isExternallyComputable()) {
            return func.apply(parentsValues); /* may throw {@link ComputableNodeFunction.ParentValueNotFoundException} */
        } else { /* externally computable node */
            throw new ExternallyComputableNodeValueUnavailableException(getKey());
        }
    }

    @Override
    public ComputableCacheNode duplicate() {
        if (isStoredValueAvailable()) {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, true, cachedValue.duplicate(), isCacheCurrent);
        } else {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, cacheEvals, null, isCacheCurrent);
        }
    }

    /**
     * Clones the node with a new value. The new value is not duplicated and is stored by reference.
     *
     * @param newValue the cache value to be replaced with the old value
     * @return a new instance of {@link ComputableCacheNode}
     */
    public ComputableCacheNode duplicateWithUpdatedValue(final Duplicable newValue) {
        if (cacheEvals && newValue != null && !newValue.hasValue()) {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, true, newValue, true);
        } else {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, cacheEvals, null, false);
        }
    }

    /**
     * Returns the computation function
     *
     * @return a function
     */
    public ComputableNodeFunction getFunction() {
        return func;
    }

    /**
     * Duplicates the cache node with outdated cache status. Since we do not ever need to access an
     * outdated cached value, we nullify the cache reference.
     *
     * @return a new instance of {@link ComputableCacheNode}
     */
    public ComputableCacheNode duplicateWithOutdatedCacheStatus() {
        return new ComputableCacheNode(getKey(), getTags(), getParents(), func, cacheEvals, null, false);
    }

    /**
     * This exception will be thrown if a computable function can not be computed
     */
    static final class ExternallyComputableNodeValueUnavailableException extends RuntimeException
            implements Serializable {
        private static final long serialVersionUID = 9056196660803073912L;

        private ExternallyComputableNodeValueUnavailableException(final String nodeKey) {
            super(String.format("Either the externally mutable node \"%s\" is not initialized or is outdated",
                    nodeKey));
        }
    }
}