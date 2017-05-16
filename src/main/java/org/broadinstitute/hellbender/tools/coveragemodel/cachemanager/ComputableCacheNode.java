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
final class ComputableCacheNode extends CacheNode {

    private final ComputableNodeFunction func;
    private Duplicable cachedValue = null;
    private final boolean isCaching;
    private boolean isCacheCurrent;

    /**
     * Public constructor
     *
     * @param key the key of the node
     * @param parents parents of the node
     * @param func a function from a map that (at least) contains parents data to the computed value of this node
     * @param isCaching does it store the value or not
     */
    ComputableCacheNode(@Nonnull final NodeKey key,
                        @Nonnull final Collection<NodeTag> tags,
                        @Nonnull final Collection<NodeKey> parents,
                        @Nullable final ComputableNodeFunction func,
                        final boolean isCaching) {
        super(key, tags, parents);
        this.func = func;
        this.isCaching = isCaching;
        Utils.validateArg(func != null || isCaching, "A computable node with null evaluation function is externally" +
                " mutable and must cache its values");
        isCacheCurrent = false;
    }

    private ComputableCacheNode(@Nonnull final NodeKey key,
                                @Nonnull final Collection<NodeTag> tags,
                                @Nonnull final Collection<NodeKey> parents,
                                @Nullable final ComputableNodeFunction func,
                                final boolean isCaching,
                                final Duplicable cachedValue,
                                final boolean isCacheCurrent) {
        super(key, tags, parents);
        this.func = func;
        this.isCaching = isCaching;
        this.isCacheCurrent = isCacheCurrent;
        this.cachedValue = cachedValue;
    }

    @Override
    boolean isPrimitive() { return false; }

    @Override
    boolean isExternallyComputed() { return func == null; }

    boolean isCaching() { return isCaching; }

    /**
     * @return true if the node is caching, has a non-null {@link Duplicable}, the cache is up to date, and the
     * duplicable has a non-null value stored in it
     */
    @Override
    boolean hasValue() {
        return isCaching && isCacheCurrent && cachedValue != null && cachedValue.hasValue();
    }

    @Override
    void set(@Nullable final Duplicable val) {
        if (isExternallyComputed()) {
            cachedValue = val;
            isCacheCurrent = true;
        } else {
            throw new UnsupportedOperationException("Can not explicitly set the value of a computable cache node with" +
                    " non-null function");
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
    Duplicable get(@Nonnull final Map<NodeKey, Duplicable> parentsValues)
            throws ComputableNodeFunction.ParentValueNotFoundException, ExternallyComputableNodeValueUnavailableException {
        if (hasValue()) {
            return cachedValue;
        } else if (!isExternallyComputed()) {
            return func.apply(parentsValues); /* may throw {@link ComputableNodeFunction.ParentValueNotFoundException} */
        } else { /* externally computable node */
            throw new ExternallyComputableNodeValueUnavailableException(getKey());
        }
    }

    @Override
    ComputableCacheNode duplicate() {
        if (hasValue()) {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, true, cachedValue.duplicate(), isCacheCurrent);
        } else {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, isCaching, null, isCacheCurrent);
        }
    }

    /**
     * Clones the node with a new value. The new value is not duplicated and is stored by reference.
     *
     * @param newValue the cache value to be replaced with the old value
     * @return a new instance of {@link ComputableCacheNode}
     */
    @Override
    ComputableCacheNode duplicateWithUpdatedValue(final Duplicable newValue) {
        if (isCaching && newValue != null && newValue.hasValue()) {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, true, newValue, true);
        } else {
            return new ComputableCacheNode(getKey(), getTags(), getParents(), func, isCaching, null, false);
        }
    }

    /**
     * Returns the computation function
     *
     * @return a function
     */
    ComputableNodeFunction getFunction() {
        return func;
    }

    /**
     * Duplicates the cache node with outdated cache status. Since we do not ever need to access an
     * outdated cached value, we nullify the cache reference.
     *
     * @return a new instance of {@link ComputableCacheNode}
     */
    ComputableCacheNode duplicateWithOutdatedCacheStatus() {
        return new ComputableCacheNode(getKey(), getTags(), getParents(), func, isCaching, null, false);
    }

    /**
     * This exception will be thrown if a computable function can not be computed
     */
    static final class ExternallyComputableNodeValueUnavailableException extends RuntimeException
            implements Serializable {
        private static final long serialVersionUID = 9056196660803073912L;

        private ExternallyComputableNodeValueUnavailableException(final NodeKey nodeKey) {
            super(String.format("Either the externally mutable node \"%s\" is not initialized or is outdated",
                    nodeKey));
        }
    }
}