package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;

/**
 * The base class for all cache nodes in an {@link ImmutableComputableGraph}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class CacheNode {

    private final String key;

    private final Collection<String> parents;

    private final Collection<String> tags;

    /**
     * Public constructor
     *
     * @param key string identifier of the cache node
     * @param tags the tags associated to this cache node
     * @param parents immediate parents of this cache node
     */
    public CacheNode(@Nonnull final String key,
                     @Nonnull final Collection<String> tags,
                     @Nonnull final Collection<String> parents) {
        this.key = Utils.nonNull(key, "The key of a cache node can not be null");
        this.tags = Collections.unmodifiableCollection(Utils.nonNull(tags, "The tag collection of a cache node can not be null"));
        this.parents = Collections.unmodifiableCollection(Utils.nonNull(parents, "The immediate parents of a cache node can not be null"));
    }

    public abstract Duplicable get(@Nonnull final Map<String, Duplicable> dict);

    public abstract boolean isPrimitive();

    public abstract boolean isStoredValueAvailable();

    public abstract void set(@Nullable final Duplicable val);

    public abstract boolean isExternallyComputable();

    public CacheNode duplicateWithUpdatedValue(final Duplicable newValue)
            throws UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    public CacheNode duplicate()
            throws UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    public String getKey() {
        return key;
    }

    public Collection<String> getParents() {
        return parents;
    }

    public Collection<String> getTags() {
        return tags;
    }

    @Override
    public String toString() {
        return key;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CacheNode cacheNode = (CacheNode) o;

        if (!key.equals(cacheNode.key)) return false;
        if (!parents.equals(cacheNode.parents)) return false;
        return tags.equals(cacheNode.tags);
    }

    @Override
    public int hashCode() {
        return key.hashCode();
    }
}
