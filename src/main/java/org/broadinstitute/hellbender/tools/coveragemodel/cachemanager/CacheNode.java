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
    /**
     * A key for identifying this cache node
     */
    private final NodeKey key;

    /**
     * The collection of lookup keys for parents of this node (can be empty)
     */
    private final Collection<NodeKey> parents;

    /**
     * The collection of tags associated to this node (can be empty)
     */
    private final Collection<NodeTag> tags;

    /**
     * Public constructor
     *
     * @param key lookup key of the cache node
     * @param tags the tags associated to this cache node
     * @param parents lookup keys of the parents of this cache node
     */
    CacheNode(@Nonnull final NodeKey key,
              @Nonnull final Collection<NodeTag> tags,
              @Nonnull final Collection<NodeKey> parents) {
        this.key = Utils.nonNull(key, "The key of a cache node can not be null");
        this.tags = Collections.unmodifiableCollection(Utils.nonNull(tags, "The tag collection of a cache node can not be null"));
        this.parents = Collections.unmodifiableCollection(Utils.nonNull(parents, "The immediate parents of a cache node can not be null"));
    }

    /**
     * Get the value stored in the node
     *
     * @param parents parent values (as a map from their keys to their values)
     * @return a {@link Duplicable}; possibly by reference
     */
    abstract Duplicable get(@Nonnull final Map<NodeKey, Duplicable> parents);

    /**
     * Set the value of the node
     *
     * @param newValue new value; possibly stored by reference
     * @throws UnsupportedOperationException if the node is automatically computable
     */
    abstract void set(@Nullable final Duplicable newValue) throws UnsupportedOperationException;

    /**
     * Is the node primitive?
     */
    abstract boolean isPrimitive();

    /**
     * Is the node initialized yet?
     */
    abstract boolean hasValue();

    /**
     * Is the node externally computed?
     */
    abstract boolean isExternallyComputed();

    /**
     * Duplicate the node with updated value
     *
     * @param newValue new value; possibly stored by reference
     * @return a new {@link CacheNode} with the same key, parents, and tags but with a new value
     * @throws UnsupportedOperationException if the node is automatically computable
     */
    abstract CacheNode duplicateWithUpdatedValue(final Duplicable newValue) throws UnsupportedOperationException;

    /**
     * Make a deep copy of the node
     *
     * @return a deeply copied instance of {@link CacheNode}
     */
    abstract CacheNode duplicate();

    /**
     * Get the string identifier of the node
     * @return a non-null {@link String}
     */
    final NodeKey getKey() {
        return key;
    }

    /**
     * Get the collection of keys of the parents of this node (can be empty)
     */
    final Collection<NodeKey> getParents() {
        return Collections.unmodifiableCollection(parents);
    }

    /**
     * Get the collection of tags associated to this node (can be empty)
     */
    final Collection<NodeTag> getTags() {
        return Collections.unmodifiableCollection(tags);
    }

    @Override
    public final String toString() {
        return key.toString();
    }

    /**
     * NOTE: equality comparison is done just based on the key
     * @param other another object
     */
    @Override
    public final boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;
        return (key.equals(((CacheNode) other).key));
    }

    /**
     * NOTE: hashcode is generated just based on the key
     */
    @Override
    public final int hashCode() {
        return key.hashCode();
    }

    /**
     * This class represents a node key. It is a wrapper around a String.
     */
    public static class NodeKey {
        private final String key;

        public NodeKey(final String key) {
            this.key = Utils.nonNull(key, "Node key must be non-null");
        }

        @Override
        public String toString() {
            return key;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            NodeKey nodeKey = (NodeKey) o;

            return key.equals(nodeKey.key);
        }

        @Override
        public int hashCode() {
            return key.hashCode();
        }
    }

    /**
     * This class represents a node tag. It is a wrapper around a String.
     */
    public static class NodeTag {
        private final String tag;

        public NodeTag(final String tag) {
            this.tag = Utils.nonNull(tag, "Node tag must be non-null");
        }

        @Override
        public String toString() {
            return tag;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            NodeTag nodeTag = (NodeTag) o;

            return tag.equals(nodeTag.tag);
        }

        @Override
        public int hashCode() {
            return tag.hashCode();
        }
    }
}
