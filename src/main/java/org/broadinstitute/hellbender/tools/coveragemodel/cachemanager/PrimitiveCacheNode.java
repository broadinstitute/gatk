package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;

/**
 * This class represents a primitive cache node (just stores a value and does not perform any computation).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
final class PrimitiveCacheNode extends CacheNode {

    private Duplicable value = null;

    @Override
    boolean isPrimitive() { return true; }

    @Override
    boolean isExternallyComputed() { return true; }

    @Override
    void set(@Nullable final Duplicable val) {
        value = val;
    }

    PrimitiveCacheNode(@Nonnull final NodeKey key,
                       @Nonnull final Collection<NodeTag> tags,
                       @Nullable final Duplicable val) {
        super(key, tags, Collections.emptyList());
        set(val);
    }

    @Override
    boolean hasValue() {
        return value != null && value.hasValue();
    }

    @Override
    Duplicable get(@Nullable final Map<NodeKey, Duplicable> parentsValues)
            throws PrimitiveValueNotInitializedException {
        if (hasValue()) {
            return value;
        } else {
            throw new PrimitiveValueNotInitializedException(getKey());
        }
    }

    @Override
    PrimitiveCacheNode duplicate() {
        if (hasValue()) {
            return new PrimitiveCacheNode(getKey(), getTags(), value.duplicate());
        } else {
            return new PrimitiveCacheNode(getKey(), getTags(), null);
        }
    }

    @Override
    PrimitiveCacheNode duplicateWithUpdatedValue(final Duplicable newValue) {
        return new PrimitiveCacheNode(getKey(), getTags(), newValue);
    }

    /**
     * This exception will be thrown if a primitive node is not initialized but its value is queried
     */
    static final class PrimitiveValueNotInitializedException extends RuntimeException implements Serializable {
        private static final long serialVersionUID = 6036472510998845566L;

        PrimitiveValueNotInitializedException(final NodeKey nodeKey) {
            super("The primitive cache " + ImmutableComputableGraphUtils.quote(nodeKey.toString()) +
                    " is not initialized yet");
        }
    }
}
