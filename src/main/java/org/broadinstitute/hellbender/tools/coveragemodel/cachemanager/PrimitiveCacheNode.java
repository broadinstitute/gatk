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
public final class PrimitiveCacheNode extends CacheNode {

    private Duplicable value = null;

    @Override
    public boolean isPrimitive() { return true; }

    @Override
    public boolean isExternallyComputable() { return true; }

    @Override
    public void set(@Nullable final Duplicable val) {
        value = val;
    }

    public PrimitiveCacheNode(@Nonnull final String key,
                              @Nonnull final Collection<String> tags,
                              @Nullable final Duplicable val) {
        super(key, tags, Collections.emptyList());
        set(val);
    }

    @Override
    public boolean isStoredValueAvailable() {
        return value != null && !value.hasValue();
    }

    @Override
    public Duplicable get(@Nullable final Map<String, Duplicable> parentsValues)
            throws PrimitiveValueNotInitializedException {
        if (isStoredValueAvailable()) {
            return value;
        } else {
            throw new PrimitiveValueNotInitializedException(String.format(
                    "The primitive cache \"%s\" is not initialized yet", getKey()));
        }
    }

    @Override
    public PrimitiveCacheNode duplicate() {
        if (value != null && !value.hasValue()) {
            return new PrimitiveCacheNode(getKey(), getTags(), value.duplicate());
        } else {
            return new PrimitiveCacheNode(getKey(), getTags(), null);
        }
    }

    @Override
    public PrimitiveCacheNode duplicateWithUpdatedValue(final Duplicable newValue) {
        return new PrimitiveCacheNode(getKey(), getTags(), newValue);
    }

    /**
     * This exception will be thrown if a primitive node is not initialized but its value is queried
     */
    static final class PrimitiveValueNotInitializedException extends RuntimeException implements Serializable {
        private static final long serialVersionUID = 6036472510998845566L;

        private PrimitiveValueNotInitializedException(String s) {
            super(s);
        }
    }
}
