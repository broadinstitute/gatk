package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * This class represents a primitive cache node (just stores a value and does not perform any computation)
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class PrimitiveCacheNode extends CacheNode {

    private static List<String> EMPTY_LIST = new ArrayList<>();

    private Duplicable value = null;

    @Override
    public boolean isPrimitive() { return true; }

    @Override
    public boolean isExternallyMutable() { return true; }

    @Override
    public void setValue(@Nullable final Duplicable val) {
        value = val;
    }

    public PrimitiveCacheNode(@Nonnull final String key,
                              @Nonnull final Collection<String> tags,
                              @Nullable final Duplicable val) {
        super(key, tags, EMPTY_LIST);
        setValue(val);
    }

    @Override
    public boolean isStoredValueAvailable() {
        return value != null && !value.isNull();
    }

    @Override
    public Duplicable getValue(@Nullable final Map<String, ? extends Duplicable> parentsValues) {
        if (value == null || value.isNull()) {
            throw new IllegalStateException("The value for primitive cache (" + toString() + ") is not set yet.");
        } else {
            return value;
        }
    }

    @Override
    public PrimitiveCacheNode duplicate() {
        if (value != null && !value.isNull()) {
            return new PrimitiveCacheNode(getKey(), getTags(), value.deepCopy());
        } else {
            return new PrimitiveCacheNode(getKey(), getTags(), null);
        }
    }

    @Override
    public PrimitiveCacheNode duplicateWithUpdatedValue(final Duplicable newValue) {
        return new PrimitiveCacheNode(getKey(), getTags(), newValue);
    }
}
