package org.broadinstitute.hellbender.utils.collections;

import java.util.Collection;
import java.util.Objects;
import java.util.Set;
import java.util.function.BiPredicate;

/**
 * Implements Set by imposing a unique-element constraint on HopscotchCollection.
 * Also implements equals and hashCode to be consistent with the documented requirements of the java Set interface.
 */
public class HopscotchSet<T> extends HopscotchCollection<T> implements Set<T> {
    public HopscotchSet() {}
    public HopscotchSet( final int capacity ) { super(capacity); }
    public HopscotchSet( final Collection<? extends T> collection ) { super(collection); }

    @Override
    public final boolean equals( final Object obj ) {
        if ( this == obj ) return true;
        if ( !(obj instanceof Set) ) return false;
        @SuppressWarnings("rawtypes") final Set that = (Set)obj;
        return this.size() == that.size() && this.containsAll(that);
    }

    @Override
    public final int hashCode() { return stream().mapToInt(Objects::hashCode).sum(); }

    /**
     * in a set, uniqueness is on the entry, so we just compare entries to see if we already have an identical one
     */
    @Override
    protected BiPredicate<T, T> entryCollides() { return Object::equals; }
}
