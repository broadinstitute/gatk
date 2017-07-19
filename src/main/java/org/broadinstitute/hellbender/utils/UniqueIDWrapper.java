package org.broadinstitute.hellbender.utils;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Create a unique ID for an arbitrary object and wrap it.
 *
 * Warning:  This will not work as desired on a Spark cluster (local spark should be fine)
 */
public class UniqueIDWrapper<A> {
    private static final AtomicLong counter = new AtomicLong();

    private final long id;
    private final A wrapped;

    public UniqueIDWrapper(A toWrap) {
        id = counter.getAndIncrement();
        wrapped = toWrap;
    }

    public long getId() {
        return id;
    }

    public A getWrapped() {
        return wrapped;
    }
}
