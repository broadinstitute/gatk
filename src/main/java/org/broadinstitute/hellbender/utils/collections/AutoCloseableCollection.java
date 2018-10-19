package org.broadinstitute.hellbender.utils.collections;

import java.util.Collection;

/**
 * An {@link AutoCloseable} collection that will automatically close all of its elements.
 */
public class AutoCloseableCollection implements AutoCloseable {
    private Collection<? extends AutoCloseable> autoCloseables;

    public AutoCloseableCollection(final Collection<? extends AutoCloseable> autoCloseables) {
        this.autoCloseables = autoCloseables;
    }

    @Override
    public void close() {
        // throw a RuntimeException due to unsuppressable warning when declaring 'throws Exception': https://bugs.openjdk.java.net/browse/JDK-8155591
        RuntimeException exception = null;
        for (AutoCloseable closeable : autoCloseables) {
            try {
                if (closeable != null) {
                    closeable.close();
                }
            } catch (Exception e) {
                if (exception == null) {
                    exception = new RuntimeException("Problem closing AutoCloseable object(s)");
                }
                exception.addSuppressed(e);
            }
        }
        if (exception != null) {
            throw exception;
        }
    }
}
