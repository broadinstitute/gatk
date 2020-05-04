package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.exceptions.GATKException;
import java.util.function.Consumer;

/**
 * Reference to an object that perform some action when closed.
 * <p>The object refereed, or <em>subject</em>, may or may not be itself closeable and close on the reference does not necessary mean that
 * the subject would also be closed in it happens to be a closeable one. That would be fully determined by the override close action.
 * @param <T>
 */
public abstract class AutoCloseableReference<T> implements AutoCloseable {

    private T subject;

    public static <T> AutoCloseableReference<T> of(final T subject, final Consumer<T> closeAction) {
        Utils.nonNull(closeAction);
        return new AutoCloseableReference<T>(subject) {
            @Override
            protected void close(T subject) {
                closeAction.accept(subject);
            }
        };
    }

    protected AutoCloseableReference(final T subject) {
        Utils.nonNull(subject);
        this.subject = subject;
    }

    /**
     * Returns the refereed object before close is called. After the close action, it will return {@code null}.
     * @return {@code null} if this reference has already being closed.
     */
    public T get() {
        if (subject == null) {
            throw new GATKException("already closed");
        }
        return subject;
    }

    public void close() {
        if (subject != null) {
            close(subject);
            subject = null;
        }
    }

    public boolean isClosed() {
        return subject == null;
    }

    protected abstract void close(final T subject);
}
