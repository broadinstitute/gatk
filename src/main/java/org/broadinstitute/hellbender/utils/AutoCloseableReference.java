package org.broadinstitute.hellbender.utils;

import java.util.function.Consumer;

/**
 * Reference to another object that perform some action when closed.
 *
 * <p>The object referred, the <em>subject</em>, may or may not be itself closeable and close invoke on this reference
 * does not necessary mean that the subject would also be closed; that depends on whether this reference closing action
 * actually propagates it to the referred object.
 * </p>
 * <p>This is useful for example in the following sceneraios:
 *  <ul>
 *      <li>
 *          the referred object does not or cannot extends/implement {@link AutoCloseable} but we want to use
 *          in <i>try-with-resources</i>.
 *      </li>
 *      <li>
 *          the referred object is auto-closeable but its close action is not the one we want to trigger
 *          at the end of a <i>try-with-resource</i>; is a totally different one or needs to be augmented.
 *      </li>
 *  </ul>
 * </p>
 * <p>
 *     This class is not thread-safe.
 * </p>
 * @param <T> the subject's type.
 */
public abstract class AutoCloseableReference<T> implements AutoCloseable {

    /**
     * Holds the reference to the subject.
     * <p>
     *     Becomes {@code null} once the the close action has been invoked.
     * </p>
     */
    private T subject;

    /**
     * Creates a new reference given the subject object an the close action.
     * @param subject the subject/referred object.
     * @param closeAction the close action.
     * @param <T> type parameter for the {@code subject}.
     * @return never {@code null}.
     * @throws IllegalArgumentException if either subject or close-action is {@code null}.
     */
    public static <T> AutoCloseableReference<T> of(final T subject, final Consumer<T> closeAction) {
        Utils.nonNull(closeAction);
        return new AutoCloseableReference<T>(subject) {
            @Override
            protected void close(T subject) {
                closeAction.accept(subject);
            }
        };
    }

    /**
     * Extending classes need to pass the subject reference using this constructor.
     * @param subject the subject reference.
     * @throws IllegalArgumentException if the subject provided is null.
     */
    protected AutoCloseableReference(final T subject) {
        this.subject = Utils.nonNull(subject);
    }

    /**
     * Returns the referred object.
     * <p>Trying to retrieve the subject after wrapper has been closed results in an exception}</p>
     * @return never {@code null}.
     */
    public final T get() {
        if (subject == null) {
            throw new IllegalStateException("reference already closed");
        }
        return subject;
    }

    @Override
    public final void close() {
        final T subject = this.subject;
        if (subject != null) {
            this.subject = null; // this way we avoid calling close twice if there is
                                 // some issues within the action.
            close(subject);
        }
    }

    /**
     * Checks whether this reference has already been closed.
     * @return true iff closed was invoked on it.
     */
    public final boolean isClosed() {
        return subject == null;
    }

    /**
     * Method that performs additional closing actions given this
     * wrapper subject.
     *
     * @param subject this wrapper subject.
     */
    protected abstract void close(final T subject);
}
