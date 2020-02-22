package org.broadinstitute.hellbender.utils.collections;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntFunction;

/**
 * List of autocloseables that when closed it closes all its elements.
 * @param <E> the element auto-closeable type.
 */
public final class AutoCloseableList<E extends AutoCloseable> extends AutoCloseableCollection<List<E>> {

    private AutoCloseableList(final List<E> elements) {
        super(elements);
    }

    public E get(final int index) {
        return autoCloseables.get(index);
    }

    public static <E extends AutoCloseable> AutoCloseableList<E> of(final List<E> elements) {
        final List<E> clone = new ArrayList<>(elements);
        return new AutoCloseableList<>(clone);
    }

    public static  <E extends AutoCloseable> AutoCloseableList<E> of(final int size, final IntFunction<E> constr) {
        final List<E> elements = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            elements.add(constr.apply(i));
        }
        return new AutoCloseableList<>(elements);
    }
}
