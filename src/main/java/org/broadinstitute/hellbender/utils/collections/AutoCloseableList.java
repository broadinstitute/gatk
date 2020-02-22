package org.broadinstitute.hellbender.utils.collections;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntFunction;

public class AutoCloseableList<E extends AutoCloseable> implements AutoCloseable {

    private final List<E> elements;

    private AutoCloseableList(final List<E> elements) {
        this.elements = elements;
    }

    public E get(final int index) {
        return elements.get(index);
    }

    public int size() {
        return elements.size();
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

    @Override
    public void close() {
        RuntimeException first = null;
        for(final E element : elements) {
            try {
                element.close();
            } catch (final Exception ex) {
                if (first == null) {
                    first = new RuntimeException(ex);
                } else {
                    first.addSuppressed(ex);
                }
            }
        }
        if (first != null) {
            throw first;
        }
    }
}
