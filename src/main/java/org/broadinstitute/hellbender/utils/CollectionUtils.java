package org.broadinstitute.hellbender.utils;

import java.util.*;
import java.util.function.IntFunction;
import java.util.function.Supplier;

/**
 * Created by valentin on 3/24/17.
 */
public final class CollectionUtils {

    //Prevents instantiation.
    private CollectionUtils() {}

    /**
     * Heap-sorts the elements of an {@link Iterable} into a collection.
     *
     * @param <E> type parameter for the element.
     * @param <C> type parameter for the result collection.
     * @param elements the input iterable with the elements to sort.
     * @param comparator the comparator for the elements in the input iterable.
     * @param resultSupplier the result collection supplier function which takes the total number of elements.
     *
     * @throws IllegalArgumentException if any of the arguments in {@code null} or {@code elements} contains {@code null}s.
     *
     * @return never {@code null}.
     */
    public static <E, C extends Collection<? super E>> C heapSort(final Iterable<? extends E> elements,
                                                          final Comparator<? super E> comparator,
                                                          final IntFunction<C> resultSupplier) {
        Utils.nonNull(comparator, "the input comparator cannot be null");
        Utils.nonNull(resultSupplier, "the input result supplier cannot be null");
        Utils.nonNull(elements);

        final PriorityQueue<E> heap = new PriorityQueue<>(comparator);
        for (final E element : elements) {
            if (element == null)
                throw new IllegalArgumentException("the input iterable cannot contain nulls");
            heap.add(element);
        }
        final C result = Utils.nonNull(resultSupplier.apply(heap.size()), "the result collection supplier cannot return a null");
        for (int i = heap.size(); i > 0; --i) {
            result.add(heap.remove());
        }
        return result;
    }

    /**
     * Heap-sorts the elements of an {@link Iterable} into a collection.
     * <p>
     *     It is guaranteed that the input result supplier won't be call until all
     *     the elements have been read from the input iterable.
     * </p>
     * <p>
     *     Therefore if one want to use the same collection as the iterable input and output,
     *     it only need so return its reference after clearing it content.
     * </p>
     *
     * @param <E> type parameter for the element.
     * @param <C> type parameter for the result collection.
     * @param elements the input iterable with the elements to sort.
     * @param resultSupplier the result collection supplier function which takes the total number of elements.
     *
     * @throws IllegalArgumentException if any of the arguments in {@code null} or {@code elements} contains {@code nulls}.
     *
     * @return never {@code null}.
     */
    public static <E extends Comparable<E>, C extends Collection<? super E>> C heapSort(
            final Iterable<? extends E> elements, final IntFunction<C> resultSupplier) {
        return heapSort(elements, Comparable<E>::compareTo, resultSupplier);
    }

    /**
     * Collects the elements in an iterable into a collection.
     *
     *  <p>
     *     The elements in the {@code source} iterable will be added in the iterating order in the result collection
     *     provided by {@code resultSupplier} using the {@link Collection#add} operation.
     * </p>
     * <p>
     *     This means that for certain collection result types, a subsequent iteration through its elements might in fact
     *     yield a different element sequence than an iteration in the source iterable.
     *     (E.g. {@link Set Sets} would remove duplicated elements, the iterating order might be arbitrary).
     * </p>
     *
     * @param elements the elements to collect.
     * @param resultSupplier the result collection factory
     * @param <E> the element type.
     * @param <C> the result collection type.
     * @return never {@code null}.
     * @throws IllegalArgumentException if the result supplier is {@code null} or if it returns a {@code null}.
     */
    public static <E, C extends Collection<E>> C collect(final Iterable<E> elements, final Supplier<C> resultSupplier) {
        final C result = Utils.nonNull(Utils.nonNull(resultSupplier).get(),
                "the result collection returned by the supplier cannot be null");
        for (final E element : Utils.nonNull(elements)) {
            result.add(element);
        }
        return result;
    }

    /**
     * Returns a collection out of the elements of an iterable in the reverse order.
     * <p>
     *     The elements in the {@code source} iterable will be added in the reverse order in the result collection
     *     provided by {@code resultSupplier} using the {@link Collection#add} operation.
     * </p>
     * <p>
     *     This means that for certain collection result types, a subsequent iteration through its elements might not
     *     be in fact the reverse order of the input iterable depending on the semantics of its {@link Collection#add}
     *     implementation. (E.g. {@link Set Sets} would remove duplicated elements).
     * </p>
     * @param source the element source iterable.
     * @param resultSupplier the result collection supplier function that takes on the size of the input.
     * @param <T> the element type.
     * @param <C> the result collection type.
     * @return never {@code null}.
     */
    public static <T, C extends Collection<? super T>> C reverse(final Iterable<? extends T> source, final IntFunction<C> resultSupplier) {
        final Deque<T> stack = new ArrayDeque<>();
        for (final T element : Utils.nonNull(source)) {
            stack.push(element);
        }
        final C result = Utils.nonNull(Utils.nonNull(resultSupplier,
                "the result supplier cannot be null").apply(stack.size()),
                "the result collection supplied cannot be null");
        while (!stack.isEmpty()) {
            result.add(stack.pop());
        }
        return result;
    }
}
