package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.List;

/**
 * A class that receives a stream of elements and transforms or filters them in some way, such as by downsampling with
 * a {@link org.broadinstitute.hellbender.utils.downsampling.Downsampler}. Elements are submitted in a push-style model,
 * in contrast to Java's pull-style {@link java.util.Iterator}. A transformer may be used to transform an iterator of
 * elements using {@link PushToPullIterator}.
 *
 * @param <T> type of items to be submitted
 * @see PushToPullIterator
 * @see org.broadinstitute.hellbender.utils.downsampling.Downsampler
 */
public interface PushPullTransformer<T> {
    /**
     * Submit one item to the transformer for consideration. Some transformers will be able to determine
     * immediately whether the item survives the transformation process, while others will need to see
     * more items before making that determination.
     *
     * @param item the individual item to submit to the transformer for consideration
     */
    void submit(T item);

    /**
     * Are there items that have survived the transformation process waiting to be retrieved?
     *
     * @return true if this transformer has > 0 finalized items, otherwise false
     */
    boolean hasFinalizedItems();

    /**
     * Return (and *remove*) all items that have survived transformation and are waiting to be retrieved.
     *
     * @return a list of all finalized items this transformer contains, or an empty list if there are none
     */
    List<T> consumeFinalizedItems();

    /**
     * Used to tell the transformer that no more items will be submitted to it, and that it should
     * finalize any pending items.
     */
    void signalEndOfInput();

    default void submit(final Collection<T> items) {
        Utils.nonNull(items, "submitted items must not be null");
        items.forEach(this::submit);
    }
}
