package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.List;

/**
 *
 * @param <T>
 */
public interface PushPullTransformer<T> {
    /**
     * Submit one item to the downsampler for consideration. Some downsamplers will be able to determine
     * immediately whether the item survives the downsampling process, while others will need to see
     * more items before making that determination.
     *
     * @param item the individual item to submit to the downsampler for consideration
     */
    void submit(T item);

    /**
     * Are there items that have survived the downsampling process waiting to be retrieved?
     *
     * @return true if this downsampler has > 0 finalized items, otherwise false
     */
    boolean hasFinalizedItems();

    /**
     * Return (and *remove*) all items that have survived downsampling and are waiting to be retrieved.
     *
     * @return a list of all finalized items this downsampler contains, or an empty list if there are none
     */
    List<T> consumeFinalizedItems();

    /**
     * Used to tell the downsampler that no more items will be submitted to it, and that it should
     * finalize any pending items.
     */
    void signalEndOfInput();

    default void submit(final Collection<T> items) {
        Utils.nonNull(items, "submitted items must not be null");
        items.forEach(this::submit);
    }
}
