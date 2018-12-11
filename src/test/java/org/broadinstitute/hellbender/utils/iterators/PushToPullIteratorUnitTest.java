package org.broadinstitute.hellbender.utils.iterators;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Predicate;

public class PushToPullIteratorUnitTest extends GATKBaseTest {

    @Test
    public void testPassThrough() {
        List<String> input = ImmutableList.of("a", "b", "c");
        PushToPullIterator<String> it = new PushToPullIterator<>(input.iterator(), new PassThroughTransformer<>());
        Assert.assertEquals(Lists.newArrayList((Iterator<String>) it), input);
    }

    @Test
    public void testPassThroughEmpty() {
        List<String> input = ImmutableList.of();
        PushToPullIterator<String> it = new PushToPullIterator<>(input.iterator(), new PassThroughTransformer<>());
        Assert.assertEquals(Lists.newArrayList((Iterator<String>) it), input);
    }

    @Test
    public void testFiltering() {
        List<String> input = ImmutableList.of("a", "B", "c");
        PushToPullIterator<String> it = new PushToPullIterator<>(input.iterator(),
                new FilteringTransformer<>(s -> s.toLowerCase().equals(s)));
        Assert.assertEquals(Lists.newArrayList((Iterator<String>) it), ImmutableList.of("a", "c"));
    }

    @Test
    public void testFilteringEmpty() {
        List<String> input = ImmutableList.of();
        PushToPullIterator<String> it = new PushToPullIterator<>(input.iterator(),
                new FilteringTransformer<>(s -> s.toLowerCase().equals(s)));
        Assert.assertEquals(Lists.newArrayList((Iterator<String>) it), input);
    }

    // A transformer that passes through all items
    private static class PassThroughTransformer<T> implements PushPullTransformer<T> {

        List<T> items = new LinkedList<>();

        @Override
        public void submit(T item) {
            items.add(item);
        }

        @Override
        public boolean hasFinalizedItems() {
            return !items.isEmpty();
        }

        @Override
        public List<T> consumeFinalizedItems() {
            List<T> items = this.items;
            this.items = new LinkedList<>();
            return items;
        }

        @Override
        public void signalEndOfInput() {
            // noop
        }
    }

    // A transformer that filters according to a given predicate
    private static class FilteringTransformer<T> extends PassThroughTransformer<T> {
        private Predicate<T> predicate;

        public FilteringTransformer(Predicate<T> predicate) {
            this.predicate = predicate;
        }

        @Override
        public void submit(T item) {
            if (predicate.test(item)) {
                super.submit(item);
            }
        }
    }
}
