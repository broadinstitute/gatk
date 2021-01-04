package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SVDeduplicator<T extends SVCallRecord> {

    private final Function<Collection<T>,T> collapser;
    private final SAMSequenceDictionary dictionary;

    public SVDeduplicator(final Function<Collection<T>,T> collapser, final SAMSequenceDictionary dictionary) {
        this.collapser = Utils.nonNull(collapser);
        this.dictionary = Utils.nonNull(dictionary);
    }

    public List<T> deduplicateSortedItems(final List<T> items) {
        if (items.size() < 2) return items;
        final List<T> result = new ArrayList<>();
        int lastPosition = 0;
        String lastContig = null;
        int i = 0;
        final Set<Integer> deduplicatedItems = new HashSet<>();
        final Collection<Integer> identicalItemIndexes = new ArrayList<>();
        while (i < items.size()) {
            if (deduplicatedItems.contains(i)) continue;
            final T record = items.get(i);
            if (!record.getContigA().equals(lastContig)) {
                lastPosition = 0;
                lastContig = record.getContigA();
            }
            if (record.getPositionA() < lastPosition) {
                throw new IllegalArgumentException("Input is not sorted by position");
            }
            int j = i + 1;
            int lastHit = i;
            identicalItemIndexes.clear();
            while (j < items.size() && record.getPositionA() == items.get(j).getPositionA()) {
                final T other = items.get(j);
                if (itemsAreIdentical(record, other)) {
                    identicalItemIndexes.add(j);
                    lastHit = j;
                }
                j++;
            }
            identicalItemIndexes.add(i);
            result.add(collapser.apply(identicalItemIndexes.stream().map(items::get).collect(Collectors.toList())));
            deduplicatedItems.addAll(identicalItemIndexes);
            i = lastHit + 1;
        }
        return result;
    }

    protected boolean itemsAreIdentical(final SVCallRecord a, final SVCallRecord b) {
        return a.getType().equals(b.getType())
                && a.getLength() == b.getLength()
                && a.getContigA().equals(b.getContigA())
                && a.getPositionA() == b.getPositionA()
                && a.getContigB().equals(b.getContigB())
                && a.getPositionB() == b.getPositionB();
    }
}
