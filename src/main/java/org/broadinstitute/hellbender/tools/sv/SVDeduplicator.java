package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SVDeduplicator<T extends SVLocatable> {

    final Function<Collection<T>,T> collapser;
    final SAMSequenceDictionary dictionary;

    public SVDeduplicator(final Function<Collection<T>,T> collapser, final SAMSequenceDictionary dictionary) {
        this.collapser = Utils.nonNull(collapser);
        this.dictionary = Utils.nonNull(dictionary);
    }

    public List<T> deduplicateItems(final List<T> items) {
        Utils.nonNull(dictionary);
        final List<T> sortedItems = items.stream().sorted(SVCallRecordUtils.getSVLocatableComparator(dictionary)).collect(Collectors.toList());
        final List<T> deduplicatedList = new ArrayList<>();
        int i = 0;
        while (i < sortedItems.size()) {
            final T record = sortedItems.get(i);
            int j = i + 1;
            final Collection<Integer> identicalItemIndexes = new ArrayList<>();
            while (j < sortedItems.size() && record.getPositionA() == sortedItems.get(j).getPositionA()) {
                final T other = sortedItems.get(j);
                if (itemsAreIdentical(record, other)) {
                    identicalItemIndexes.add(j);
                }
                j++;
            }
            if (identicalItemIndexes.isEmpty()) {
                deduplicatedList.add(record);
                i++;
            } else {
                identicalItemIndexes.add(i);
                final List<T> identicalItems = identicalItemIndexes.stream().map(sortedItems::get).collect(Collectors.toList());
                deduplicatedList.add(collapser.apply(identicalItems));
                i = j;
            }
        }
        return deduplicatedList;
    }

    public boolean itemsAreIdentical(final T a, final T b) {
        return a.getContigA().equals(b.getContigA())
                && a.getPositionA() == b.getPositionA()
                && a.getContigB().equals(b.getContigB())
                && a.getPositionB() == b.getPositionB();
    }
}
