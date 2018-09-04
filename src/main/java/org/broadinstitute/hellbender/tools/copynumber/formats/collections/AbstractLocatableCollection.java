package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import javax.annotation.Nonnull;
import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a sequence dictionary, an immutable, coordinate-sorted (with no overlaps allowed) collection of records
 * that extend {@link Locatable} (although contigs are assumed to be non-null when writing to file),
 * a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records. Records are sorted using the sequence dictionary.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public abstract class AbstractLocatableCollection<METADATA extends LocatableMetadata, RECORD extends Locatable> extends AbstractRecordCollection<METADATA, RECORD> {
    private final Lazy<OverlapDetector<RECORD>> overlapDetector;

    /**
     * @param metadata records are sorted using the contained {@link SAMSequenceDictionary}
     */
    AbstractLocatableCollection(final METADATA metadata,
                                final List<RECORD> records,
                                final TableColumnCollection mandatoryColumns,
                                final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(metadata, sortRecords(records, metadata.getSequenceDictionary()), mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
        CopyNumberArgumentValidationUtils.validateIntervals(getRecords(), metadata.getSequenceDictionary());
        this.overlapDetector = new Lazy<>(() -> OverlapDetector.create(getRecords()));
    }

    /**
     * @throws IllegalArgumentException if records are not sorted according to the {@link SAMSequenceDictionary} contained in the input file
     */
    AbstractLocatableCollection(final File inputFile,
                                final TableColumnCollection mandatoryColumns,
                                final Function<DataLine, RECORD> recordFromDataLineDecoder,
                                final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
        CopyNumberArgumentValidationUtils.validateIntervals(getRecords(), getMetadata().getSequenceDictionary());
        this.overlapDetector = new Lazy<>(() -> OverlapDetector.create(getRecords()));
    }

    private static <T extends Locatable> List<T> sortRecords(final List<T> records,
                                                             final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonNull(records);
        Utils.nonNull(sequenceDictionary);
        return records.stream()
                .sorted(IntervalUtils.getDictionaryOrderComparator(sequenceDictionary))
                .collect(Collectors.toList());
    }

    /**
     * @return  a new modifiable list of {@link SimpleInterval}s corresponding to the {@link Locatable}s
     *          for each record contained in the collection
     */
    public List<SimpleInterval> getIntervals() {
        return getRecords().stream()
                .map(r -> new SimpleInterval(r.getContig(), r.getStart(), r.getEnd()))
                .collect(Collectors.toList());
    }

    /**
     * @return a lazily-created {@link OverlapDetector} of the contained records
     */
    public OverlapDetector<RECORD> getOverlapDetector() {
        return overlapDetector.get();
    }

    public Comparator<Locatable> getComparator() {
        return IntervalUtils.getDictionaryOrderComparator(getMetadata().getSequenceDictionary());
    }

    @Override
    protected Metadata.Type getMetadataType() {
        return Metadata.Type.LOCATABLE;
    }

    /**
     * Gets a list of locatable collections and returns the ascending sort order of collections.
     * It is also asserted that the collections have the same metadata, each collection is composed of sorted
     * intervals, and the collections span non-overlapping regions.
     *
     *  Example:
     *
     *  (ascending order -->)
     *
     *  collection 0:                       |----|  |--|  |-----| |-|
     *  collection 1: |--| |----|   \----|
     *  collection 2:                                                 |----| |--| |-|
     *
     *  output: [1, 0, 2]
     *
     * @param collections a list of collections
     * @return sort order as a list of integers
     */
    public static <METADATA extends LocatableMetadata, RECORD extends Locatable> List<Integer> getShardedCollectionSortOrder(
            @Nonnull final List<? extends AbstractLocatableCollection<METADATA, RECORD>> collections) {
        Utils.nonNull(collections);
        collections.forEach(Utils::nonNull);
        Utils.validateArg(!collections.isEmpty(), "The list must contain at least one collection.");
        final METADATA firstMetadata = collections.get(0).getMetadata();
        Utils.validateArg(collections.stream()
                .map(AbstractLocatableCollection<METADATA, RECORD>::getMetadata)
                        .allMatch(metadata -> metadata.equals(firstMetadata)),
                "The collections must have the same metadata.");

        final int numShards = collections.size();
        final List<Integer> sortOrder = IntStream.range(0, numShards)
                .mapToObj(shardIndex ->
                        ImmutablePair.of(shardIndex, collections.get(shardIndex).getRecords().get(0)))
                .sorted((p1, p2) -> IntervalUtils.compareLocatables(p1.getValue(), p2.getValue(),
                        firstMetadata.getSequenceDictionary()))
                .map(ImmutablePair::getLeft)
                .collect(Collectors.toList());

        /* assert that the final result is valid (will fail if any the two requirements mentioned
         * in the javadoc is violated) */
        final List<RECORD> sortedAssembledRecordList = sortOrder.stream()
                .map(idx -> collections.get(idx).getRecords())
                .flatMap(List::stream)
                .collect(Collectors.toList());
        CopyNumberArgumentValidationUtils.validateIntervals(sortedAssembledRecordList,
                firstMetadata.getSequenceDictionary());

        return sortOrder;
    }
}