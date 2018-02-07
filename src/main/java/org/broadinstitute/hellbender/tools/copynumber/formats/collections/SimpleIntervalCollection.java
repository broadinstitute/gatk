package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.netflix.servo.util.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a collection of {@link SimpleInterval} associated to a sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SimpleIntervalCollection extends AbstractLocatableCollection<LocatableMetadata, SimpleInterval> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END
     */
    enum SimpleIntervalTableColumn {
        CONTIG,
        START,
        END;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, SimpleInterval> SIMPLE_INTERVAL_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(SimpleIntervalTableColumn.CONTIG);
        final int start = dataLine.getInt(SimpleIntervalTableColumn.START);
        final int end = dataLine.getInt(SimpleIntervalTableColumn.END);
        return new SimpleInterval(contig, start, end);
    };

    private static final BiConsumer<SimpleInterval, DataLine> SIMPLE_INTERVAL_RECORD_TO_DATA_LINE_ENCODER = (simpleInterval, dataLine) ->
            dataLine.append(simpleInterval.getContig())
                    .append(simpleInterval.getStart())
                    .append(simpleInterval.getEnd());

    public SimpleIntervalCollection(final File inputFile) {
        super(inputFile, SimpleIntervalTableColumn.COLUMNS, SIMPLE_INTERVAL_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }

    public SimpleIntervalCollection(final LocatableMetadata metadata,
                                    final List<SimpleInterval> simpleIntervals) {
        super(metadata, simpleIntervals, SimpleIntervalTableColumn.COLUMNS, SIMPLE_INTERVAL_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }

    /**
     * Gets a list of {@link SimpleIntervalCollection}s and a {@link SAMSequenceDictionary}, asserts that (1) each
     * collection comprises intervals of ascending order, and (2) the collections span different genomic intervals,
     * and finally returns an int array that yields the ascending order of the collections:
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
     * @param intervalCollections a list of individually ordered {@link SimpleIntervalCollection}s
     * @param samSequenceDictionary a {@link SAMSequenceDictionary} for comparing the order of intervals
     * @param validateResult if true, the result will be comprehensively validated
     * @return an int array
     */
    public static int[] getSimpleIntervalCollectionSortedOrder(final List<SimpleIntervalCollection> intervalCollections,
                                                               final SAMSequenceDictionary samSequenceDictionary,
                                                               final boolean validateResult) {
        Utils.nonNull(intervalCollections);
        Utils.nonNull(samSequenceDictionary);
        intervalCollections.forEach(Utils::nonNull);

        final int numShards = intervalCollections.size();

        final int[] sortOrder = IntStream.range(0, numShards)
                .mapToObj(shardIndex ->
                        ImmutablePair.of(shardIndex, intervalCollections.get(shardIndex).getRecords().get(0)))
                .sorted((p1, p2) -> IntervalUtils.compareLocatables(p1.getValue(), p2.getValue(), samSequenceDictionary))
                .map(ImmutablePair::getLeft)
                .mapToInt(i -> i)
                .toArray();

        if (validateResult) {
            /* assert that the final result is valid (will fail if any the two requirements mentioned
             * in the javadoc is violated) */
            final List<SimpleInterval> sortedAssembledIntervalList = Arrays.stream(sortOrder)
                    .mapToObj(idx -> intervalCollections.get(idx).getRecords())
                    .flatMap(List::stream)
                    .collect(Collectors.toList());
            final boolean assembledArrayIsSorted = IntStream.range(0, sortedAssembledIntervalList.size() - 1)
                    .allMatch(idx -> IntervalUtils.compareLocatables(
                            sortedAssembledIntervalList.get(idx + 1),
                            sortedAssembledIntervalList.get(idx),
                            samSequenceDictionary) > 0);
            Utils.validateArg(assembledArrayIsSorted,
                    "The provided simple interval collections either span overlapping genomic intervals or " +
                            "are not individually sorted.");
        }

        return sortOrder;
    }
}
