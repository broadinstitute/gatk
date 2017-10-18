package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.Ordering;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.OptionalInt;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents an immutable, coordinate-sorted (with no overlaps allowed) collection of records
 * that extend {@link Locatable} (although contigs are assumed to be non-null when writing to file),
 * a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records. Records are sorted using
 * {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}. See the unit test
 * for a simple example of a subclass.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class LocatableCollection<RECORD extends Locatable> extends RecordCollection<RECORD> {
    public static final Comparator<Locatable> LEXICOGRAPHICAL_ORDER_COMPARATOR = IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

    /**
     * Records are sorted using {@code LEXICOGRAPHICAL_ORDER_COMPARATOR}.
     */
    public LocatableCollection(final List<RECORD> records,
                               final TableColumnCollection mandatoryColumns,
                               final Function<DataLine, RECORD> recordFromDataLineDecoder,
                               final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(
                Utils.nonNull(records).stream().sorted(LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList()),
                mandatoryColumns,
                recordFromDataLineDecoder,
                recordToDataLineEncoder);
        validateIntervals(getRecords());
    }

    /**
     * @throws IllegalArgumentException if records are not sorted using {@code LEXICOGRAPHICAL_ORDER_COMPARATOR}
     */
    public LocatableCollection(final File inputFile,
                               final TableColumnCollection mandatoryColumns,
                               final Function<DataLine, RECORD> recordFromDataLineDecoder,
                               final BiConsumer<RECORD, DataLine> recordToDataLineEncoder) {
        super(inputFile, mandatoryColumns, recordFromDataLineDecoder, recordToDataLineEncoder);
        validateIntervals(getRecords());
    }

    //check for ordering, duplicates, and overlaps
    private static <T extends Locatable> void validateIntervals(final List<T> records) {
        if (!Ordering.from(LEXICOGRAPHICAL_ORDER_COMPARATOR).isStrictlyOrdered(records)) {
            throw new IllegalArgumentException("Records fwere not strictly sorted in lexicographical order.");
        }
        final OptionalInt failureIndex = IntStream.range(1, records.size())
                .filter(i -> IntervalUtils.overlaps(records.get(i - 1), records.get(i)))
                .findFirst();
        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("Records contain at least two overlapping intervals: %s and %s",
                            records.get(index - 1), records.get(index)));
        }
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

    public OverlapDetector<RECORD> getOverlapDetector() {
        return OverlapDetector.create(getRecords());
    }
}