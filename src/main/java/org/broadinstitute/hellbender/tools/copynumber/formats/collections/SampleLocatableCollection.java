package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.Ordering;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
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
 * that extend {@link Locatable} (although contigs are assumed to be non-null when writing to file)
 * associated with a sample, a set of mandatory column headers given by a {@link TableColumnCollection},
 * and lambdas for reading and writing records. Records are sorted using
 * {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}. See the unit test
 * for a simple example of a subclass.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class SampleLocatableCollection<T extends Locatable> extends SampleRecordCollection<T> {
    public static final Comparator<Locatable> LEXICOGRAPHICAL_ORDER_COMPARATOR = IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

    /**
     * Records are sorted using {@code LEXICOGRAPHICAL_ORDER_COMPARATOR}.
     */
    public SampleLocatableCollection(final SampleMetadata sampleMetadata,
                                     final List<T> records,
                                     final TableColumnCollection mandatoryColumns,
                                     final Function<DataLine, T> dataLineToRecordFunction,
                                     final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        super(
                sampleMetadata,
                Utils.nonNull(records).stream().sorted(LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList()),
                mandatoryColumns,
                dataLineToRecordFunction,
                recordAndDataLineBiConsumer);
        validateIntervals(getSampleName(), getIntervals());
    }

    /**
     * @throws IllegalArgumentException if records are not sorted using {@code LEXICOGRAPHICAL_ORDER_COMPARATOR}
     */
    public SampleLocatableCollection(final File inputFile,
                                     final TableColumnCollection mandatoryColumns,
                                     final Function<DataLine, T> dataLineToRecordFunction,
                                     final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        super(inputFile, mandatoryColumns, dataLineToRecordFunction, recordAndDataLineBiConsumer);
        validateIntervals(getSampleName(), getIntervals());
    }

    //check for ordering, duplicates, and overlaps
    private static void validateIntervals(final String sampleName,
                                          final List<SimpleInterval> intervals) {
        if (!Ordering.from(LEXICOGRAPHICAL_ORDER_COMPARATOR).isStrictlyOrdered(intervals)) {
            throw new IllegalArgumentException(String.format("Records for sample %s were not strictly sorted in lexicographical order.", sampleName));
        }
        final OptionalInt failureIndex = IntStream.range(1, intervals.size())
                .filter(i -> IntervalUtils.overlaps(intervals.get(i - 1), intervals.get(i)))
                .findFirst();
        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("Records for sample %s contain at least two overlapping intervals: %s and %s",
                            sampleName, intervals.get(index - 1), intervals.get(index)));
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

    public OverlapDetector<T> getOverlapDetector() {
        return OverlapDetector.create(getRecords());
    }
}