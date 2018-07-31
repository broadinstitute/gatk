package org.broadinstitute.hellbender.tools.copynumber.arguments;

import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyNumberArgumentValidationUtils {
    private CopyNumberArgumentValidationUtils() {}

    /**
     * Validate that the interval-argument collection parameters minimally modify the input intervals.
     */
    public static void validateIntervalArgumentCollection(final IntervalArgumentCollection intervalArgumentCollection) {
        Utils.validateArg(intervalArgumentCollection.getIntervalSetRule() == IntervalSetRule.UNION,
                "Interval set rule must be set to UNION.");
        Utils.validateArg(intervalArgumentCollection.getIntervalExclusionPadding() == 0,
                "Interval exclusion padding must be set to 0.");
        Utils.validateArg(intervalArgumentCollection.getIntervalPadding() == 0,
                "Interval padding must be set to 0.");
        Utils.validateArg(intervalArgumentCollection.getIntervalMergingRule() == IntervalMergingRule.OVERLAPPING_ONLY,
                "Interval merging rule must be set to OVERLAPPING_ONLY.");
    }

    /**
     * Validate that a list of locatables is valid and sorted according to a sequence dictionary and contains no duplicates or overlaps.
     */
    public static <T extends Locatable> void validateIntervals(final List<T> intervals,
                                                               final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonNull(intervals);
        Utils.nonNull(sequenceDictionary);
        Utils.validateArg(intervals.stream().allMatch(i -> IntervalUtils.intervalIsOnDictionaryContig(new SimpleInterval(i), sequenceDictionary)),
                "Records contained at least one interval that did not validate against the sequence dictionary.");
        if (!Ordering.from(IntervalUtils.getDictionaryOrderComparator(sequenceDictionary)).isStrictlyOrdered(intervals)) {
            throw new IllegalArgumentException("Records were not strictly sorted in dictionary order.");
        }
        final OptionalInt failureIndex = IntStream.range(1, intervals.size())
                .filter(i -> IntervalUtils.overlaps(intervals.get(i - 1), intervals.get(i)))
                .findFirst();
        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("Records contain at least two overlapping intervals: %s and %s",
                            intervals.get(index - 1), intervals.get(index)));
        }
    }

    /**
     * Compares two non-null sequence dictionaries using sequence index, name, and length only.
     * Less stringent than {@link SAMSequenceDictionary#isSameDictionary}.
     */
    public static boolean isSameDictionary(final SAMSequenceDictionary dictionary1,
                                           final SAMSequenceDictionary dictionary2) {
        Utils.nonNull(dictionary1);
        Utils.nonNull(dictionary2);
        if (dictionary1 == dictionary2) {
            return true;
        }

        final Iterator<SAMSequenceRecord> dictionary1Sequences = dictionary1.getSequences().iterator();
        for (final SAMSequenceRecord dictionary2Sequence : dictionary2.getSequences()) {
            if (!dictionary1Sequences.hasNext()) {
                return false;
            } else {
                final SAMSequenceRecord dictionary1Sequence = dictionary1Sequences.next();
                if (!isSameSequence(dictionary1Sequence, dictionary2Sequence)) {
                    return false;
                }
            }
        }
        return !dictionary1Sequences.hasNext();
    }

    private static boolean isSameSequence(final SAMSequenceRecord sequence1,
                                          final SAMSequenceRecord sequence2) {
        return sequence1 == sequence2 ||
                !(sequence1 == null || sequence2 == null) &&
                        sequence1.getSequenceIndex() == sequence2.getSequenceIndex() &&
                        sequence1.getSequenceName() == sequence2.getSequenceName() &&       // Compare using == since we intern() the Strings
                        !(sequence1.getSequenceLength() != SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH &&
                                sequence2.getSequenceLength() != SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH &&
                                sequence1.getSequenceLength() != sequence2.getSequenceLength());
    }

    /**
     * Checks equality of the sequence dictionary and intervals contained in an {@code locatableCollection}
     * against those contained in an {@link AnnotatedIntervalCollection} represented by {@code annotatedIntervalsFile}.
     * If the latter is {@code null}, then {@code null} is returned; otherwise,
     * the {@link AnnotatedIntervalCollection} represented by {@code inputAnnotatedIntervalsFile} is returned
     * if the intervals are equal, and an exception is thrown if they are not.
     */
    public static AnnotatedIntervalCollection validateAnnotatedIntervals(final File annotatedIntervalsFile,
                                                                         final AbstractLocatableCollection<?, ?> locatableCollection,
                                                                         final Logger logger) {
        Utils.nonNull(locatableCollection);
        Utils.nonNull(logger);
        if (annotatedIntervalsFile == null) {
            logger.info("No GC-content annotations for intervals found; explicit GC-bias correction will not be performed...");
            return null;
        }
        logger.info("Reading and validating GC-content annotations for intervals...");
        final AnnotatedIntervalCollection annotatedIntervals = new AnnotatedIntervalCollection(annotatedIntervalsFile);
        final SAMSequenceDictionary sequenceDictionary = locatableCollection.getMetadata().getSequenceDictionary();
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(annotatedIntervals.getMetadata().getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in annotated-intervals file does not match the master sequence dictionary.");
        }
        Utils.validateArg(annotatedIntervals.getIntervals().equals(locatableCollection.getIntervals()),
                "Annotated intervals do not match provided intervals.");
        return annotatedIntervals;
    }

    /**
     * Same as {@link #validateAnnotatedIntervals}, except we only require that {@code annotatedIntervalsFile}
     * contains as a subset all the intervals contained in {@code locatableCollection} along with equality of the sequence dictionaries.
     * The corresponding subset of annotated intervals is returned if appropriate.
     */
    public static AnnotatedIntervalCollection validateAnnotatedIntervalsSubset(final File annotatedIntervalsFile,
                                                                               final AbstractLocatableCollection<?, ?> locatableCollection,
                                                                               final Logger logger) {
        Utils.nonNull(locatableCollection);
        Utils.nonNull(logger);
        if (annotatedIntervalsFile == null) {
            logger.info("No GC-content annotations for intervals found; explicit GC-bias correction will not be performed...");
            return null;
        }
        logger.info("Reading and validating GC-content annotations for intervals...");
        IOUtils.canReadFile(annotatedIntervalsFile);
        final AnnotatedIntervalCollection annotatedIntervals = new AnnotatedIntervalCollection(annotatedIntervalsFile);
        final SAMSequenceDictionary sequenceDictionary = locatableCollection.getMetadata().getSequenceDictionary();
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(annotatedIntervals.getMetadata().getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in annotated-intervals file does not match the master sequence dictionary.");
        }
        final Set<SimpleInterval> intervalsSubset = new HashSet<>(locatableCollection.getIntervals());
        Utils.validateArg(annotatedIntervals.getIntervals().containsAll(intervalsSubset),
                "Annotated intervals do not contain all specified intervals.");
        return new AnnotatedIntervalCollection(
                locatableCollection.getMetadata(),
                annotatedIntervals.getRecords().stream()
                        .filter(i -> intervalsSubset.contains(i.getInterval()))
                        .collect(Collectors.toList()));
    }
}
