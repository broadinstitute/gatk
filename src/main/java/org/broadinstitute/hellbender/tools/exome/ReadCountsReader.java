package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Reader for a read counts table.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadCountsReader extends TableReader<ReadCountRecord> {

    private final Function<DataLine, ReadCountRecord> recordExtractor;

    private final List<String> countColumnNames;

    private final TargetCollection<Target> targets;

    private final boolean ignoreMissingTargets;

    public ReadCountsReader(final String sourceName, final Reader sourceReader,
                            final TargetCollection<Target> targets,
                            final boolean ignoreMissingTargets) throws IOException {
        super(sourceName, sourceReader);
        this.targets = targets;
        this.ignoreMissingTargets = ignoreMissingTargets;
        if (targets == null && ignoreMissingTargets) {
            throw new IllegalArgumentException("When ignore missing targets is true, targets cannot be null");
        }
        countColumnNames = super.columns().names().stream()
                .filter(name -> !TargetTableColumn.isStandardTargetColumnName(name))
                .collect(Collectors.toList());
        final Function<DataLine, SimpleInterval> intervalExtractor = intervalExtractor(super.columns(),
                (message) -> formatException(message));
        final Function<DataLine, String> targetNameExtractor = targetNameExtractor(super.columns());
        final Function<DataLine, double[]> countExtractor = countExtractor(super.columns());
        recordExtractor = composeRecordExtractor(intervalExtractor, targetNameExtractor, countExtractor, targets);
    }

    /**
     * Returns the list with the count column names.
     * @return unmodifiable list with the count column names in the same order as they appear in the records.
     */
    public List<String> getCountColumnNames() {
        return Collections.unmodifiableList(countColumnNames);
    }

    public ReadCountsReader(final File file, final TargetCollection<Target> targets, final boolean ignoreMissingTargets)
            throws IOException {
        super(file);
        this.targets = targets;
        this.ignoreMissingTargets = ignoreMissingTargets;
        if (targets == null && ignoreMissingTargets) {
            throw new IllegalArgumentException("When ignore missing targets is true, targets cannot be null");
        }
        countColumnNames = super.columns().names().stream()
                .filter(name -> !TargetTableColumn.isStandardTargetColumnName(name))
                .collect(Collectors.toList());
        final Function<DataLine, SimpleInterval> intervalExtractor = intervalExtractor(super.columns(),
                (message) -> formatException(message));
        final Function<DataLine, String> targetNameExtractor = targetNameExtractor(super.columns());
        final Function<DataLine, double[]> countExtractor = countExtractor(super.columns());
        recordExtractor = composeRecordExtractor(intervalExtractor, targetNameExtractor, countExtractor, targets);
    }

    public ReadCountsReader(final Reader sourceReader, final TargetCollection<Target> targets,
                            final boolean ignoreMissingTargets) throws IOException {
        super(sourceReader);
        this.targets = targets;
        this.ignoreMissingTargets = ignoreMissingTargets;
        if (targets == null && ignoreMissingTargets) {
            throw new IllegalArgumentException("When ignore missing targets is true, targets cannot be null");
        }
        countColumnNames = super.columns().names().stream()
                .filter(name -> !TargetTableColumn.isStandardTargetColumnName(name))
                .collect(Collectors.toList());
        final Function<DataLine, SimpleInterval> intervalExtractor = intervalExtractor(super.columns(),
                (message) -> formatException(message));
        final Function<DataLine, String> targetNameExtractor = targetNameExtractor(super.columns());
        final Function<DataLine, double[]> countExtractor = countExtractor(super.columns());
        recordExtractor = composeRecordExtractor(intervalExtractor, targetNameExtractor, countExtractor, targets);
    }

    protected ReadCountsReader(final String sourceName, final Reader sourceReader) throws IOException {
        this(sourceName, sourceReader, null, false);
    }

    public ReadCountsReader(final File file) throws IOException {
        this(file, null, false);
    }

    public ReadCountsReader(final Reader sourceReader) throws IOException {
        this(sourceReader, null, false);
    }

    /**
     * Returns a function that translate an source line string value array into
     * into a interval.
     *
     * @param contigColumnNumber    the number of the input column that contains the
     *                              contig name. {@code -1} if missing.
     * @param startColumnNumber     the number of the input column that contains the
     *                              start position. {@code -1} if missing.
     * @param endColumnNumber       the number of the input column that contains the
     *                              end position. {@code -1} if missing.
     * @param errorExceptionFactory instantiates the exception to thrown in case
     *                              of a formatting error. cannot be {@code null}.
     * @return not a {@code null} if there is enough information to find out the
     * sample intervals, {@code null} if it is insufficient.
     */
    private static Function<DataLine, SimpleInterval> composeIntervalBuilder(
            final int contigColumnNumber,final int startColumnNumber, final int endColumnNumber,
            final Function<String, RuntimeException> errorExceptionFactory) {
        if (contigColumnNumber == -1 || startColumnNumber == -1 || endColumnNumber == -1) {
            return null;
        }

        return (v) -> {
            final String contig = v.get(contigColumnNumber);
            final int start = v.getInt(startColumnNumber);
            final int end = v.getInt(endColumnNumber);
            if (start <= 0) {
                throw errorExceptionFactory.apply(String.format("start position must be greater than 0: %d", start));
            } else if (start > end) {
                throw errorExceptionFactory.apply(String.format("end position '%d' must equal or greater than the start position '%d'", end, start));
            } else {
                return new SimpleInterval(contig, start, end);
            }
        };
    }

    /**
     * Composes a lambda to extract the name of the target given a row of values from the input read-count file.
     * <p>
     * This method will return {@code null} if it is not possible to extract the target name from the input directly; for
     * example the input only contain the coordinates of the target and not the target name itself (
     * (i.e. the {@link TargetTableColumn#NAME NAME} column is missing).
     * </p>
     *
     * @param columns the column-name array for that file.
     * @return non-{@code null} iff is not possible to extract the target name from the input directly.
     */
    private static Function<DataLine, String> targetNameExtractor(final TableColumnCollection columns) {
        final int nameColumnIndex = columns.indexOf(TargetTableColumn.NAME.toString());
        return nameColumnIndex < 0 ? null : (v) -> v.get(nameColumnIndex);
    }

    /**
     * Constructs an per line interval extractor given the header column names.
     *
     * @param columns               the header column names.
     * @param errorExceptionFactory the error handler to be called when there is any problem resoling the interval.
     * @return never {@code null} if there is enough columns to extract the coordinate information, {@code null} otherwise.
     */
    private static Function<DataLine, SimpleInterval> intervalExtractor(
            final TableColumnCollection columns, final Function<String, RuntimeException> errorExceptionFactory) {

        final int contigColumnNumber = columns.indexOf(TargetTableColumn.CONTIG.toString());
        final int startColumnNumber = columns.indexOf(TargetTableColumn.START.toString());
        final int endColumnNumber = columns.indexOf(TargetTableColumn.END.toString());
        return composeIntervalBuilder(contigColumnNumber, startColumnNumber, endColumnNumber, errorExceptionFactory);
    }


    private <T extends Target> Function<DataLine, ReadCountRecord> composeRecordExtractor(
            final Function<DataLine, SimpleInterval> intervalExtractor,
            final Function<DataLine, String> targetNameExtractor, final Function<DataLine, double[]> countExtractor,
            final TargetCollection<?> targets) {
        if (targetNameExtractor == null && (targets == null || intervalExtractor == null)) {
            throw formatException("the input files does not contain a target name column and no target file was provided");
        } else if (targetNameExtractor != null && intervalExtractor != null) {
            if (targets == null) {
                return (v) -> new ReadCountRecord(new Target(targetNameExtractor.apply(v), intervalExtractor.apply(v)), countExtractor.apply(v));
            } else {
                return (v) -> extractRecordWithTargetNameAndIntervalExtractorsAndTargetCollection(v, intervalExtractor, targetNameExtractor, countExtractor);
            }
        } else if (targetNameExtractor != null) {
            if (targets == null) {
                return (v) -> new ReadCountRecord(new Target(targetNameExtractor.apply(v), null), countExtractor.apply(v));
            } else {
                return (v) -> extractRecordWithTargetNameExtractorAndTargetCollection(v, targetNameExtractor, countExtractor);
            }
        } else { // at this point targets != null && intervalExtractor != null.
            return (v) -> extractRecordWithIntervalExtractorAndTargetCollection(v, intervalExtractor, countExtractor);
        }
    }

    private ReadCountRecord extractRecordWithTargetNameExtractorAndTargetCollection(
            final DataLine v, final Function<DataLine, String> targetNameExtractor,
            final Function<DataLine, double[]> countExtractor) {
        final String name = targetNameExtractor.apply(v);
        final Target target = createRecordTarget(targets, name);
        if (target == null) {
            if (ignoreMissingTargets) {
                return null;
            } else {
                throw formatException(String.format("unknown target '%s' not present in the target collection", name));
            }
        } else {
            return new ReadCountRecord(target, countExtractor.apply(v));
        }
    }

    /**
     * Creates a record target using the targets collection and location.
     * @param targets targets collection.
     * @param where query interval.
     * @throws UserException.BadInput if the overlapped target's location is not exactly the one provided in the input.
     * @return {@code null} if there is no target at {@code where}.
     */
    private Target createRecordTarget(final TargetCollection<Target> targets, final SimpleInterval where) {
        final Target target = targets.target(where);
        if (target == null) {
            return null;
        } else if (target.getInterval().equals(where)) {
            return new Target(target.getName(), target.getInterval());
        } else {
            throw formatException(String.format("mismatching yet overlapping intervals in the input (%s) and the target collection (%s)", where, targets.location(target)));
        }
    }

    /**
     * Creates a record target using the targets collection and location.
     * @param targets targets collection.
     * @param name query interval.
     * @param <E> the target type.
     * @throws UserException.BadInput if the overlapped target's location is not exactly the one provided in the input.
     * @return {@code null} if there is no target at {@code where}.
     */
    private <E> Target createRecordTarget(final TargetCollection<E> targets, final String name) {
        final E target = targets.target(name);
        if (target == null) {
            return null;
        } else {
            return new Target(name, targets.location(target));
        }
    }

    private ReadCountRecord extractRecordWithIntervalExtractorAndTargetCollection(
            final DataLine v, final Function<DataLine, SimpleInterval> intervalExtractor,
            final Function<DataLine, double[]> countExtractor) {
        final SimpleInterval interval = intervalExtractor.apply(v);
        final Target target = createRecordTarget(targets, interval);
        if (target == null) {
            if (ignoreMissingTargets) {
                return null;
            } else {
                throw formatException(String.format("unknown target with interval %s", interval));
            }
        } else {
            return new ReadCountRecord(target, countExtractor.apply(v));
        }
    }

    private ReadCountRecord extractRecordWithTargetNameAndIntervalExtractorsAndTargetCollection(
            final DataLine data, final Function<DataLine, SimpleInterval> intervalExtractor,
            final Function<DataLine, String> targetNameExtractor, final Function<DataLine, double[]> countExtractor) {
        final String name = targetNameExtractor.apply(data);
        final SimpleInterval interval = intervalExtractor.apply(data);
        final Target target = createRecordTarget(targets, interval);
        if (target == null) {
            return ignoreMissingTargets ? null : new ReadCountRecord(new Target(name, interval), countExtractor.apply(data));
        } else if (!target.getName().equals(name)){
            throw formatException(String.format("conflicting target resolution from the name (%s) and interval (%s) provided", name, interval));
        } else {
            return new ReadCountRecord(target, countExtractor.apply(data));
        }
    }

    @Override
    protected ReadCountRecord createRecord(final DataLine dataLine) {
        return recordExtractor.apply(dataLine);
    }

    private Function<DataLine, double[]> countExtractor(final TableColumnCollection columns) {
        final int[] countColumnIndexes = IntStream.range(0, columns.columnCount())
                .filter(i -> !TargetTableColumn.isStandardTargetColumnName(columns.nameAt(i))).toArray();
        return (v) -> {
            final double[] result = new double[countColumnIndexes.length];
            for (int i = 0; i < countColumnIndexes.length; i++) {
                result[i] = v.getDouble(countColumnIndexes[i]);
            }
            return result;
        };
    }

    public boolean hasTargetIntervals() {
        return targets != null
                || columns().containsAll(
                TargetTableColumn.CONTIG.toString(),
                TargetTableColumn.START.toString(),
                TargetTableColumn.END.toString());
    }
}
