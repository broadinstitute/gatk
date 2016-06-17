package org.broadinstitute.hellbender.utils.tsv;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.HashSet;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Common constants for table readers and writers.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TableUtils {

    /**
     * Column separator {@value #COLUMN_SEPARATOR_STRING}.
     */
    public static final char COLUMN_SEPARATOR = '\t';

    /**
     * Column separator as an string.
     */
    public static final String COLUMN_SEPARATOR_STRING = String.valueOf(COLUMN_SEPARATOR);

    /**
     * Comment line prefix string {@value}.
     * <p>
     * Lines that start with this prefix (spaces are not ignored), will be considered comment
     * lines (neither a header line nor data line).
     * </p>
     */
    public static final String COMMENT_PREFIX = "#";

    /**
     * Quote character {@value #QUOTE_STRING}.
     * <p>
     * Character used to quote table values that contain special characters.
     * </p>
     */
    public static final char QUOTE_CHARACTER = '\"';

    /**
     * Quote character as a string.
     */
    public static final String QUOTE_STRING = String.valueOf(QUOTE_CHARACTER);

    /**
     * Escape character {@value #ESCAPE_STRING}.
     * <p>
     * Within {@value #QUOTE_STRING quotes},
     * the user must prepend this character when including {@value #QUOTE_STRING quotes}
     * or the {@value #ESCAPE_STRING escape}
     * character itself as part of the string.
     * </p>
     */
    public static final char ESCAPE_CHARACTER = '\\';

    /**
     * Escape character as a string.
     */
    public static final String ESCAPE_STRING = String.valueOf(ESCAPE_CHARACTER);

    /**
     * Creates a new table reader given an record extractor factory based from the columns found in the input.
     * <p>
     * The record extractor factory would take on to arguments where the first one is a {@link TableColumnCollection}
     * with the input columns and the second a exception factory to be used in order to compose the exception
     * to throw when there is a formatting error.
     * </p>
     * <p>
     * This factory must return a function that maps {@link DataLine} instances into record instances
     * (generic type {@link R}).
     * </p>
     * <p>
     * The record factory must never return {@code null} and it must use the format exception factory provided
     * to extractor factory in order to compose the exception to throw in case that there is a formatting
     * error.
     * </p>
     *
     * @param file                   the input file
     * @param recordExtractorFactory the record extractor function factory.
     * @param <R>                    the end record type.
     * @return never {@code null}.
     * @throws IOException              if any took place while instantiating the reader.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     */
    public static <R> TableReader<R> reader(final File file,
                                            final BiFunction<TableColumnCollection, Function<String, RuntimeException>, Function<DataLine, R>> recordExtractorFactory)
            throws IOException {
        Utils.nonNull(recordExtractorFactory,"the record extractor factory cannot be null");
        return new TableReader<R>(file) {
            private Function<DataLine,R> recordExtractor;

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                recordExtractor = recordExtractorFactory.apply(columns,this::formatException);
                if (recordExtractor == null) {
                    throw new IllegalStateException("the record extractor function cannot be null");
                }
            }

            @Override
            protected R createRecord(DataLine dataLine) {
                return recordExtractor.apply(dataLine);
            }
        };
    }

    /**
     * Creates a new table reader given an record extractor factory based from the columns found in the input.
     * <p>
     * The record extractor factory would take on to arguments where the first one is a {@link TableColumnCollection}
     * with the input columns and the second a exception factory to be used in order to compose the exception
     * to throw when there is a formatting error.
     * </p>
     * <p>
     * This factory must return a function that maps {@link DataLine} instances into record instances
     * (generic type {@link R}).
     * </p>
     * <p>
     * The record factory must never return {@code null} and it must use the format exception factory provided
     * to extractor factory in order to compose the exception to throw in case that there is a formatting
     * error.
     * </p>
     *
     * @param reader                 the input reader.
     * @param recordExtractorFactory the record extractor function factory.
     * @param <R>                    the end record type.
     * @return never {@code null}.
     * @throws IOException if any took place while instantiating the reader.
     */
    public static <R> TableReader<R> reader(final Reader reader,
                                            final BiFunction<TableColumnCollection, Function<String, RuntimeException>, Function<DataLine, R>> recordExtractorFactory)
            throws IOException {
        Utils.nonNull(recordExtractorFactory,"the record extractor factory cannot be null");
        return new TableReader<R>(reader) {
            private Function<DataLine,R> recordExtractor;

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                recordExtractor = recordExtractorFactory.apply(columns,this::formatException);
                if (recordExtractor == null) {
                    throw new IllegalStateException("the record extractor function cannot be null");
                }
            }

            @Override
            protected R createRecord(DataLine dataLine) {
                return recordExtractor.apply(dataLine);
            }
        };
    }

    /**
     * Creates a new table reader given an record extractor factory based from the columns found in the input.
     * <p>
     * The record extractor factory would take on to arguments where the first one is a {@link TableColumnCollection}
     * with the input columns and the second a exception factory to be used in order to compose the exception
     * to throw when there is a formatting error.
     * </p>
     * <p>
     * This factory must return a function that maps {@link DataLine} instances into record instances
     * (generic type {@link R}).
     * </p>
     * <p>
     * The record factory must never return {@code null} and it must use the format exception factory provided
     * to extractor factory in order to compose the exception to throw in case that there is a formatting
     * error.
     * </p>
     *
     * @param sourceName             the source name, {@code null} indicates that the source is anonymous.
     * @param reader                 the input reader.
     * @param recordExtractorFactory the record extractor function factory.
     * @param <R>                    the end record type.
     * @return never {@code null}.
     * @throws IOException if any took place while instantiating the reader.
     */
    public static <R> TableReader<R> reader(final String sourceName,
                                            final Reader reader,
                                            final BiFunction<TableColumnCollection, Function<String, RuntimeException>, Function<DataLine, R>> recordExtractorFactory)
            throws IOException {
        Utils.nonNull(recordExtractorFactory,"the record extractor factory cannot be null");
        return new TableReader<R>(sourceName, reader) {
            private Function<DataLine,R> recordExtractor;

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                recordExtractor = recordExtractorFactory.apply(columns,this::formatException);
                if (recordExtractor == null) {
                    throw new IllegalStateException("the record extractor function cannot be null");
                }
            }

            @Override
            protected R createRecord(DataLine dataLine) {
                return recordExtractor.apply(dataLine);
            }
        };
    }

    /**
     * Creates a new table writer given the destination file, columns and the data-line composer.
     * @param file the destination file.
     * @param columns the output columns.
     * @param dataLineComposer the data-line composer given the record object.
     * @param <R> the record type.
     * @return never {@code null}.
     * @throws IllegalArgumentException if any, {@code file}, {@code columns} or {@code dataLineComposer}, is {@code null}.
     * @throws IOException if any was thrown when instantiating the writer.
     */
    public static <R> TableWriter<R> writer(final File file, final TableColumnCollection columns, final BiConsumer<R, DataLine> dataLineComposer) throws IOException {
        return new DataLineComposerBasedTableWriter<>(file, columns, dataLineComposer);
    }

    /**
     * Creates a new table writer given the destination writer, columns and the data-line composer.
     * @param writer the destination writer.
     * @param columns the output columns.
     * @param dataLineComposer the data-line composer given the record object.
     * @param <R> the record type.
     * @return never {@code null}.
     * @throws IllegalArgumentException if any, {@code writer}, {@code columns} or {@code dataLineComposer}, is {@code null}.
     * @throws IOException if any was thrown when instantiating the writer.
     */
    public static <R> TableWriter<R> writer(final Writer writer, final TableColumnCollection columns, final BiConsumer<R, DataLine> dataLineComposer) throws IOException {
        return new DataLineComposerBasedTableWriter<>(writer, columns, dataLineComposer);
    }

    /**
     * Checks if all mandatory columns are present in a {@link TableColumnCollection}.
     * @param columns                   the TableColumnCollection of columns to check
     * @param mandatoryColumns          the TableColumnCollection of mandatory columns
     * @param formatExceptionFactory    the format exception function factory
     * @throws UserException.BadInput   if any mandatory columns are missing
     */
    public static void checkMandatoryColumns(final TableColumnCollection columns, final TableColumnCollection mandatoryColumns,
                                             final Function<String, RuntimeException> formatExceptionFactory) {
        if (!columns.containsAll(mandatoryColumns.names())) {
            final Set<String> missingColumns = Sets.difference(new HashSet<>(columns.names()), new HashSet<>(mandatoryColumns.names()));
            throw formatExceptionFactory.apply("Bad header in file.  Not all mandatory columns are present.  Missing: " + StringUtils.join(missingColumns, ", "));
        }
    }

    private static final class DataLineComposerBasedTableWriter<R> extends TableWriter<R> {
        private final BiConsumer<R, DataLine> dataLineComposer;

        private DataLineComposerBasedTableWriter(final File file, final TableColumnCollection columns, final BiConsumer<R, DataLine> dataLineComposer) throws IOException {
            super(file, columns);
            this.dataLineComposer = Utils.nonNull(dataLineComposer, "the data-line composer cannot be null");
        }

        private DataLineComposerBasedTableWriter(final Writer writer, final TableColumnCollection columns, final BiConsumer<R, DataLine> dataLineComposer) throws IOException {
            super(writer, columns);
            this.dataLineComposer = Utils.nonNull(dataLineComposer, "the data-line composer cannot be null");
        }

        @Override
        protected void composeLine(R record, DataLine dataLine) {
            dataLineComposer.accept(record, dataLine);
        }
    }

    /**
     * Declared to make instantiation impossible.
     */
    private TableUtils() {
        throw new UnsupportedOperationException();
    }
}
