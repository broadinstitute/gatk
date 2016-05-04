package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Table columns of an allelic count tab-separated file at different verbosity levels.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public enum AllelicCountTableColumns {
    CONTIG,
    POSITION,
    REF_COUNT,
    ALT_COUNT,
    REF_NUCLEOTIDE,
    ALT_NUCLEOTIDE,
    READ_DEPTH,
    HET_LOG_ODDS;

    public enum AllelicCountTableVerbosity {
        BASIC, INTERMEDIATE, FULL
    }

    public static final List<String> BASIC_COLUMN_NAME_ARRAY = createUnmodifiableList(
            CONTIG, POSITION, REF_COUNT, ALT_COUNT);

    public static final List<String> INTERMEDIATE_COLUMN_NAME_ARRAY = createUnmodifiableList(
            CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE, READ_DEPTH);

    public static final List<String> FULL_COLUMN_NAME_ARRAY = createUnmodifiableList(
            CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE, READ_DEPTH, HET_LOG_ODDS);

    /**
     * Get {@link AllelicCount} table columns at a given {@link AllelicCountTableVerbosity} level
     *
     * @param verbosity verbosity level
     * @return list of column names
     * @throws UserException.BadInput
     */
    public static List<String> getColumns(final AllelicCountTableVerbosity verbosity) {
        switch (verbosity) {
            case BASIC:
                return BASIC_COLUMN_NAME_ARRAY;
            case INTERMEDIATE:
                return INTERMEDIATE_COLUMN_NAME_ARRAY;
            case FULL:
                return FULL_COLUMN_NAME_ARRAY;
            default:
                throw new UserException.BadInput("The AllelicCount verbosity is invalid.");
        }
    }

    private static List<String> createUnmodifiableList(final AllelicCountTableColumns ... columns) {
        return Collections.unmodifiableList(Arrays.asList(columns).stream()
                .map(AllelicCountTableColumns::name).collect(Collectors.toList()));
    }
}
