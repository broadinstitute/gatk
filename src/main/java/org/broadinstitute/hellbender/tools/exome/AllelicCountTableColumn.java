package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * Table columns of an allelic count tab-separated file at different verbosity levels.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum AllelicCountTableColumn {
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

    public static final TableColumnCollection BASIC_COLUMNS = new TableColumnCollection(
            CONTIG, POSITION, REF_COUNT, ALT_COUNT);

    public static final TableColumnCollection INTERMEDIATE_COLUMNS = new TableColumnCollection(
            CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE, READ_DEPTH);

    public static final TableColumnCollection FULL_COLUMNS = new TableColumnCollection((Object[]) values());

    /**
     * Get {@link AllelicCount} table columns at a given {@link AllelicCountTableVerbosity} level
     *
     * @param verbosity verbosity level
     * @return {@link TableColumnCollection} of column names
     * @throws UserException.BadInput
     */
    public static TableColumnCollection getColumns(final AllelicCountTableVerbosity verbosity) {
        switch (verbosity) {
            case BASIC:
                return BASIC_COLUMNS;
            case INTERMEDIATE:
                return INTERMEDIATE_COLUMNS;
            case FULL:
                return FULL_COLUMNS;
            default:
                throw new UserException.BadInput("The AllelicCount verbosity is invalid.");
        }
    }
}
