package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Collection of integer copy-number states.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class IntegerCopyNumberStateCollection {

    private final List<IntegerCopyNumberState> copyNumberStates;
    private final TableColumnCollection columnCollection;

    public IntegerCopyNumberStateCollection(final List<String> copyNumberStatesColumns) {
        Utils.nonEmpty(copyNumberStatesColumns);
        this.columnCollection = new TableColumnCollection(copyNumberStatesColumns);
        this.copyNumberStates = new ArrayList<>();
        copyNumberStatesColumns
                .forEach(copyNumberString -> copyNumberStates.add(parseIntegerCopyNumber(copyNumberString)));
    }

    /**
     * Get list of all copy number states in this collection
     */
    public List<IntegerCopyNumberState> getCopyNumberStates() {
        return copyNumberStates;
    }

    /**
     * Get copy number state given an its index
     */
    public IntegerCopyNumberState get(final int index) {
        return copyNumberStates.get(index);
    }

    /**
     * Get the number of states in this collection
     */
    public int size() {
        return columnCollection.columnCount();
    }

    /**
     * Get the correspodning {@link TableColumnCollection} for this copy number state collection
     */
    public TableColumnCollection getTableColumnCollection() {
        return columnCollection;
    }

    private IntegerCopyNumberState parseIntegerCopyNumber(final String copyNumberStateString) {
        if (!copyNumberStateString.startsWith(GermlineCNVNamingConstants.COPY_NUMBER_TABLE_COLUMN_PREFIX)) {
            throw new UserException.BadInput(String.format("The column names of the copy number posterior file must " +
                    "start with %s", GermlineCNVNamingConstants.COPY_NUMBER_TABLE_COLUMN_PREFIX));
        }
        final String integerCopyNumberStateString = copyNumberStateString.substring(GermlineCNVNamingConstants.COPY_NUMBER_TABLE_COLUMN_PREFIX.length());
        try {
            final int copyNumber = Integer.parseInt(integerCopyNumberStateString);
            return new IntegerCopyNumberState(copyNumber);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput(e.getMessage()); //TODO write a better exception message
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final IntegerCopyNumberStateCollection that = (IntegerCopyNumberStateCollection) o;

        return copyNumberStates != null ? copyNumberStates.equals(that.copyNumberStates) : that.copyNumberStates == null;
    }

    @Override
    public int hashCode() {
        return copyNumberStates != null ? copyNumberStates.hashCode() : 0;
    }

    @Override
    public String toString() {
        return copyNumberStates.stream()
                .map(IntegerCopyNumberState::toString)
                .collect(Collectors.joining(", ", "[", "]"));
    }
}
