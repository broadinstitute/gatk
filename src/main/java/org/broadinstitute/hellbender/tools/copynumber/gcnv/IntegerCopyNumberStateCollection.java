package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Collection of copy number states considered in the germline CNV calling pipeline
 */
public class IntegerCopyNumberStateCollection {

    private final List<IntegerCopyNumberState> copyNumberStates;
    private final TableColumnCollection columnCollection;
    private final List<Allele> alleles;

    /**
     * The state with this copy number represents a reference allele
     */
    private static final int REFERENCE_COPY_NUMBER_VALUE = 2;

    public IntegerCopyNumberStateCollection(final List<String> copyNumberStatesColumns) {
        Utils.nonEmpty(copyNumberStatesColumns);
        this.columnCollection = new TableColumnCollection(copyNumberStatesColumns);
        this.copyNumberStates = new ArrayList<>();
        copyNumberStatesColumns
                .forEach(copyNumberString -> copyNumberStates.add(parseIntegerCopyNumber(copyNumberString)));
        this.alleles = copyNumberStates.stream()
                .map(state -> state.toAllele(REFERENCE_COPY_NUMBER_VALUE))
                .collect(Collectors.toList());
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

    /**
     * Get corresponding alleles for this set of copy number states
     */
    public List<Allele> getAlleles() {
        return alleles;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || o.getClass() != this.getClass()) {
            return false;
        }

        final IntegerCopyNumberStateCollection that = (IntegerCopyNumberStateCollection) o;

        return that.copyNumberStates.equals(this.copyNumberStates);
    }

    @Override
    public int hashCode() {
        return copyNumberStates.stream().map(IntegerCopyNumberState::hashCode).reduce(1, (a, b) -> a * 31 + b);
    }
}
