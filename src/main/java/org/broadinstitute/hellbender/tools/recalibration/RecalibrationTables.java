package org.broadinstitute.hellbender.tools.recalibration;

import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

import java.util.ArrayList;

/**
 * Utility class to facilitate base quality score recalibration.
 */

public final class RecalibrationTables {
    public enum TableType {
        READ_GROUP_TABLE,
        QUALITY_SCORE_TABLE,
        OPTIONAL_COVARIATE_TABLES_START;
    }

    private final ArrayList<NestedIntegerArray<RecalDatum>> tables;
    private final int qualDimension;
    private final int eventDimension = EventType.values().length;
    private final int numReadGroups;

    public RecalibrationTables(final Covariate[] covariates) {
        this(covariates, covariates[TableType.READ_GROUP_TABLE.ordinal()].maximumKeyValue() + 1);
    }

    public RecalibrationTables(final Covariate[] covariates, final int numReadGroups) {
        tables = new ArrayList<>(covariates.length);
        for ( int i = 0; i < covariates.length; i++ )
            tables.add(i, null); // initialize so we can set below

        qualDimension = covariates[TableType.QUALITY_SCORE_TABLE.ordinal()].maximumKeyValue() + 1;
        this.numReadGroups = numReadGroups;

        tables.set(TableType.READ_GROUP_TABLE.ordinal(),
                 new NestedIntegerArray<>(numReadGroups, eventDimension));

        tables.set(TableType.QUALITY_SCORE_TABLE.ordinal(), makeQualityScoreTable());

        for (int i = TableType.OPTIONAL_COVARIATE_TABLES_START.ordinal(); i < covariates.length; i++)
            tables.set(i, new NestedIntegerArray<>(numReadGroups, qualDimension, covariates[i].maximumKeyValue()+1, eventDimension));
    }

    public NestedIntegerArray<RecalDatum> getReadGroupTable() {
        return getTable(TableType.READ_GROUP_TABLE.ordinal());
    }

    public NestedIntegerArray<RecalDatum> getQualityScoreTable() {
        return getTable(TableType.QUALITY_SCORE_TABLE.ordinal());
    }

    public NestedIntegerArray<RecalDatum> getTable(final int index) {
        return tables.get(index);
    }

    public int numTables() {
        return tables.size();
    }

    /**
     * @return true if all the tables contain no RecalDatums
     */
    public boolean isEmpty() {
        for( final NestedIntegerArray<RecalDatum> table : tables ) {
            if( !table.getAllValues().isEmpty() ) { return false; }
        }
        return true;
    }

    /**
     * Allocate a new quality score table, based on requested parameters
     * in this set of tables, without any data in it.  The return result
     * of this table is suitable for acting as a thread-local cache
     * for quality score values
     * @return a newly allocated, empty read group x quality score table
     */
    public NestedIntegerArray<RecalDatum> makeQualityScoreTable() {
        return new NestedIntegerArray<>(numReadGroups, qualDimension, eventDimension);
    }

    /**
     * Merge all of the tables from toMerge into into this set of tables
     */
    public void combine(final RecalibrationTables toMerge) {
        if ( numTables() != toMerge.numTables() )
            throw new IllegalArgumentException("Attempting to merge RecalibrationTables with different sizes");

        for ( int i = 0; i < numTables(); i++ ) {
            final NestedIntegerArray<RecalDatum> myTable = this.getTable(i);
            final NestedIntegerArray<RecalDatum> otherTable = toMerge.getTable(i);
            RecalUtils.combineTables(myTable, otherTable);
        }
    }
}
