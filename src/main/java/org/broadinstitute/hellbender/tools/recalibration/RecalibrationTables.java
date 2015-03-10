package org.broadinstitute.hellbender.tools.recalibration;

import org.broadinstitute.hellbender.tools.recalibration.covariates.Covariate;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * Utility class to facilitate base quality score recalibration.
 */

public final class RecalibrationTables implements Iterable<NestedIntegerArray<RecalDatum>>{
    private final List<NestedIntegerArray<RecalDatum>> tables;
    private final int qualDimension;
    private final int eventDimension = EventType.values().length;
    private final int numReadGroups;

    private final NestedIntegerArray<RecalDatum> readGroupTable;
    private final NestedIntegerArray<RecalDatum> qualityScoreTable;
    private final NestedIntegerArray<RecalDatum> contextTable;
    private final NestedIntegerArray<RecalDatum> cycleTable;

    public RecalibrationTables(final StandardCovariateList covariates) {
        this(covariates, covariates.getReadGroupCovariate().maximumKeyValue() + 1);
    }

    public RecalibrationTables(StandardCovariateList covariates, final int numReadGroups) {
        this.tables = new ArrayList<>(covariates.size());

        this.qualDimension = covariates.getQualityScoreCovariate().maximumKeyValue() + 1;
        this.numReadGroups = numReadGroups;

        this.readGroupTable = new NestedIntegerArray<>(numReadGroups, eventDimension);
        tables.add(readGroupTable);

        this.qualityScoreTable = makeQualityScoreTable();
        tables.add(qualityScoreTable);

        final Covariate contextCovariate = covariates.getContextCovariate();
        this.contextTable = new NestedIntegerArray<>(numReadGroups, qualDimension, contextCovariate.maximumKeyValue() + 1, eventDimension);
        tables.add(contextTable);

        final Covariate cycleCovariate = covariates.getCycleCovariate();
        this.cycleTable = new NestedIntegerArray<>(numReadGroups, qualDimension, cycleCovariate.maximumKeyValue() + 1, eventDimension);
        tables.add(cycleTable);
    }

    public boolean isReadGroupTable(NestedIntegerArray<RecalDatum> table) {
        return table.equals(getReadGroupTable());
    }

    public boolean isQualityScoreTable(NestedIntegerArray<RecalDatum> table) {
        return table.equals(getQualityScoreTable());
    }

    public boolean isContextTable(NestedIntegerArray<RecalDatum> table) {
        return table.equals(getContextTable());
    }

    public boolean isCycleTable(NestedIntegerArray<RecalDatum> table) {
        return table.equals(getCycleTable());
    }

    public List<NestedIntegerArray<RecalDatum>> getOptionalTables() {
        return Arrays.asList(contextTable, cycleTable);
    }

    public NestedIntegerArray<RecalDatum> getReadGroupTable() {
        return readGroupTable;
    }

    public NestedIntegerArray<RecalDatum> getQualityScoreTable() {
        return qualityScoreTable;
    }

    public NestedIntegerArray<RecalDatum> getContextTable() {
        return contextTable;
    }

    public NestedIntegerArray<RecalDatum> getCycleTable() {
        return cycleTable;
    }

    public NestedIntegerArray<RecalDatum> getTable(final int index) {
        return tables.get(index);
    }

    public int numTables() {
        return tables.size();
    }

    @Override
    public Iterator<NestedIntegerArray<RecalDatum>> iterator() {
        return tables.iterator();
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
