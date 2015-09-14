package org.broadinstitute.hellbender.tools.walkers.bqsr;

import org.broadinstitute.hellbender.tools.recalibration.ReadCovariates;
import org.broadinstitute.hellbender.tools.recalibration.RecalDatum;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.EventType;

import java.io.Serializable;

public final class RecalibrationEngine implements Serializable {
    private static final long serialVersionUID = 1L;

    protected final StandardCovariateList covariates;

    /**
     * Has finalizeData() been called?
     */
    private boolean finalized = false;

    private final RecalibrationTables tables;

    /**
     * Initialize the recalibration engine
     *
     * Called once before any calls to updateDataForRead are made.  The engine should prepare itself
     * to handle any number of updateDataForRead calls containing ReadRecalibrationInfo containing
     * keys for each of the covariates provided.
     *
     * The engine should collect match and mismatch data into the recalibrationTables data.
     *
     * @param covariates an array of the covariates we'll be using in this engine, order matters
     * @param numReadGroups the number of read groups we should use for the recalibration tables
     */
    public RecalibrationEngine(final StandardCovariateList covariates, final int numReadGroups) {
        if ( covariates == null ) {
            throw new IllegalArgumentException("Covariates cannot be null");
        }
        if ( numReadGroups < 1 ) {
            throw new IllegalArgumentException("numReadGroups must be >= 1 but got " + numReadGroups);
        }

        this.covariates = covariates;
        this.tables = new RecalibrationTables(covariates, numReadGroups);
    }

    /**
     * Update the recalibration statistics using the information in recalInfo
     * @param recalInfo data structure holding information about the recalibration values for a single read
     */
    public void updateDataForRead( final ReadRecalibrationInfo recalInfo ) {
        if ( finalized ) {
            throw new IllegalStateException("FinalizeData() has already been called");
        }

        final GATKRead read = recalInfo.getRead();
        final ReadCovariates readCovariates = recalInfo.getCovariatesValues();
        final NestedIntegerArray<RecalDatum> qualityScoreTable = tables.getQualityScoreTable();

        for( int offset = 0; offset < read.getBases().length; offset++ ) {
            if( ! recalInfo.skip(offset) ) {

                for (final EventType eventType : EventType.values()) {
                    final int[] keys = readCovariates.getKeySet(offset, eventType);
                    final int eventIndex = eventType.ordinal();
                    final byte qual = recalInfo.getQual(eventType, offset);
                    final double isError = recalInfo.getErrorFraction(eventType, offset);

                    RecalUtils.incrementDatumOrPutIfNecessary(qualityScoreTable, qual, isError, keys[0], keys[1], eventIndex);

                    for (int i = 2; i < covariates.size(); i++) { //XXX the 2 is hard-wired here as the number of special covariates
                        if (keys[i] < 0) {
                            continue;
                        }

                        RecalUtils.incrementDatumOrPutIfNecessary(tables.getTable(i), qual, isError, keys[0], keys[1], keys[i], eventIndex);
                    }
                }
            }
        }
    }


    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to updateDataForRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public void finalizeData() {
        if ( finalized ) {
            throw new IllegalStateException("FinalizeData() has already been called");
        }

        final NestedIntegerArray<RecalDatum> byReadGroupTable = tables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> byQualTable = tables.getQualityScoreTable();

        // iterate over all values in the qual table
        for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[2];
            final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
            final RecalDatum qualDatum = leaf.value;

            if ( rgDatum == null ) {
                // create a copy of qualDatum, and initialize byReadGroup table with it
                byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            } else {
                // combine the qual datum with the existing datum in the byReadGroup table
                rgDatum.combine(qualDatum);
            }
        }

        finalized = true;
    }

    /**
     * Finalize, if appropriate, all derived data in recalibrationTables.
     *
     * Called once after all calls to updateDataForRead have been issued.
     *
     * Assumes that all of the principal tables (by quality score) have been completely updated,
     * and walks over this data to create summary data tables like by read group table.
     */
    public static void finalizeRecalibrationTables(final RecalibrationTables tables) {
        Utils.nonNull(tables);
        final NestedIntegerArray<RecalDatum> byReadGroupTable = tables.getReadGroupTable();
        final NestedIntegerArray<RecalDatum> byQualTable = tables.getQualityScoreTable();

        // iterate over all values in the qual table
        for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {
            final int rgKey = leaf.keys[0];
            final int eventIndex = leaf.keys[2];
            final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
            final RecalDatum qualDatum = leaf.value;

            if ( rgDatum == null ) {
                // create a copy of qualDatum, and initialize byReadGroup table with it
                byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            } else {
                // combine the qual datum with the existing datum in the byReadGroup table
                rgDatum.combine(qualDatum);
            }
        }
    }

    /**
     * Get a possibly not-final recalibration table, to deal with distributed execution.
     */
    public RecalibrationTables getRecalibrationTables() {
        return tables;
    }

   /**
     * Get the final recalibration tables, after finalizeData() has been called
     *
     * This returns the finalized recalibration table collected by this engine.
     *
     * It is an error to call this function before finalizeData has been called
     *
     * @return the finalized recalibration table collected by this engine
     */
    public RecalibrationTables getFinalRecalibrationTables() {
        if ( ! finalized ) {
            throw new IllegalStateException("Cannot get final recalibration tables until finalizeData() has been called");
        }
        return tables;
    }
}
