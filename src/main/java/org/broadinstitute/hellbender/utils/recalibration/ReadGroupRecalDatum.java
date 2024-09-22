package org.broadinstitute.hellbender.utils.recalibration;

import ngs.ReadGroup;

//public class ReadGroupRecalDatum extends RecalDatum {
//    /**
//     * estimated reported quality score based on combined data's individual q-reporteds and number of observations
//     */ // tsato: Explain
//    private double estimatedQReported; // tsato: estimated by Illumina? Should be renamed to OriginalQuality
//    // tsato: This column only occurs in ReadGroupTable ..... so consider inheriting RecalDatum.
//
//
//    public ReadGroupRecalDatum(final ReadGroup readGroup) {
//        super(readGroup);
//    }
//
//    public RecalDatum(final long _numObservations, final double _numMismatches, final byte reportedQuality) {
//        if ( _numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
//        if ( _numMismatches < 0.0 ) throw new IllegalArgumentException("numMismatches < 0");
//        if ( reportedQuality < 0 ) throw new IllegalArgumentException("reportedQuality < 0");
//
//        numObservations = _numObservations;
//        numMismatches = (_numMismatches*MULTIPLIER);
//        estimatedQReported = reportedQuality; // tsato: estimated is very confusing here
//        empiricalQuality = UNINITIALIZED;
//    }
//
//    /**
//     * Copy copy into this recal datum, overwriting all of this objects data
//     * @param copy  RecalDatum to copy
//     */
//    public RecalDatum(final RecalDatum copy) {
//        this.numObservations = copy.numObservations;
//        this.numMismatches = copy.numMismatches;
//        this.estimatedQReported = copy.estimatedQReported;
//        this.empiricalQuality = copy.empiricalQuality;
//    }
//}
