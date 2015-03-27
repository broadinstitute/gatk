package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.Comparator;
import java.util.List;

/*
 * Represents a data item for VQSR (a site).
 * Package private because it's not usable outside of VQSR.
 */
final class VariantDatum {

    public double[] annotations; //values of the annotations
    public boolean[] isNull;     //is any of these values "null" (empty) - used for marginalization for plotting
    public boolean isKnown;      // is this a known site
    public double lod;            //lod score from the model
    public boolean atTruthSite;
    public boolean atTrainingSite;
    public boolean atAntiTrainingSite;
    public boolean isTransition;
    public boolean isSNP;
    public boolean failingSTDThreshold;
    public double prior;                 //per-site prior
    public GenomeLoc loc;
    public int worstAnnotation;
    public boolean isAggregate; // this datum was provided to aid in modeling but isn't part of the input callset

    public static final Comparator<VariantDatum> VariantDatumLODComparator = (datum1, datum2) -> Double.compare(datum1.lod, datum2.lod);

    public static int countCallsAtTruth(final List<VariantDatum> data, double minLOD ) {
        return (int)data.stream().filter(d -> (d.atTruthSite && d.lod >= minLOD)).count(); //XXX cast to int for compatibility
    }
}