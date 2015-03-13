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

    /**
     * Computes and sets the index of the worst performing annotation.
     * That's the annotation with the lowest ratio of likelihoods from the good and bad models.
     * @return the index of the worst annotation or -1 if probabilities for all dimensions are null.
     */
    public void setWorstPerformingAnnotation( final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel ) {
        int worstAnnotation = -1;
        double minProb = Double.MAX_VALUE;
        for( int i = 0; i < this.annotations.length; i++ ) {
            final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(this, i);
            final Double badProbLog10 = badModel.evaluateDatumInOneDimension(this, i);
            if( goodProbLog10 != null && badProbLog10 != null ) {
                final double prob = goodProbLog10 - badProbLog10;
                if (prob < minProb) {
                    minProb = prob;
                    worstAnnotation = i;
                }
            }
        }
        this.worstAnnotation = worstAnnotation;
    }
}
