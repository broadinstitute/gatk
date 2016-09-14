package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Comparator;
import java.util.List;

/*
 * Represents a data item for VQSR (a site).
 * Package private because it's not usable outside of VQSR.
 */
final class VariantDatum {

    public double[] annotations;
    public boolean[] isNull;
    public boolean isKnown;
    public double lod;
    public boolean atTruthSite;
    public boolean atTrainingSite;
    public boolean atAntiTrainingSite;
    public boolean isTransition;
    public boolean isSNP;
    public boolean failingSTDThreshold;
    public double originalQual;
    public double prior;
    public SimpleInterval loc;
    public int worstAnnotation;
    public double worstValue;
    public MultivariateGaussian assignment; // used in K-means implementation
    public boolean isAggregate; // this datum was provided to aid in modeling but isn't part of the input callset
    public Allele referenceAllele;
    public Allele alternateAllele;

    public static final Comparator<VariantDatum> VariantDatumLODComparator = (datum1, datum2) -> Double.compare(datum1.lod, datum2.lod);

    public static int countCallsAtTruth(final List<VariantDatum> data, double minLOD ) {
        return (int)data.stream().filter(d -> (d.atTruthSite && d.lod >= minLOD)).count(); //XXX cast to int for compatibility
    }

    /**
     * Return a comparator for VariantDatums, given a sequence Dictionary.
     * @param seqDictionary
     * @return a lambda closure comparator that uses the provided sequence dictionary
     */
    public static Comparator<VariantDatum> getComparator(final SAMSequenceDictionary seqDictionary) {
        return (vd1, vd2) -> IntervalUtils.compareLocatables(vd1.loc, vd2.loc, seqDictionary);
    }

}