package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Comparator;
import java.util.List;

final class VariantDatum {

    public double[] annotations;
    public boolean[] isNull;
    public boolean atTruthSite;
    public boolean atTrainingSite;
    public boolean isSNP;
    public SimpleInterval loc;
    public Allele referenceAllele;
    public Allele alternateAllele;

    public static int countCallsAtTruth(final List<VariantDatum> data) {
        return (int) data.stream().filter(d -> d.atTruthSite).count(); //XXX cast to int for compatibility
    }

    /**
     * Return a comparator for VariantDatums, given a sequence Dictionary.
     * @return a lambda closure comparator that uses the provided sequence dictionary
     */
    public static Comparator<VariantDatum> getComparator(final SAMSequenceDictionary seqDictionary) {
        return (vd1, vd2) -> IntervalUtils.compareLocatables(vd1.loc, vd2.loc, seqDictionary);
    }

}