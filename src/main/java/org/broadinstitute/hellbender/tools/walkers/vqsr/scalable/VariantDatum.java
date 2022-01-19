package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Comparator;

// TODO refactor so VariantDatums can essentially be final; they should also extend some Locatable
final class VariantDatum {
    public double[] annotations;
    public boolean[] isNull;
    public boolean atTruthSite;
    public boolean atTrainingSite;
    public boolean isTransition;
    public SimpleInterval loc;
    public Allele referenceAllele;
    public Allele alternateAllele;
    public double score;
}