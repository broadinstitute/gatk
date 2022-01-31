package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Set;

final class VariantDatum implements Locatable {
    public double[] annotations;
    public Set<String> labels;
    public SimpleInterval loc;
    public Allele referenceAllele;
    public Allele alternateAllele;
    public double score;

    @Override
    public String getContig() {
        return loc.getContig();
    }

    @Override
    public int getStart() {
        return loc.getStart();
    }

    @Override
    public int getEnd() {
        return loc.getEnd();
    }
}