package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Set;

// TODO refactor so VariantDatums can essentially be final; they should also extend some Locatable
final class VariantDatum implements Locatable {
    public double[] annotations; // TODO consider HashMap with keys
    public Set<String> labels;     // TODO consider HashMap with keys
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