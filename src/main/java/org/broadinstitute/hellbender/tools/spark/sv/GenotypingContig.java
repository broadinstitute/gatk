package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.Map;

/**
 * Created by valentin on 8/25/17.
 */
public class GenotypingContig {

    public final Map<Haplotype, AlignedContig> haplotypeAlignments;

    public final Map<Haplotype, AlignedContigScore> haplotypeScores;

    public final Haplotype call;

    public final int qual;

}
