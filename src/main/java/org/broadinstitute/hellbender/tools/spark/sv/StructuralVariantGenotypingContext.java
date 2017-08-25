package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.List;

/**
 * Created by valentin on 8/25/17.
 */
public class StructuralVariantGenotypingContext implements Serializable {

    private static final long serialVersionUID = 1L;

    public final StructuralVariantContext variant;

    public final List<Tuple2<GATKRead, GATKRead>> readPairs;

    public final List<Haplotype> haplotypes;

    public final List<GenotypingContig> contigs;
}
