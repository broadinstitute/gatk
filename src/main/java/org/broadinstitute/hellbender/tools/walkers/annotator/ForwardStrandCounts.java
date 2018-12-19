package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticGenotypingEngine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.OptionalInt;

public class ForwardStrandCounts extends PerAlleleAnnotation implements StandardMutectAnnotation {
    @Override
    protected int aggregate(final List<Integer> values) {
        return values.size();
    }

    @Override
    protected String getVcfKey() { return GATKVCFConstants.FORWARD_STRAND_COUNT_KEY; }

    @Override
    protected String getDescription() { return "Forward strand read counts for each allele"; }

    @Override
    protected boolean includeRefAllele() { return true; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        return OptionalInt.of(read.isReverseStrand() ? 0 : 1);
    }
}
